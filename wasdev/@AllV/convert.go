package main

import (
	"fmt"
	"github.com/wsxiaoys/terminal/color"
	"io/ioutil"
	"log"
	"os"
	"regexp"
	"sort"
	"strings"
)

const (
	FUNCTION_NAME_REGEX       = `^[ \t]*function[A-Z,a-z,\,,\[,\],=, ,\t]* ([a-z,0-9,A-Z]+)[ ]*\(`
	FUNCTION_DEFINITION_REGEX = `^[ \t]*function(.*)`
)

func reversSortStringSlice(s []string) []string {
	sort.Strings(s)
	var ss sort.StringSlice
	ss = s
	sort.Sort(sort.Reverse(ss))
	return ss
}

func getFuncNames(filecontent []string) ([]string, []string) {
	r := make([]string, 0)
	r2 := make([]string, 0)
	var funcCounter int64
	functionRx := regexp.MustCompile(FUNCTION_NAME_REGEX)
	functionDefRx := regexp.MustCompile(FUNCTION_DEFINITION_REGEX)
	for lineno, line := range filecontent {
		if functionRx.MatchString(line) {
			funcnameB := functionRx.FindSubmatch([]byte(line))
			funcname := string(funcnameB[1])
			funcCounter = funcCounter + 1
			fmt.Printf("%v %d ", funcCounter, lineno)
			color.Print("@g", funcname)
			color.Println("@y|", line)
			r = append(r, funcname)

			if functionDefRx.MatchString(line) {
				funcdefB := functionDefRx.FindSubmatch([]byte(line))
				funcdef := string(funcdefB[1])
				r2 = append(r2, funcdef)
			} else {
				log.Fatal("Can detect the funciton name but not the function definition! That is not acceptable")
			}
		}
	}
	return r, r2
}

func renameFuncs(classname string, funclines []string, funcnames []string) []string {
	r := funclines
	declaration := funclines[0]
	for _, funcname := range funcnames {
		//Test cases:
		// funcname -> AliFunc
		/*
			if length(a)
						AliFuncs
						AliFunc ()
						sonthing
						AliFunc()
						somthindsfjf
						AliFunc(124314,34124,234,213)
						AliFunc (124314,34124,234,213)
						a = AliFunc()
						[a, b]=AliFunc()
						a = PlotRateSequence(aExpts,'color',colors{j},scaling,'offset',(j-1)*2,'bytime','callback',{@AliFunc, cellids(j)});
						text(max(a.AliFunc)*1.01,a.meanrates(AliFunc),AliFunc(cellids(j)),'color',colors{j},'horizontalalignment','left');
					end

		*/
		rx := ""
		rx = rx + `(?:[ \[\t\,\(=@\s\n\r]+|^)(`
		rx = rx + funcname
		rx = rx + `)[ ,\t,\(,\),\s]+`
		funcRx := regexp.MustCompile(rx)
		funcnameRx := regexp.MustCompile("(" + funcname + ")")
		for i, _ := range r[1:] { // dont change the function declarations
			r[i] = funcRx.ReplaceAllStringFunc(r[i], func(s string) string {
				return funcnameRx.ReplaceAllString(s, classname+"."+funcname)
			})
		}

	}
	r[0] = declaration // this embarrassing
	return r
}

func main() {
	//originalFile := "./OriginalPlotClusters.m"
	//className := "PC"

	originalFile := "./OriginalAllVPCs.m"
	className := "AllV"

	b, err := ioutil.ReadFile(originalFile)
	if err != nil {
		log.Fatal("Cant read the original,", err)
	}

	lines := strings.Split(string(b), "\n")

	funcnames, funcdefs := getFuncNames(lines)
	//revese the funcnames, to reduce the change of collision, like a funcname being embded inside another function name
	funcnames = reversSortStringSlice(funcnames)

	//TODO: in all the function declaration lines, replace `=([A-z])` with '= \1'
	functionRx := regexp.MustCompile(FUNCTION_NAME_REGEX)

	var funcCounter int64
	//insidefunc := false
	//TODO: the assuption is the original file STARTS with a function declaration line, no comments of white space before it
	funclines := make([]string, 0)
	funcname := ""
	lineno := 0
	line := ""
	for lineno, line = range lines {
		if functionRx.MatchString(line) {
			if len(funclines) != 0 {
				patchedfunclines := renameFuncs(className, funclines, funcnames)
				color.Println("@rWriting: ", funcname)
				err = ioutil.WriteFile(funcname+".m", []byte(strings.Join(patchedfunclines, "")), os.ModePerm)
				if err != nil {
					log.Fatal("Error writing the file", err)
				}
				funclines = make([]string, 0)
			}
			funcnameB := functionRx.FindSubmatch([]byte(line))
			funcname = string(funcnameB[1])
			funcCounter = funcCounter + 1
			fmt.Printf("%v %d ", funcCounter, lineno)
			color.Print("@g", funcname)
			color.Println("@y|", line)
		}
		funclines = append(funclines, line)
	}
	//write the last function
	if len(funclines) != 0 {
		patchedfunclines := renameFuncs(className, funclines, funcnames)
		color.Println("@rWriting: ", funcname)
		err = ioutil.WriteFile(funcname+".m", []byte(strings.Join(patchedfunclines, "")), os.ModePerm)
		if err != nil {
			log.Fatal("Error writing the file", err)
		}
		funclines = make([]string, 0)
	}

	//TODO: do the main class def file
	classdeflines := ""
	classdeflines = classdeflines + "classdef " + className + "\r"
	classdeflines = classdeflines + `properties
        Version = '1'
    end
    methods (Static)        
`
	classdeflines = classdeflines + strings.Join(funcdefs, "\r")
	classdeflines = classdeflines + `

    end
end
`
	err = ioutil.WriteFile(className+".m", []byte(classdeflines), os.ModePerm)
	if err != nil {
		log.Fatal("Error writing the Classdefninition file", err)
	}
}
