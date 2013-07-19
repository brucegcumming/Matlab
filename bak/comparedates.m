function list = comparedates(dir1,dir2,whichdir)

%  compares files in two directories.  For files with the same name, it
%  reports if the file in the specified directory 'whichdir' (either 1 or
%  2) is the newest.

A=dir(dir1);
B=dir(dir2);

count=0;
list={};
for i=1:length(A);
    for j=1:length(B)
        if strcmp(A(i).name,B(j).name)
            if whichdir==1;
                if A(i).datenum > B(j).datenum
                    count=count+1;
                    list{count}= B(j).name;
                end
            elseif whichdir==2;
                if A(i).datenum < B(j).datenum
                    count=count+1;
                    list{count}= B(j).name;
                end
            end
        end
    end
end
