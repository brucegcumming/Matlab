function outname = Expt2FileName(DATA,Expt, cluster)

[expname, extype, suff] = Expt2Name(Expt);
it = strmatch(extype,DATA.expstrs,'exact');
if length(it)
expname = strrep(expname,extype,DATA.expnames{it(1)});
end
if DATA.state.includeprobename
cs = ['.p' num2str(DATA.probe) 'c'];
else
cs = '.c';
end
outname = [strrep(DATA.datafilename,'.mat',cs) num2str(cluster) '.' expname suff '.mat'];


