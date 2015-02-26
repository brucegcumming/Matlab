function ExcludeTrials(DATA, mousept)
if DATA.state.plotseq == 2
id = find([DATA.Expt.Trials.Trial] > mousept.start(1,1) & ...
[DATA.Expt.Trials.Trial] < mousept.finish(1,1));
ex = [DATA.Expt.Trials(id).Trial];
if isfield(DATA.Expt,'ExcludeCluster')
DATA.Expt.ExcludeCluster{1} = union(DATA.Expt.ExcludeCluster{1}, ex);
else
DATA.Expt.ExcludeCluster{1} = ex;
end
end


hold off;    
DATA.Expt = cmb.PlotCombined(DATA, DATA.Expt);
set(DATA.toplevel,'UserData',DATA);

