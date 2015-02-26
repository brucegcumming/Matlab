function SetClusters(a,b,tag)
DATA = get(findobj('Tag',tag),'UserData');
DATA.spikelist = cmb.WhichClusters(DATA.toplevel);
if DATA.state.online == 2
DATA.Expts = cmb.CountTxtSpikes(DATA.Expts,DATA.probe,DATA.spikelist);
end
if DATA.state.autoreplotgraph
id = union(DATA.currentexpt(1),DATA.exabsid);
for j = id
DATA = cmb.CountSpikes(DATA,j);
end
cmb.combine('setexpplot',DATA);
end
set(DATA.toplevel,'UserData',DATA);

