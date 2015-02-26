function LoadSpikeTimes(a,b, varargin)

DATA = GetDataFromFig(a);
eid = DATA.currentexpt(1);
ns = ones(size(DATA.probelist));
for j = 1:length(DATA.probelist)
DATA.AllSpikes{j}.times = [];
end
for e = 1:length(DATA.Expts)
DATA.currentexpt = e;
sfile = cmb.ClusterFile(DATA, e, 'codes');
if exist(sfile,'file')
load(sfile);
for j = 1:length(DATA.probelist)
n = length(codes{j});
if n
DATA.AllSpikes{j}.codes(ns(j):ns(j)+n-1,2) = codes{j};
DATA.AllSpikes{j}.times(ns(j):ns(j)+n-1,1) = times{j};
ns(j) = ns(j)+n;
end
end
end
DATA.Expts{e}.gui.classified = 1;
end
for j = 1:length(DATA.probelist)
if isfield(DATA.AllSpikes{j},'codes')
ns(j) = sum(DATA.AllSpikes{j}.codes(:,2) > 0);
else
ns(j) = NaN;
end
end
DATA.currentexpt = eid;
DATA.state.nospikes = 1;
DATA.state.recount = 1;
DATA.state.showspikes = 0;
DATA.state.showspkxy = 0;
cmb.SetGui(DATA);
set(DATA.toplevel,'UserData',DATA);

