function ReloadClusters(a,b, varargin)
DATA = GetDataFromFig(a);
probes = DATA.probe;
probe = DATA.probe;
reclassify = 0;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'all',3)
probes = DATA.probelist;
elseif strncmpi(varargin{j},'relcassify',5)
reclassify = 1;
end
j = j+1;
end
for j = 1:length(probes)
DATA.probe = probes(j);
DATA = cmb.LoadClusters(DATA,cmb.ClusterFile(DATA));
if reclassify
DATA = cmb.ReClassifyAll(DATA,'probes',probes(j));
cmb.SetExptClusters(DATA);
end
end
DATA.probe = probe;
set(DATA.toplevel,'UserData',DATA);


