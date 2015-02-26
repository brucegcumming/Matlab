function LoadAllProbeParams(a,b, varargin)
DATA = GetDataFromFig(a);
expid = DATA.currentexpt(1);
j = 1; 
while j <= length(varargin)
j = j+1;
end
DATA.allexp = expid;
Spikes = DATA.AllData.Spikes;
tic;
DATA = cmb.GetAllProbeFig(DATA);
if ~isfield(DATA,'ptsize')
DATA.ptsize = 1;
end
for j = 1:length(DATA.probelist)
set(DATA.toplevel,'name',sprintf('Loading %d',DATA.probelist(j)));
drawnow;
if DATA.state.online == 0
DATA.AllData.Spikes = cmb.GetProbeFiles(DATA, DATA.probelist(j),DATA.subprobe);
else
filename = ['C:' DATA.Expts{DATA.currentexpt(1)}.Header.Name];
if DATA.probelist(j) > 16
filename = strrep(filename,'/Expt','A/Expt');
end
[DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, filename , DATA.probevars{j},[DATA.probe DATA.subprobe]);
end
nspk = length(DATA.AllData.Spikes.times);
DATA = CalcClusterVars(DATA,  1:nspk,'force');
DATA.AllClusters(j).cx = DATA.Spikes.cx(1:nspk);
DATA.AllClusters(j).cy = DATA.Spikes.cy(1:nspk);
DATA.AllClusters(j).times = DATA.AllData.Spikes.times;
DATA.AllClusters(j).codes = DATA.AllData.Spikes.codes(:,2);
end
DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
set(DATA.toplevel,'UserData',DATA);
cmb.NotBusy(DATA);

