function LoadAllProbes(a,b, varargin)
DATA = GetDataFromFig(a);
expid = DATA.currentexpt(1);
j = 1; 
while j <= length(varargin)
j = j+1;
end
times(1) = DATA.Expts{expid}.Trials(1).Start(1)-10000;
times(2) = DATA.Expts{expid}.Trials(end).End(end) + 10000;
DATA.allexp = expid;
Spikes = DATA.AllData.Spikes;
tic;
DATA = cmb.GetAllProbeFig(DATA);
if ~isfield(DATA,'ptsize')
DATA.ptsize = 1;
end
for j = 1:length(DATA.probelist)
if DATA.state.online == 0
DATA.AllData.Spikes = cmb.GetProbeFiles(DATA, DATA.probelist(j),DATA.subprobe,'trange',times/10000);
else
filename = DATA.Expts{DATA.currentexpt(1)}.Header.loadname;
if DATA.probelist(j) > 16 || DATA.probesource(j) == 1
filename = strrep(filename,'/Expt','A/Expt');
end
[DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, filename , DATA.probevars{j},[DATA.probe DATA.subprobe]);
end
if isempty(DATA.AllData.Spikes) | isempty(DATA.Expts{expid}.gui.spks)
DATA.clustervals(j).x = [];
DATA.clustervals(j).y = [];
DATA.clustervals(j).spkrange = [0 0 ];
else
DATA = SetExptSpikes(DATA, DATA.currentexpt(1), 0);
ispk = DATA.Expts{expid}.gui.spks;
DATA.clustervals(j).x = DATA.Spikes.cx(ispk);
DATA.clustervals(j).y = DATA.Spikes.cy(ispk);
DATA.clustervals(j).spkrange = [ispk(1) ispk(end)];
subplot(6,4,j);
if DATA.densityplot
cmb.PlotXYDensity(DATA.Spikes.cx(ispk),DATA.Spikes.cy(ispk));slave

else
plot(DATA.Spikes.cx(ispk),DATA.Spikes.cy(ispk),'.','markersize',DATA.ptsize);
end
title(sprintf('%d: %d',j,length(ispk)));
end
end
subplot(6,4,2);
text(1,1.5,sprintf('Expt %d',DATA.currentexpt(1)),'units','norm');
toc
DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
set(DATA.toplevel,'UserData',DATA);

