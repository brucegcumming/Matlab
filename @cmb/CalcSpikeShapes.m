function CalcSpikeShapes(a,b, varargin)
%
% caluclate mean spike waveform for biggest spikes
% for each expt,probe
DATA = GetDataFromFig(a);
rebuild = DATA.state.forcebuild;

probelist = DATA.probelist;
TrialVar.cx = [];
outname = strrep(DATA.datafilename,'.mat','spks.mat');
j = 1;
while j < length(varargin)
if strncmpi(varargin{j},'online',5)
outname = strrep(DATA.datafilename,'.mat','spks.mat');
end
end
if exist(outname,'file') & ~rebuild
load(outname);
DATA.MeanSpike = MeanSpike;
DATA.TrialVar = TrialVar;
DATA.TemplateScores = TemplateScores;
if exist('TemplateInfo','var')
DATA.TemplateInfo = TemplateInfo;
clear TemplateInfo;
end
else
if ~isfield(DATA,'TemlateScores')
DATA.TemplateScores = [];
end
oldstate = DATA.state.autoplotnewprobe;
DATA.state.autoplotnewprobe = 1;
for p = 1:length(probelist)
if p ~= DATA.probe
DATA = cmb.SetProbe(DATA, probelist(p));
end
DATA = cmb.CalcMeanSpike(DATA,1:length(DATA.Expts));
end
DATA.state.autoplotnewprobe = oldstate;
DATA = cmb.SaveSpikeShape(DATA,outname);
end

GetFigure('SpikeShape');
subplot(2,1,1);
imagesc(std(DATA.MeanSpike.v,[],3)');
subplot(2,1,2);
id = find(sum(DATA.TrialVar.cx) > 0);
imagesc(abs(DATA.TrialVar.cx(:,id)));
set(DATA.toplevel,'UserData',DATA);
GetFigure('ClusterShape');
PlotSpikeShapes(DATA.MeanSpike,'auto');
GetFigure('TemplateScore')
cmb.PlotTemplateScores(DATA,2);


