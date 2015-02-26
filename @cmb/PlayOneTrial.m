function PlayOneTrial(DATA, a, b, varargin)
setgui = 1;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'nosetgui',6)
setgui = 0;
end
j =j+1;
end

%DATA = cmb.combine('getstate');
expid = DATA.currentexpt(1);

if b < 0  & DATA.currenttrial > 1 %step back
DATA.currenttrial = DATA.currenttrial-1;
while DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).Trial < 0 & DATA.currenttrial > 1 
DATA.currenttrial = DATA.currenttrial-1;
end
elseif b > 0 & DATA.currenttrial < length(DATA.Expts{DATA.currentexpt(1)}.Trials) %step back
DATA.currenttrial = DATA.currenttrial+1;
while DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).Trial < 0 &  ...
DATA.currenttrial < length(DATA.Expts{DATA.currentexpt(1)}.Trials) 
DATA.currenttrial = DATA.currenttrial+1;
end
elseif b == 0 % step to a trial on list
DATA.currenttrial = find([DATA.Expts{DATA.currentexpt(1)}.Trials.Trial] == a);    
DATA.currenttrial = a;    
end
if DATA.currenttrial <= 0
DATA.currenttrial = 1;
end
Trial = DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial);
if Trial.Trial < 0 && b == 0 && setgui %%manually select -Trial = include again
DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).Trial = abs(Trial.Trial);
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
end
itrial = find(DATA.AllData.Trialids == abs(Trial.Trial));
set(0,'CurrentFigure',DATA.svfig);
hold off;
%DATA = cmb.PlotTrialSpikes(DATA, itrial, mycolors, DATA.clusters);
if DATA.state.recut
nc = size(DATA.cluster,1)+1;
else
nc = DATA.s2clusters;
end
Trial.ed = GetEval(DATA.Expts{DATA.currentexpt(1)},'ed',DATA.currenttrial);
probes =  [DATA.probe DATA.xprobes];
if (length(probes) > 1 & DATA.plot.syncoverlay == 0) || DATA.plot.showwave
Aargs = {'timemode'};
if  ~isfield(DATA,'timefig')
DATA.timefig = GetFigure('SpikeTime');
end
set(0,'CurrentFigure',DATA.timefig);
else
Aargs = {};
end


if DATA.plot.voltxy
set(gca,'xlim',[-5 5],'ylim',[-5 5]);
end
if DATA.syncsign < 2
cmb.PlotTrialSyncSpikes(DATA, [Trial.Start(1) Trial.End(end)], [DATA.probe DATA.xprobes], mycolors,'Trial',Trial);
else
for j = 2:length(probes)
DATA = cmb.APlotTrialSpikes(DATA, [Trial.Start(1)-DATA.state.preperiod Trial.End(end)+DATA.state.postperiod], mycolors, nc, 1,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(DATA.svhn), Aargs{:});
end
DATA = cmb.APlotTrialSpikes(DATA, [Trial.Start(1)-DATA.state.preperiod Trial.End(end)+DATA.state.postperiod], mycolors, nc, 1,'Trial',Trial,'probe',[probes(1) 1 length(probes)],Aargs{:});
set(0,'CurrentFigure',DATA.svfig);
if isfield(DATA,'CellList') && size(DATA.CellList,2) >= Trial.Trial
for j = 1:length(probes)
id = find(DATA.CellList(:,Trial.Trial) == probes(j));
if id
text(10,DATA.plot.currentoffset(j),['cell' num2str(id(1))]);
end
end
end
end
if length(probes) > 1 & isfield(DATA,'timefig') %also show spikes over time
set(0,'CurrentFigure',DATA.timefig);
vh = DATA.tvh;
if DATA.plot.timebyspikeprobe > 0
p = DATA.probe;
tw = 100;
set(gca,'xlim',[0 tw*2]);
[sid, stimes, codes] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], p,[]);
sid = sid(codes ==1);
stimes = stimes(codes ==1);
dt = [-110:10:110];
dts = [];
for j = 1:length(probes)
[sids{j}, alltimes{j} allcodes{j}] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], probes(j),[]);
voff(j) = DATA.plot.SpikeVsep*(j-1);
diffs{j} = [];
for k = 1:length(stimes)
dts(k,:) = hist(stimes(k) - alltimes{j},dt);
end
corrs(j,:) = sum(dts);
end
nsmp = size(DATA.AllSpikes{probes(1)}.values,2);
hscale = DATA.AllSpikes{probes(1)}.interval * 10000;
hold off;
for k = 1:length(sid)
ts = stimes(k);
for j = 1:length(probes)
x = [];
y = [];
id = find(alltimes{j} > ts-tw & alltimes{j} < ts+tw);
for t = 1:length(id)
x = (alltimes{j}(id(t))-ts+tw) +[1:nsmp]*hscale;
y = DATA.AllSpikes{probes(j)}.values(sids{j}(id(t)),:)+voff(j);
plot(x,y,'color',DATA.spkcolor{allcodes{j}(id(t))+1});
diffs{j} = [diffs{j} alltimes{j}(id(t))-ts];
hold on;
end
end
%        hold off;
set(gca,'xlim',[0 tw*2],'ylim',[-DATA.plot.SpikeMaxV max(voff)+DATA.plot.SpikeMaxV]);
drawnow;
pause(0.05);
end
DATA.currentspike = sid(1);
else
for j = 1:length(probes)
DATA = cmb.APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 0,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(1+DATA.svhn),'timemode');
end
end
end
if DATA.state.uselfp
GetFigure('TrialLFP');
cmb.PlotLFPRaw(DATA.state,Trial,DATA.Expts{expid}.Header.LFPsamplerate);
end
cmb.CheckSpoolButton(DATA);
set(DATA.toplevel,'UserData',DATA);
set(0,'CurrentFigure',DATA.svfig);
tid = DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).id;
if setgui
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'value',DATA.currenttrial);
end

