function DATA = PlotAllProbeXY(a,b)

if isstruct(a)
DATA = a;
else
DATA = GetDataFromFig(a);
end
cx = [];
cy = [];
DATA = cmb.GetAllProbeFig(DATA);
set(gcf, 'KeyPressFcn',@cmb.KeyPressed);
set(gcf, 'WindowButtonDownFcn',@cmb.ButtonPressed);
set(gcf, 'WindowButtonMotionFcn',@cmb.ButtonDragged);
set(gcf, 'WindowButtonUpFcn',@cmb.ButtonReleased);
densityplot = cmb.GetCheck('AllDensity',DATA.figs.allprobes);
DATA.alldensityplot = densityplot;
if isfield(DATA.plot,'useprobe') & sum(DATA.plot.useprobe) > 1
probelist = find(DATA.plot.useprobe);
else
probelist = DATA.probe;
end
[nr,nc] = Nsubplots(length(probelist));
if length(probelist) == 24
nr = 6;
nc = 4;
end
times = [DATA.Expts{DATA.allexp}.Header.trange];
DATA.ptsize = 1;
oldprobe = DATA.probe;
if isfield(DATA.plot,'useprobe')
probes = find(DATA.plot.useprobe);
else
probes = DATA.probe;
end
for j = 1:length(probes)
DATA.probe = probes(j);
p = DATA.probe;
if isfield(DATA,'AllSpikes') 
ispk = FindSpikes(DATA, times, p,[]);
if length(ispk)
DATA = CalcClusterVars(DATA,  ispk,'probe',p);
%            DATA.AllSpikes{p}.cx = DATA.Spikes.cx;
%            DATA.AllSpikes{p}.cy = DATA.Spikes.cy;
cx = DATA.Spikes.cx(ispk);
cy = DATA.Spikes.cy(ispk);
DATA.clustervals(p).spkrange = [ispk(1) ispk(end)];
DATA.AllSpikes{p}.spklist = ispk;
end
elseif isfield(DATA,'AllClusters')
if iscell(DATA.AllClusters)
e = DATA.currentexpt(1);
ispk = find(DATA.AllClusters{e}(j).times > times(1) &...
DATA.AllClusters{e}(j).times < times(2));
cx = DATA.AllClusters{e}(j).cx(ispk);
cy = DATA.AllClusters{e}(j).cy(ispk);
else
ispk = find(DATA.AllClusters(j).times > times(1) &...
DATA.AllClusters(j).times < times(2));
cx = DATA.AllClusters(j).cx(ispk);
cy = DATA.AllClusters(j).cy(ispk);
DATA.AllClusters(j).spklist = ispk;
end
else
ispk = DATA.clustervals(j).spkrange(1):DATA.clustervals(j).spkrange(2);
cx = DATA.clustervals(j).x;
cy = DATA.clustervals(j).y;
end
subplot(nr,nc,j);
hold off;
dprime = 0;
if densityplot
cmb.PlotXYDensity(cx,cy);
else
[a, dprime] = SetSpkCodes(DATA,ispk,DATA.probe,2);
%            plot(cx,cy,'.','markersize',1);
end
set(gca,'UserData',DATA.probelist(j));
hold on;
if cmb.iscluster(DATA.cluster,1,p)
DATA = cmb.DrawClusters(DATA, DATA.cluster,0);
[xr, yr] = cmb.ClusterRange(DATA.cluster,p);
else
xr = [0 0];
yr = [0 0];
end
if DATA.plot.autoscale && length(cx)
%or max cluster
DATA = cmb.SetXYRanges(DATA,cx, cy);
end
set(gca,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
if DATA.plot.prettyfigs
title(sprintf('Probe %d',p));
else
title(sprintf('%d: %d,%.2f',p,length(ispk),dprime));
end
end
DATA.probe = oldprobe;
subplot(nr,nc,1);
text(2,1.4,sprintf('Expt %d ids %d-%d, Trials %d-%d',DATA.allexp,DATA.Expts{DATA.allexp}.Trials(1).id,...
DATA.Expts{DATA.allexp}.Trials(end).id,DATA.Expts{DATA.allexp}.Trials(1).Trial,...
DATA.Expts{DATA.allexp}.Trials(end).Trial),'units','norm');

