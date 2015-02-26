function [ispk, sspk,cx, cy] = PlotTrialSyncSpikes(DATA, times, probes, colors, varargin)

SPKMIN=22;
SPKPEAK=7;

nv = 32;
dt = 2;
j=1;
ctype = 2;
step = DATA.plot.SpikeMaxV;
if length(probes) > 3
step = 3 * DATA.plot.SpikeMaxV/4;
end
if DATA.plot.SpikeVsep > 0
step = DATA.plot.SpikeVsep;
end
if DATA.plot.syncoverlay
timemode = 0;
else
timemode = 1;
end
while j <= length(varargin)
if strncmpi(varargin{j},'lineoff',6)
j = j+1;
lineoff = varargin{j};
elseif strncmpi(varargin{j},'Probe',4)
j = j+1;
probe = varargin{j};
elseif strncmpi(varargin{j},'Trial',4)
j = j+1;
Trial = varargin{j};
end
j = j+1;
end

if isfield(DATA,'sids') && isfield(DATA,'AllSpikes')
ai = find(DATA.AllSpikes{probes(1)}.times(DATA.sids{1}) > times(1) ...
& DATA.AllSpikes{probes(1)}.times(DATA.sids{1}) < times(2));
for j = 1:length(probes)
sids{j} = DATA.sids{j}(ai);
ispks{j} = find(DATA.AllSpikes{probes(j)}.times > times(1) ...
& DATA.AllSpikes{probes(j)}.times < times(2));
itimes{j} = DATA.AllSpikes{probes(j)}.times(sids{j});
end
ispk = ispks{1};
else
for j = 1:length(probes)
if DATA.spikelist == -1
spikelist = [0 1 2 3 4];
else
spikelist = DATA.spikelist;
end
ispks{j} = find(DATA.AllSpikes{probes(j)}.times > times(1) &...
DATA.AllSpikes{probes(j)}.times < times(2) & ...
ismember(DATA.AllSpikes{probes(j)}.codes(:,2),spikelist));
itimes{j} = DATA.AllSpikes{probes(j)}.times(ispks{j});
end
iids = [];
dts{1} = zeros(size(itimes{1}));
for p = 2:length(ispks)
ids = [];
for j = 1:length(itimes{1})
id = find(abs(itimes{p}-itimes{1}(j)) < dt);
if length(id)
ids = [ids id(1)];
dts{p}(length(ids)) = (itimes{p}(id(1))-itimes{1}(j)) * nv/10;
iids = [iids j];
end
end
sids{p} = ispks{p}(ids);
end
sids{1} = ispks{1}(iids);
end
if ~isempty(sids{1})
if timemode
set(0,'CurrentFigure',DATA.timefig);
vh = DATA.tvh;
else
set(0,'CurrentFigure',DATA.svfig);
vh = DATA.svh;
end
nc = (length(probes)-1)/2;
n = DATA.svhn;
for j = 1:length(probes)
for k = 1:n
vs{(j-1)*n+k} = [];
xs{(j-1)*n+k} = [];
end
dts{j} = itimes{j} - itimes{1};
for k = 1:length(sids{j});
spk = sids{j}(k);
nl = (j-1)*n + DATA.AllSpikes{probes(j)}.codes(spk,2)+1;
if DATA.plot.syncoverlay
vs{nl} = [vs{nl} DATA.AllSpikes{probes(j)}.values(spk,:)+(j-nc-1)*step NaN];
xs{nl} = [xs{nl} [1:nv]-dts{j}(k) NaN];
else
vs{nl} = [vs{nl} DATA.AllSpikes{probes(j)}.values(spk,:) NaN];
xs{nl} = [xs{nl} (k-1)*nv + [1:nv]-dts{j}(k) NaN];
end
end


for k = 1:n
nl = (j-1)*n + k;
if nl > length(DATA.svh)
vh(nl) = line('Xdata' , xs{nl}, 'Ydata', vs{nl});
else
if ishandle(DATA.svh(j))
set(vh(nl),'Xdata' , xs{nl}, 'Ydata', vs{nl});
end
end
end
end
if length(xs{1})
set(gca,'xlim',[1 max(xs{1})]);
end
end
sspk = cat(2,sids{:});
nspk = length(sids{1});

title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f %d/%d spks',abs(Trial.Trial),...
Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,nspk,length(ispk)));

cx = [];
cy = [];
if DATA.plot.synccluster == 1 %min/min
[cx, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),SPKMIN, 0,[]);
[cy, DATA] = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),SPKMIN, 0,[]);
elseif DATA.plot.synccluster == 4 %PCA1 spk1 vs 2
cx = DATA.pca(3,ai);
cy = DATA.pca(4,ai);
elseif DATA.plot.synccluster == 5 %PCA2 1 vs 2
cx = DATA.pca(5,ai);
cy = DATA.pca(6,ai);
elseif DATA.plot.synccluster == 6 %PCA 1 vs PCA2
cx = DATA.pca(1,ai);
cy = DATA.pca(2,ai);
elseif DATA.plot.synccluster == 7 %PCA vs xcorr
cx = DATA.pca(1,ai);
cy = DATA.AllSpikes{probes(1)}.values(sids{1},:) .* DATA.AllSpikes{probes(2)}.values(sids{2},:);
cy = sum(cy');
elseif DATA.plot.synccluster == 8 %Current Cluster X, for each probe
cx  = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterX, 1,[]);
cy  = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),DATA.plot.clusterX, 1,[]);
elseif DATA.plot.synccluster == 9 %Current Cluster Y, for each probe
cx  = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterY, 1,[]);
cy  = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),DATA.plot.clusterY, 1,[]);
elseif DATA.plot.synccluster == 10 %Current Cluster Plot for sum
V = (DATA.AllSpikes{probes(1)}.values(sids{1},:)+DATA.AllSpikes{probes(2)}.values(sids{2},:))/2;
dV= ( DATA.AllSpikes{probes(1)}.dVdt(sids{1},:)+ DATA.AllSpikes{probes(2)}.dVdt(sids{2},:))/2;
cx  = GetSpikeVals(DATA,sids{1}, V, dV, DATA.plot.clusterX, 1);
cy  = GetSpikeVals(DATA,sids{1}, V, dV, DATA.plot.clusterY, 1);
else
[cx, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterX, 0,[]);
[cy, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterY, 0,[]);
end
if length(cx) & length(cy)
set(0,'CurrentFigure',DATA.xyfig);
for j = 0:length(DATA.cluster{probes(1)})+1
sp = find(DATA.AllSpikes{probes(1)}.codes(sids{1}, ctype) == j);
plot(cx(sp),cy(sp),...
'.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
hold on; %% need this when called from PlotOneTrial
end
end
drawnow;

