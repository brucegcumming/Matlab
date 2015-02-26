function DATA = PlotSpikes(DATA, ispk)

classify = 1;
ctype = 2;
[cx, DATA] = GetSpikeVals(DATA,ispk, DATA.AllSpikes.values(ispk,:),DATA.AllSpikes.dVdt(ispk,:), DATA.plot.clusterX, classify, []);
DATA.Spikes.cx(ispk) = cx;
%      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
nc =  length(DATA.cluster)+1;
cmb.DrawSpikeWaves(DATA, ispk, nc, 2);
nc =  length(DATA.cluster)+1;
cmb.DrawSpikeWaves(DATA, ispk, nc, 2);
[cy, DATA] = GetSpikeVals(DATA,ispk, DATA.AllSpikes.values(ispk,:),DATA.AllSpikes.dVdt(ispk,:),DATA.plot.clusterY, classify, []);
DATA.Spikes.cy(ispk) = cy;
DATA.currentspike = ispk(1);
DATA = SetSpkCodes(DATA,ispk,DATA.probe,0);
nc =  length(DATA.cluster)+1;
cmb.DrawSpikeWaves(DATA, ispk, nc, 2);

set(0,'CurrentFigure',DATA.xyfig);
for j = 0:nc
sp = find(DATA.AllData.Spikes.codes(ispk, ctype) == j);
plot(cx(sp),cy(sp),...
'.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
hold on; %% need this when called from PlotOneTrial
end

