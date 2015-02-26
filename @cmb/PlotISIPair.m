function PlotISIPair(DATA, pair)

GetFigure(DATA.tag.spikev);
hold off;
spk = find(DATA.AllData.Spikes.times >= floor(DATA.isis(pair)));
spk = spk(1);
tsmp = [1:size(DATA.AllData.Spikes.values,2)] * DATA.AllData.Spikes.interval * 10000;
plot(tsmp,DATA.AllData.Spikes.values(spk-1,:));
[a,ta] = max(DATA.AllData.Spikes.values(spk-1,:));
[b,tb] = max(DATA.AllData.Spikes.values(spk,:));
hold on;
isi = diff(DATA.AllData.Spikes.times([spk-1 spk]));
times = isi+tsmp+20000*DATA.AllData.Spikes.interval;
isi = isi + (tb-ta) * DATA.AllData.Spikes.interval * 10000;
plot(times,DATA.AllData.Spikes.values(spk,:));
title(sprintf('Time %.0f ISI %.1f',DATA.AllData.Spikes.times(spk-1),isi));

