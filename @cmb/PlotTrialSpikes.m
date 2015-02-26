function DATA = PlotTrialSpikes(DATA, itrial, colors, clusters)
spka = DATA.AllData.SpikeIdx(itrial,1);
spkb = DATA.AllData.SpikeIdx(itrial,2);
Spks = DATA.AllData.Spikes;
if spka
for spk = spka:spkb;
adc = Spks.values(spk,:);
plot(adc,'color',colors{Spks.codes(spk,1)+1});
hold on;
energy  = sum(diff(adc).^2);
DATA.Spikes.energy(spk)= energy;
svar(spk) = var(adc);
DATA.Spikes.vw(spk) = svar(spk)/energy;
end
title(sprintf('Trial %d (id%d)',DATA.AllData.Trialids(itrial),...
DATA.AllData.Trialids(itrial)));

drawnow;


set(0,'CurrentFigure',DATA.xyfig);
for j = 1:length(clusters)
sp = intersect(clusters{j},[spka:spkb]);
plot(DATA.Spikes.energy(sp),DATA.Spikes.vw(sp),'.','color',colors{j});
end
end


