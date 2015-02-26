function Expt = ReCountExptSpikes(Expt, latency, varargin)





for j = 1:length(Expt.Trials)
    durs(j) = Expt.Trials(j).End(end) - Expt.Trials(j).Start(1);
end

duration = prctile(durs,95);

for j = 1:length(Expt.Trials)
    Expt.Trials(j).count = sum(Expt.Trials(j).Spikes > latency & Expt.Trials(j).Spikes < latency+duration);
end
