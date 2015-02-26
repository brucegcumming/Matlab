function [Spikes, Expt] = GetExptSpikeTimes(path)

d = mydir([path '*.mat']);

for j = 1:length(d)
    load(d(j).name);
    spiket = [];
    for t = 1:length(Expt.Trials)
        T = Expt.Trials(t);
        spiket = cat(1,spiket, T.Spikes+T.Start(1));
    end
    p = GetProbeFromName(d(j).name);
    Spikes{p} = spiket;
end


