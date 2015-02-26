function [allisi, hist] = SpikeISI(Expt)

allisi = [];
for j = 1:length(Expt.Trials);
    isis = diff(Expt.Trials(j).Spikes)';
    allisi = [allisi isis];
end

[hist.y,hist.x] = smhist(allisi);