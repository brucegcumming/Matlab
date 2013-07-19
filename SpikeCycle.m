function [spikes, nt, psth]  = SpikeCycle(Trials, tmin, tmax, period)
%[spikes, nt]  = SpikeCycle(Trials, tmin, tmax, period)
%Builds a list of spike times foldeded periodically into 'period'
%Only spikes after time tmin, and before tmax (relative to Trial.Trigger) are
%allowed in.

spikes = [];
nt = 0;
for j = 1:length([Trials.Trial])
  for t = [Trials(j).Trigger];
    spk = [Trials(j).Spikes'-t];
    spks = spk(find(spk > tmin & spk < tmax));
    spikes = [spikes mod(spks, period)];
    nt = nt + 1;
  end
end
spikes = round(spikes);
for j = 1:period
    psth(j) = sum(spikes == j);
end
