function [ac, aceach] = autocorrelate_bins(Trials, latency, binw)


% [ACF ACFeach] = autocorrelate(Trials, latency, binw)
%
% autocorrelate calculates the autocorreltion function for
% Spike trains in a set of Trials. Spikes are counted from start
%
% Returns ACF: Sum of ACs across trials.
%         ACFeach array with one ACF for each trial
%
% Arguments:  Trials Trial array
%             Latency count spikes shifted by this (in 0.1 ms
%                units)
%             binw bin-width

duration = min([Trials.End] - [Trials.Start]);
duration = floor(duration/binw)
latency = floor(latency/binw)

k = 1;
for trial = 1:length(Trials)
  atimes = floor(Trials(trial).Spikes ./ binw);
  btimes = atimes(find(atimes > latency & atimes < duration + latency));

  len = length(atimes);
  for j = 1:duration/2
    aceach(j,k) = length(intersect(btimes, mod(btimes+j, duration+latency)+latency));
  end
  k = k+1;
end

ac = sum(aceach,2) ./ (k-1);