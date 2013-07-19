function [isis, trials, spkids] = CalcISI(Trials, varargin)
%
%[isis, trials, spkids] = CalcISI(Trials, varargin)
%spkids are actually times
latency = 500;

isis = [];
trials = [];
spkids = [];
for j = 1:length(Trials)
    duration = Trials(j).End(end) - Trials(j).Start(1);
 spks = find(Trials(j).Spikes > latency & ...
     Trials(j).Spikes < duration+latency);
 isis = [isis diff(Trials(j).Spikes(spks)')];
 trials = [trials ones(size(spks(2:end)')) * j];
 spkids = [spkids Trials(j).Spikes(spks(2:end))'+Trials(j).Start(1)];
end