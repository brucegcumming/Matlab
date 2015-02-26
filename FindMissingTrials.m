function [badtrial,badid] = FindMissingTrials(Expt, t, varargin)
%[badtrial,badid] = FindMissingTrials(Expt, t, varargin)
% find trials in Expt for which no elements in t fall in the trial

badtrial = [];
badid = [];

trange = Expt.Trials(end).Start(1) - Expt.Trials(1).Start(1);
if diff(minmax(t(t>0))) < trange/1000 %wrong units
    t = t .* 10000;
end

for j = 1:length(Expt.Trials)
    good(j) = sum(t> Expt.Trials(j).Start(1) & t < Expt.Trials(j).End(end));
    if good(j) == 0
        badtrial(end+1) = j;
        badid(end+1) = Expt.Trials(j).id;
    end
end
    