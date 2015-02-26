function t = GetSpkTimes(E, varargin)
%extract spike times from and expt struct.

T = E.Trials;
t = [];
for j = 1:length(T)
    t = [t T(j).Spikes'+T(j).Start(1)];
end