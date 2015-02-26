function Expt = CountEvents(Expt, t, latency, varargin);
%Expt = CountEvents(Expt, t, latency) counts number of t that fall in each
%Trial of expt
dt = t - latency;

for j = 1:length(Expt.Trials)
   id = find(dt > Expt.Trials(j).Start(1) & dt < Expt.Trials(j).End(end));
   Expt.Trials(j).count = length(id);
end