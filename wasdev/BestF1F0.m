function [ratio, details] = BestF1F0(Expt)

tf = GetEval(Expt,'tf','mode');
period = 10000/tf;

tres = PlotExpt(Expt,'noplot');
[maxrate, bestid] = max(tres.means(:));
trials = tres.ids{bestid};
durs = [Expt.Trials.End] - [Expt.Trials.Start];
dur = prctile(durs,50);
details.dur = dur;
endt = dur+500;

ncycles = floor((dur - 500)/period);
if ncycles ==0 & dur > period
    ncycyles = 1;
elseif ncycles == 0
    ratio = NaN;
    details.period = period;
    return;
end
start = endt - period * ncycles;
for j = 1:length(Expt.Trials)
    Expt.Trials(j).Trigger = 0;
end

[spks, nt, psth] = SpikeCycle(Expt.Trials(trials), start, endt, period);
f1 = famp(1:length(psth),psth,1/length(psth));
f0 = mean(psth);
ratio = f1/f0;
details.stimx = tres.x(bestid);
details.stimy = tres.y(bestid);
details.ncycles = ncycles;
details.t = 1:length(psth);
details.psth = psth;


