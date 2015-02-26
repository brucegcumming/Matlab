function avg = EmAvg(Expt, varargin)


j = 1;
while j < nargin;
    j = j+1;
end

buflen = 1e10;
for trial = 1:length(Expt.Trials)
    buflen = min([buflen length(Expt.Trials(trial).Eyevals.lv)]);
end

lvs = [];
rvs = [];
lhs = [];
rhs = [];
for trial = 1:length(Expt.Trials)
    lvs(trial,:) = Expt.Trials(trial).Eyevals.lv(1:buflen)';
    rvs(trial,:) = Expt.Trials(trial).Eyevals.rv(1:buflen)';
    rhs(trial,:) = Expt.Trials(trial).Eyevals.rh(1:buflen)';
    lhs(trial,:) = Expt.Trials(trial).Eyevals.lh(1:buflen)';
end

avg(1,:) = mean(lhs,1);
avg(2,:) = mean(rhs,1);
avg(3,:) = mean(lvs,1);
avg(4,:) = mean(rvs,1);

if showplot
    plot(avgs);
end
