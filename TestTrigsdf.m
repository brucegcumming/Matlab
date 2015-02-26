function TestTrigSdf(Expt)

nloops = 10;
ts = now;
for n = 1:nloops
for j = 1:length(Expt.Trials)
    Expt.Trials(j).Trigger = n+[100:200:20000];
    Expt.Trials(j).Trigger = 0;
end
[sdf,nspk] = trigsdfc(Expt.Trials, 100, [0:10:20000]);
end
fprintf('took %2f\n',mytoc(ts));