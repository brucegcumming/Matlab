function cp = ROCsims(ntrials);
signal = 0.4455;
nreps = 1000;
for j = 1:nreps
rnd = randn(ntrials,2);
rnd(:,2) = rnd(:,2)+signal;
cps(j) = CalcCP(rnd(:,2),rnd(:,1));
end
cp = mean(cps);