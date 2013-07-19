function R = CheckExpFit(snr, varargin)

nreps = 100;
showplot = 1;
base = 1;
tau = 400;
x = 1:1000;

for k = 1:length(snr)
for j = nreps:-1:1
y = base + (exp(-x./tau) + randn(size(x)) .*snr(k));
fit = FitExp(x,y,'freebase');
R.taus(j,k) = fit.tau;
R.bases(j,k) = fit.params(3);
R.peaks(j,k) = fit.params(2);
if showplot >1 || showplot ==1 && j == nreps
    plot(x,y);
hold on;
z = FitExp(x,fit.params,'eval');
plot(x,z,'r');
hold off;
end
end
end
