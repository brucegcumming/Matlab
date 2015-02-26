function poisson(rates, meansd, vararagin)
%poisson(mean, measd, vararagin) Simulate rate varying poisson
%
%

nloops = 100;
ns = 20;
if length(meansd) ==1
    meansd = ones(size(rates)) * meansd;
end

for k = 1:length(rates)
    r = rates(k);
    for j = 1:nloops
        rs = r + randn(1,ns) .* meansd(k);
        rs(rs<0) = 0;
        resp(j) = sum(poissrnd(rs));
    end
    means(k) = mean(resp);
    vars(k) = var(resp)
end

plot(means, vars);

