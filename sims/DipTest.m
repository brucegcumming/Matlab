function DipTest(dx, n)
%Diptest, Check did statistic with Gaussaions,
nboot = 1000;
nreps = 10;
for k = 1:nreps;
for j = 1:length(dx)
    
xa = randn(n,1);
xb = randn(n,1)+dx(j);
[dip, pval(j,k), xlow,xup]=HartigansDipSignifTest([xa xb],nboot);
end
end
pval = mean(pval,2);
if length(dx) > 1
    plot(dx,pval);
else
end

fprintf('p%.5f ',pval);
