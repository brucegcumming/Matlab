
sd = 0.2;
noise = 5;
x = -1:0.2:1;
resp = 10 + 50 .* exp(-(x.^2)./ (2.* sd.^2));
%plot(x,resp);

shift = 0;
scale = 2;
nreps = 5;
data = repmat(resp,nreps,1);
data = reshape(data,nreps * length(x),1);
adata = data + randn(size(data)) .*noise;
bdata = (data + randn(size(data)) .* noise + shift) .* scale;
alldata = [adata bdata];
smp = reshape(alldata,nreps,length(x),2);
means = mean(smp);
stds = std(smp);
hold off;
errorbar(x,means(:,:,1),stds(:,:,1),'r');
hold on;
errorbar(x,means(:,:,2),stds(:,:,2),'b');
pvals = anova2(alldata,nreps)
