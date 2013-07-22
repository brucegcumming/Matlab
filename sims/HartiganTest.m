function HartiganTest(varargin)


distances = [0:1:6];
npts = [100 1000 10000 100000; 100 1000 10000 100000];
sdratio = 1;
nreps = 2;
plottype = 'dips';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sdratio',3)
        j = j+1;
        sdratio = varargin{j};
    elseif strncmpi(varargin{j},'npts',4)
        j = j+1;
        npts = varargin{j};
    elseif strmatch(varargin{j},{'bii','dips','hist'})
        plottype = varargin{j};
    end
    j = j+1;
end

for ni = 1:nreps
for k = 1:size(npts,2)
n = npts(1,k);
if size(npts,1) > 1
m = npts(2,k);
else
    m = n;
end

    for j  = 1:length(distances)
    
a = randn(1,n).*sdratio;
b = randn(1,m)+distances(j);
cs = sort(cat(2,a,b));

[dips(ni,j,k),a, b, c, d, e,f,g ] = HartigansDipTest(cs);
bii(ni,j,k) = (1 + skewness(cs).^2)./(kurtosis(cs)+3);
kurt(ni,j,k) = kurtosis(cs);
skew(ni,j,k) = skewness(cs);
nc(ni,j,k) = length(cs);
end
end
end
scaledip = 0;
if scaledip
dips = squeeze(mean(dips.*sqrt(nc)));
else
dips = squeeze(mean(dips));
end
bii = squeeze(mean(bii));
if strcmp(plottype,'dips')
plot(distances,dips);
elseif strcmp(plottype,'bii')
plot(distances,bii);
elseif strcmp(plottype,'hist')
hist(cs);
end
legend(num2str(npts(1,:)'));