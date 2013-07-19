function DirRatios(varargin)

nsd = [0.5:0.5:4];
nrep = 1000;
nsmp = 10;
DI = 1;
j =1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nsd',4)
        j = j+1;
        nsd = varargin{j};
    elseif strncmpi(varargin{j},'nsmp',4)
        j = j+1;
        nsmp = varargin{j};
    end
    j = j+1;
end

ndir = 8;
dirs = 0:2*pi/ndir:2*pi;
dirs = dirs(2:end);
dsig = abs(cos(dirs)) + DI .* cos(dirs);
sins = repmat(sin(dirs),nsmp,1);
coss = repmat(cos(dirs),nsmp,1);
dsig = repmat(dsig,nsmp,1);
for j = length(nsd):-1:1
    allr = [];
    for n = nrep:-1:1
        noise = randn(nsmp,ndir);
        resp = noise .* nsd(j) + dsig;
        r = sum(sins .* resp + i * coss.*resp,2);
        p(n) = sum(angle(r) < 0)./nsmp;
        allr = [r allr];
    end
    means(j) = mean(p);
    stds(j) = std(p);
end
binom = sqrt(means .* (1-means)./nsmp);
plot(binom,stds,'ro-');
refline(1);
%errorbar(nsd, means, stds);