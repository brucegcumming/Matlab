function GaussRatios(varargin)

nsd = [0:0.02:0.2];
nrep = 1000;
nsmp = 10;
j =1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nsd',3)
        j = j+1;
        nsd = varargin{j};
    elseif strncmpi(varargin{j},'nsmp',4)
        j = j+1;
        nsmp = varargin{j};
    elseif strncmpi(varargin{j},'nrep',4)
        j = j+1;
        nrep = varargin{j};
    end
    j = j+1;
end

simtype = 2;
if simtype == 1 %%resmpling simulation
nres = 100;
for j = length(nsd):-1:1
    for n = nrep:-1:1
        gnoise = randn(nsmp,1) + nsd(j);
        noise = Bresample(gnoise, nres);
        p(n)= sum(sum(noise,2) > 0)./nres;

    end
    means(j) = mean(p);
    stds(j) = std(p);
end
    binom = sqrt(means .* (1-means)./3);
subplot(2,1,1);
    plot(binom,stds,'ro');
refline(1);
subplot(2,1,2);
errorbar(nsd, means, stds);
else
    for j = length(nsd):-1:1;
        gnoise = randn(nsmp,nrep) + nsd(j);
        bnoise = randn(nsmp,nrep) + nsd(j);
        ps = mean(gnoise > 0);
        p(j) = mean(ps);
        sds(j) = std(ps);
        rev = mean(sign(gnoise) ~= sign(bnoise));
        prev(j) = mean(rev);
        cu(j) = prctile(rev,97.5);
        u = p2m(mean(ps));
        pu(j) = m2p(u+2);
        pl(j) = m2p(u-2);
        pcu(j) = prctile(ps,97.5);
        pcl(j) = prctile(ps,2.5);
        cl(j) = prctile(rev,2.5);
    end
    subplot(2,1,1);
    plot(p,sds,'ro-');
    subplot(2,1,2);
    hold off;
    if 0
    plot(p,prev);
    hold on;
    plot(p,cl,':');
    plot(p,cu,':');
    hold on;
    else
        plot(nsd,p);
        hold on;
        plot(nsd,pcu,':');
        plot(nsd,pu,'r:');
        plot(nsd,pcl,':');
        plot(nsd,pl,'r:');
    end
end

function u = p2m(p)

u = erfinv(2*p -1);

function p = m2p(u)
p = (1+erf(u))/2;
