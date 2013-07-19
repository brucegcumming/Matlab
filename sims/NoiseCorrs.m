function [xcs, details] = NoiseCorr(varargin)

rate = 100;
ff = 1;
ntrials = 100;
nreps = 100;
xtalk = 0.1;
mode = 1;
swapmode = 1;
corr = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'corr',4)
        j = j+1;
        corr = varargin{j};
    elseif strncmpi(varargin{j},'fano',4)
        j = j+1;
        ff = varargin{j};
    elseif strncmpi(varargin{j},'gauss',5)
        mode = 0;
    elseif strncmpi(varargin{j},'ntrials',5)
        j = j+1;
        ntrials = varargin{j};
    elseif strncmpi(varargin{j},'poiss',5)
        mode = 1;
    elseif strncmpi(varargin{j},'rate',4)
        j = j+1;
        rate = varargin{j};
    elseif strncmpi(varargin{j},'xtalk',5)
        j = j+1;
        xtalk = varargin{j};
    end
    j = j+1;
end

if mode == 1
    rate = rate/5000;
    pc = corr * rate;
allcounts = [];
allmix = [];

    
for j = 1:nreps
    if ff > 1
        psa = poissrnd(rate/2,1000,ntrials);
        psb = poissrnd(rate/2,1000,ntrials);
    else
        psa = poissrnd(rate,1000,ntrials);
        psb = poissrnd(rate,1000,ntrials);
    end
    if pc ~= 0
        psc = poissrnd(pc,1000,ntrials);
        psa(psc == 1) = 1;
        psb(psc == 1) = 1;
    end
    if ff > 1
        tb = sum(psa);
        tc = sum(psb);
        tc = tc + round(rand(size(tc))); %make sure the are not all even
        tb = tb + round(rand(size(tb))); %make sure the are not all even
        for k = 1:ntrials
            id = ceil(rand(1,tb(k)) .* 1000);
            psa(id,k) = 1;
            id = ceil(rand(1,tc(k)) .* 1000);
            psb(id,k) = 1;
        end
    end
    aid = find(psa);
    bid = find(psb);
    mixa = psa;
    mixb = psb;

    nx = round(length(bid) * xtalk);
    xid = randperm(length(bid));
    mixa(bid(xid(1:nx))) = 1;
    if swapmode == 1
        mixb(bid(xid(1:nx))) = 0;
    end

    nx = round(length(aid) * xtalk);
    xid = randperm(length(aid));
    mixb(aid(xid(1:nx))) = 1;
    if swapmode == 1
        mixa(aid(xid(1:nx))) = 0;
    end

    counts(1,:) = sum(psa);
    counts(2,:) = sum(psb);
    mix(1,:) = sum(mixa);
    mix(2,:) = sum(mixb);
    xc = corrcoef(mix');
    xcs(j) = xc(1,2);
    xc = corrcoef(counts');
    uxcs(j) = xc(1,2);
    Im = mean(mix);
    dm = abs(mix(1,:)-Im + i * (mix(2,:)-Im));
    Ic = mean(counts);
    dc = abs(counts(1,:)-Ic + i * (counts(2,:)-Ic));
    dd(j) = mean(dc)-mean(dm);
    allcounts = [allcounts counts];
    allmix = [allmix mix];
end
else
for j = 1:nreps
    counts = rate .* randn(2,ntrials) * sqrt(rate);
    mix(1,:) = counts(1,:) + xtalk .* counts(2,:);
    mix(2,:) = counts(2,:) + xtalk .* counts(1,:);
    if swapmode == 1
         mix(1,:) = mix(1,:) - xtalk .* counts(1,:)
         mix(2,:) = mix(2,:) - xtalk .* counts(2,:)
    end
    xc = corrcoef(mix');
    xcs(j) = xc(1,2);
end
end
mean(xcs);
details.uxcs = uxcs;
details.allcounts = allcounts;
details.allmix = allmix;
