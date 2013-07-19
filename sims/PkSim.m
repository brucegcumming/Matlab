function [pred, details] = PkSim(varargin)
%PkSim simulates makeing a Psych Kernel using a finite set of randomimages.
%similar to the ORBW task. 
%PkSim('npix',16,'ntrials',100000, 'oned','noise',[4 6 8 12 16 24 30 60 150 300],'bias',[0.5 0.55 0.6 0.65])
% shows the interaction between bias an predicability on the relationship
% between consistency and the prediction perfomrnace of the true kernel

npermute = 10;
bias = 0.5;
npix = 10;
ntrials = 50;
nseeds = 100;
nframes = 120;
noisesd = 1;
mode = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'bias',4)
        j = j+1;
        bias = varargin{j};
    elseif strncmpi(varargin{j},'npermute',4)
        j = j+1;
        npermute = varargin{j};
    elseif strncmpi(varargin{j},'npix',4)
        j = j+1;
        npix = varargin{j};
    elseif strncmpi(varargin{j},'ntrials',4)
        j = j+1;
        ntrials = varargin{j};
    elseif strncmpi(varargin{j},'oned',4)
         mode = 2;
    elseif strncmpi(varargin{j},'noise',4)
        j = j+1;
        noisesd = varargin{j};
    end
    j = j+1;
end
psum = zeros(npix,npix);
nsum = psum;

if mode == 2
    [pred, details] = RandPkSim(bias, npix, ntrials, noisesd);
    return;
end
tk = rand(npix, npix);
tk = tk - mean(tk(:));
for j = 1:nseeds
impwr{j} = rand(npix, npix)-0.5;
ims(j,:,:) = impwr{j};
end
ses = ceil(rand(nframes,ntrials) * nseeds);
for j = 1:ntrials
    fpwr{j} = squeeze(mean(ims(ses(:,j),:,:)));
    tresp(j) = sum(fpwr{j}(:) .* tk(:));
end


%make choices using a real kernel
choices = tresp > prctile(tresp,100*bias);
%choices = rand(1,ntrials) > bias;
pn = sum(choices);
nn= sum(~choices);

for j = 1:ntrials
    Expt.Trials(j).Seedseq = ses(:,j);
if choices(j)
    psum = psum + fpwr{j};
else
    nsum = nsum+fpwr{j};
end
end
pid = find(choices);
nid = find(~choices);
zid = 1:length(choices);
tid = [find(ismember(zid,pid)) find(ismember(zid,nid))];
[kp, a] = CalcPKernel(Expt, pid, nid, impwr, fpwr(tid));
a = PredictChoice(kp, fpwr, ses, (choices*2)-1, a);
kp = fftshift(kp);
kl = psum./pn-nsum./nn;


for j = 1:ntrials
    scores(j) = sum(fpwr{j}(:) .* kl(:));
    ascores(j) = sum(fpwr{j}(:) .* kp(:));
    if choices(j)
        pi = psum-fpwr{j};
        ni = nsum;
        ipn = pn - 1;
        inn = nn;
    else
        ni = nsum-fpwr{j};
        pi = psum;
        inn = nn - 1;
        ipn = pn;
    end
    ik = (pi./ipn)-(ni./inn);
    iscores(j) = sum(fpwr{j}(:) .* ik(:));
end

obias = nn/(pn+nn);
crit = prctile(scores,100 * obias);
icrit = prctile(iscores,100 * obias);
pc = scores >= crit;

pred(1) = sum(pc == choices)./(nn+pn);
pred(2) = obias * obias + (1-obias) * (1-obias);
pc = iscores >= icrit;
pred(3) = sum(pc == choices)./(nn+pn);


id = randperm(ntrials);
ichoices = choices(id);
psum = zeros(npix,npix);
nsum = psum;

for j = 1:ntrials
if ichoices(j)
    psum = psum + impwr{j};
else
    nsum = nsum+impwr{j};
end
end

ipid = find(ichoices);
inid = find(~ichoices);
tid = [find(ismember(zid,ipid)) find(ismember(zid,inid))];
[kp, b] = CalcPKernel(Expt, ipid, inid, impwr, fpwr(tid));
b = PredictChoice(kp, fpwr, ses, (ichoices*2)-1, b);
for j = 1:ntrials
    if ichoices(j)
        pi = psum-fpwr{j};
        ni = nsum;
        ipn = pn - 1;
        inn = nn;
    else
        ni = nsum-fpwr{j};
        pi = psum;
        inn = nn - 1;
        ipn = pn;
    end
    ik = pi./ipn-ni./inn;
    iscores(j) = sum(fpwr{j}(:) .* ik(:));
end
icrit = prctile(iscores,100 * obias);
pc = iscores >= icrit;
pred(4) = sum(pc == choices)./(nn+pn);


for k = 1:npermute
id = randperm(ntrials);
ichoices = choices(id);
psum = zeros(npix,npix);
nsum = psum;

for j = 1:ntrials
if ichoices(j)
    psum = psum + fpwr{j};
else
    nsum = nsum+fpwr{j};
end
end

ipid = find(ichoices);
inid = find(~ichoices);
tid = [find(ismember(zid,ipid)) find(ismember(zid,inid))];
[kp, b] = CalcPKernel(Expt, ipid, inid, impwr, fpwr(tid));
b = PredictChoice(kp, fpwr, ses, (ichoices*2)-1, b);
ik = psum./pn-nsum./nn;

for j = 1:ntrials
    if ichoices(j)
        pi = psum-fpwr{j};
        ni = nsum;
        ipn = pn - 1;
        inn = nn;
    else
        ni = nsum-fpwr{j};
        pi = psum;
        inn = nn - 1;
        ipn = pn;
    end
    ik = pi./ipn-ni./inn;
    iscores(j) = sum(fpwr{j}(:) .* ik(:));
end
icrit = prctile(iscores,100 * obias);
pc = iscores >= icrit;
pperf(k,1) = sum(pc == ichoices)./(nn+pn);
pperf(k,2) = b.predchoice(1);
end
pc = iscores >= icrit;
pred(4) = sum(pc == ichoices)./(nn+pn);
pred(4) = mean(pperf(:,1));
pred(5) = mean(pperf(:,2));

function [kernel, details] = CalcPKernel(Expt, pid, nid, pwr, tpwr)

checktrials = 1;

[a,b] = Counts([Expt.Trials(pid).Seedseq]);
[c,d] = Counts([Expt.Trials(nid).Seedseq]);

psum = zeros(size(pwr{1}));
nsum = zeros(size(pwr{1}));
for j = 1:length(a)
    se = b(j);
    psum = psum + pwr{se} .* a(j);
end
for j = 1:length(c)
    nsum = nsum + pwr{d(j)} .* c(j);
end
pim = psum./sum(a);
nim = nsum./sum(c);
details.pim = pim;
details.nim = nim;
details.psum = psum;
details.nsum = nsum;
details.nframes = [sum(a) sum(c)];

kernel = fftshift(pim-nim);

aim = (nsum+psum)./(sum(a)+sum(c));
details.allmean = aim;

if checktrials
psum = zeros(size(pwr{1}));
nsum = zeros(size(pwr{1}));
for j = 1:length(pid)
    psum = psum + tpwr{j}; %in order as zid = [pid nid]
end
for j = length(pid)+1:length(nid)+length(pid)
    nsum = nsum + tpwr{j}; %in order as zid = [pid nid]
end
tk = psum-nsum;    
end

function details = PredictChoice(kernel, tpwr, seeds, choices, details, varargin)
%
%tpwr is not the mean power spectrum for each trial. Avoids a lot of
%expensive arracy indexing
jacknife = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nojacknife',3)
        jacknife = 0;
    end
    j = j+1;
end

sk = fftshift(kernel);
nf = size(seeds,1);
for j = 1:length(tpwr)
    resp(j) = sum(tpwr{j}(:) .*sk(:));
end
aresp = resp; %copy so resp can be used for jacknife
if jacknife
%   fsums = sum(pwr(seeds,:,:)); Ntrials * 200 x 256 x 256 = too much
%   memory
    for j = 1:length(choices)
        if choices(j) == 1
        psum = details.psum - tpwr{j} .*nf;
        pn = details.nframes(1) - nf;
        nsum = details.nsum;
        nn = details.nframes(2);
        else
        nsum = details.nsum - tpwr{j}.*nf;
        nn = details.nframes(2)-nf;
        psum = details.psum;
        pn = details.nframes(1);
        end
        ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
        resp(j) = sum(tpwr{j}(:).*ski(:));
    end
else
%calculate projection of each frame onto kernel

end
% set a low criterion to get more + choices
p = sum(choices==-1)./length(choices);
crit = prctile(resp,p.*100);
details.criterion(1) = crit;
pchoice = ((resp >= crit) * 2)-1;
icorrect = choices == pchoice;
correct = (resp >= crit & choices ==1) | (resp < crit & choices == -1);
pc(1) = sum(correct)./length(choices);
pc(2) = p.*p + (1-p).*(1-p); %performance of random prediction
if jacknife
    crit = prctile(aresp,p.*100);
    details.criterion(2) = crit;
    details.pkaresp = aresp;
    correct = (aresp >= crit & choices ==1) | (aresp < crit & choices == -1);
    pc(3) = sum(correct)./length(choices);
end
details.predchoice = pc;
details.pkresp = resp;
pid = find(choices ==1);
nid = find(choices == -1);
psum = details.psum;
nsum = details.nsum;
ski = sk;
j = 1;
scores(1) = sum(ski(:).*sk(:));
while scores(j) > 0 & j < 100
    if length(pid) > length(nid)
        [a, id] = max(aresp(pid));
        psum = psum - tpwr{pid(id)};
        pid = setdiff(pid,pid(id));
    else
        [a, id] = min(aresp(nid));
        nsum = nsum - tpwr{nid(id)};
        nid = setdiff(nid,nid(id));
    end
    pn = length(pid)*nf;
    nn = length(nid)*nf;
    ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
    j = j+1;
    scores(j) = sum(ski(:).*sk(:));
end
details.fixpid = pid;
details.fixnid = nid;


function [pred, details] = RandPkSim(bias, npix, ntrials, noisesd)


for j = 1:length(bias)
    for k = 1:length(noisesd)
tk = rand(npix,1);
noise = rand(2,ntrials) .* noisesd(k);
ims = rand(npix,ntrials);

tresp = tk' * ims;
tresp = [tresp; tresp];
crit = prctile(tresp(:), bias(j) * 100);
choices = tresp >= crit;

sresp = tresp + noise;
crit = prctile(sresp(:), bias(j) * 100);
ichoices = sresp >= crit;
pred(j,k) = sum(choices(:) == ichoices(:))./length(ichoices(:));
details.consistency(j,k) = sum(ichoices(1,:) == ichoices(2,:))/length(choices);
details.biasc(j) = bias(j).^2 + (1-bias(j)).^2;
q = pred(j,k) - details.biasc(j);
p = bias(j);
np =1-p;
details.predconsist(j,k) = p^2 + (1-p)^2 + 2 * q^2 + 2*q*(p^2+np^2)- 4 * q * p * np ;
details.predconsist(j,k) = 0.5 + 2*q^2;
    end
end

if length(bias) > 2 & length(noise) > 2
    plot(pred',details.consistency');
    refline(1);
    hold on;
    plot(pred',details.predconsist','--');
end