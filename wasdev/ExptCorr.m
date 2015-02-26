function [exptxc, details] = ExptCorr(exa,exb,conditions,varargin)
%[exptxc, details] = ExptCorr(exa,exb,conditions,varargin)
%Calculates cross-correlation between spikes in exa and spikes in xb;
usecluster = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'hash',4)
        usecluster = 0;
    end
    j = j+1;
end

ta = exa.Trials;
tb = exb.Trials;
times = -3010:10:3010;
trials = intersect([ta.Trial],[tb.Trial]);

ai = find(ismember([ta.Trial],trials));
bi = find(ismember([tb.Trial],trials));
if isempty(bi) | isempty(ai)
    [a,b] = fileparts(exa.Header.Name);
    fprintf('No overlapping Trials %s.%s cells %d,%d\n',b,exa.Header.expname,exa.Header.cellnumber,exb.Header.cellnumber);
    exptxc = NaN;
    details.rawxc = NaN;
    details.times = [];
    details.rsc = [];
    details.counts = [];
    details.ntrials = 0;
    details.state = [];
    details.sxc = NaN;
    return;
end

for j = 1:length(ai)
    durs(j) = ta(ai(j)).End(end) - ta(ai(j)).Start(1);
end
dur = median(durs);
ptimes = 0:10:mean(durs);

if isempty(conditions)
    aid = ai;
    bid = bi;
end
for c = 1:length(conditions)
ti = eval(['find([ta(ai).' conditions{c} ')']);
aid = ai(ti);
bid = bi(ti);
end
rptshuffle = 0;
if isempty(bid) | isempty(aid)
    fprintf('No Trials satisfy %s, %s\n',sprintf('%s,',conditions{:}),exa.Header.expname);
    exptxc = NaN;
    details = [];
    return;
end
if isfield(ta,'se')
    ses = [ta(aid).se];
    [secounts, sess] = Counts(ses);
    id = find(secounts > 1);
    if length(id) > length(aid)/3;
        id = find(ismember(ses, sess(id)));
        aid = aid(id);
        bid = bid(id);
        rptshuffle = 1;
    end
end
xc = zeros(size(times));
if diff(size(xc)) < 0
    xc =xc';
end
sxc = xc;
details.trialids = [ta(aid).id];
if isfield(tb,'se')
ses = [tb(bid).se];
end
[d, zbin] = min(abs(times));
rnd = randperm(length(aid));
for j = 1:length(aid)
    if usecluster
    aspk = ta(aid(j)).Spikes;
    bspk = tb(bid(j)).Spikes;
    else
    aspk = ta(aid(j)).OSpikes;
    bspk = tb(bid(j)).OSpikes;
    end
    acounts(j) = sum(aspk > 500 & aspk < dur+500);
    bcounts(j) = sum(bspk > 500 & bspk < dur+500);
    if rptshuffle
        id = find(ses == ses(j) & [1:length(ses)] ~= j);
        if usecluster
        cspk = tb(bid(id(1))).Spikes;
        else
        cspk = tb(bid(id(1))).OSpikes;
        end
    else
        cspk = tb(bid(rnd(j))).Spikes;
    end
        ccounts(j) = sum(cspk > 500 & cspk < dur+500);
    scounts = zeros(2,length(aspk));
    for k = 1:length(aspk)
        cc = histc(aspk(k) - bspk',times); %transpose works if length bsk == 1
        scounts(1,k) = cc(zbin);
        xc = xc + cc;
        cc = histc(aspk(k) - cspk',times);
        sxc = sxc + cc;
        scounts(2,k) = cc(zbin);
    end
    synccounts(:,j) = sum(scounts,2);
    psth{1}(j,:) = histc(aspk,ptimes);
end
details.rawxc = xc(2:end-1);
xc = xc(2:end-1) - sxc(2:end-1);
exptxc = xc;
details.times = times(2:end-1);
%ccounts are shuffled. So rsc(2) and rsc(4) are shuffled controls
rsc = corrcoef(acounts,bcounts);
details.rsc(1) = rsc(1,2);
rsc = corrcoef(acounts,ccounts);
details.rsc(2) = rsc(1,2);
details.counts = cat(1,acounts,bcounts);
rsc = corrcoef(acounts-synccounts(1,:),bcounts-synccounts(1,:));
details.rsc(3) = rsc(1,2);
rsc = corrcoef(acounts-synccounts(2,:),ccounts-synccounts(2,:));
details.rsc(4) = rsc(1,2);
if isfield(exa,'probes')
    details.probes = [mean(exa.probes(aid)) mean(exb.probes(bid))];
elseif isfield(exa.Header,'probe')
    details.probes = [exa.Header.probe exb.Header.probe];
else
    details.probes = [0 0];
end
ps = GetEval(exa,'probesep');
if ps > 0
     details.probesep = diff(details.probes) * ps;
else
    details.probesep = NaN;
end
details.ntrials = [length(aid) length(ai)];
details.state.rpts = rptshuffle;
details.state.means = [mean(acounts) mean(bcounts)];
details.sxc = sxc(2:end-1);
xc = zeros(size(times));
if diff(size(xc)) < 0
    xc =xc';
end


if isfield(ta,'RespDir')
    pid = find([ta.RespDir] > 0);
    nid = find([ta.RespDir] < 0);
    apid = intersect(pid,aid);
    anid = intersect(nid,aid);
    pid = find([tb.RespDir] > 0);
    nid = find([tb.RespDir] < 0);
    bpid = intersect(pid,bid);
    bnid = intersect(nid,bid);
    for j = 1:length(apid)
        if usecluster
            aspk = ta(apid(j)).Spikes;
            bspk = tb(bpid(j)).Spikes;
        else
            aspk = ta(apid(j)).OSpikes;
            bspk = tb(bpid(j)).OSpikes;
        end
        acounts(j) = sum(aspk > 500 & aspk < dur+500);
        bcounts(j) = sum(bspk > 500 & bspk < dur+500);
    scounts = zeros(2,length(aspk));
    for k = 1:length(aspk)
        cc = histc(aspk(k) - bspk',times);
        scounts(1,k) = cc(zbin);
        xc = xc + cc;
        cc = histc(aspk(k) - cspk',times);
        sxc = sxc + cc;
        scounts(2,k) = cc(zbin);
    end
    end
    synccounts(:,j) = sum(scounts,2);
    psth{1}(j,:) = histc(aspk,ptimes);
    details.state.choicexc(1,:) = xc(2:end-1);

    for j = 1:length(anid)
        if usecluster
            aspk = ta(anid(j)).Spikes;
            bspk = tb(bnid(j)).Spikes;
        else
            aspk = ta(anid(j)).OSpikes;
            bspk = tb(bnid(j)).OSpikes;
        end
        acounts(j) = sum(aspk > 500 & aspk < dur+500);
        bcounts(j) = sum(bspk > 500 & bspk < dur+500);
    scounts = zeros(2,length(aspk));
    for k = 1:length(aspk)
        cc = histc(aspk(k) - bspk',times);
        scounts(1,k) = cc(zbin);
        xc = xc + cc;
        cc = histc(aspk(k) - cspk',times);
        sxc = sxc + cc;
        scounts(2,k) = cc(zbin);
    end
    end
    synccounts(:,j) = sum(scounts,2);
    psth{1}(j,:) = histc(aspk,ptimes);
    details.state.choicexc(2,:) = xc(2:end-1);
    details.state.choicen = [length(apid) length(anid)];
end
