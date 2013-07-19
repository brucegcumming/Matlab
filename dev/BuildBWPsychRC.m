function [kernel, details] = BuildBWPsychRC(Expt, varargin)
% [kernel, details] = BuildBWPsychRC(Expt, varargin)
% Builds Psych and STA kernels in Fourier domain.
% [kernel, details] = BuildBWPsychRC(Expt, 'sprc') to build the STA
% [kernel, details] = BuildBWPsychRC(Expt) to build psych only
%BuildBWPsychRC(details,kernel) to plot a previously calcualted result


rotate = 0;
nbins = 6;
seedoffset = 815; %works for lemM079
sprc = 0;
delays = [0 300:50:1000];
permute = 0;
ftpwr = [];
fts = [];
framerate = 166.7;
bw = 130;
or = [];
meanimage = 0;
pargs = {}; pj = 0;
plottype = 'BWRC';

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'bw',2);
        j = j+1;
        bw = varargin{j};
    elseif strncmpi(varargin{j},'delay',4);
        j = j+1;
        delays = varargin{j};
    elseif strncmpi(varargin{j},'mean',4);
        meanimage = 1;
    elseif strncmpi(varargin{j},'nojack',4);
        pj = pj+1;
        pargs{pj} = varargin{j};
    elseif strncmpi(varargin{j},'or',2);
        j = j+1;
        or = varargin{j};
    elseif strncmpi(varargin{j},'predict',4)
        j = j+1;
        kernel = varargin{j};
        j = j+1;
        details = varargin{j};
        if isempty(ftpwr)
            ftpwr = abs(details.fts);
        end
        pid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] == 130);
        nid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] == 130);
        zid = [pid nid];
        ses = [Expt.Trials(zid).Seedseq];
        choices = [Expt.Trials(zid).RespDir];
        details = PredictChoice(kernel, ftpwr, ses, choices, details, varargin{j+1:end});
        ppid = find(details.pkresp >= details.criterion(1));
        pnid = find(details.pkresp < details.criterion(1));
        details.ccp(4) = CalcCP([Expt.Trials(ppid).count],[Expt.Trials(pnid).count]);
        if isfield(details,'pkaresp') %used jacknife
        ppid = find(details.pkaresp >= details.criterion(2));
        pnid = find(details.pkaresp < details.criterion(2));
        details.ccp(5) = CalcCP([Expt.Trials(ppid).count],[Expt.Trials(pnid).count]);
        end
        return;
    elseif strncmpi(varargin{j},'permute',4)
        permute = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            permute = varargin{j};
        end
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'rotate',4);
        j = j+1;
        rotate = varargin{j} .* pi/180;
    elseif strncmpi(varargin{j},'sprc',4);
            sprc = 1;
    elseif strncmpi(varargin{j},'ftpwr',3);
        j = j+1;
        ftpwr = varargin{j};
    elseif strncmpi(varargin{j},'fts',3);
        j = j+1;
        fts = varargin{j};
        ftpwr = abs(fts);
    end
    j = j+1;
end

if isfield(Expt,'predchoice') & isnumeric(varargin{1}) %its a result file.
    [kernel, details] = PlotKernel(varargin{1}, Expt, plottype, varargin{2:end});
    return;
end
px  = 0.0292;
sf = Expt.Stimvals.sf;
wi = Expt.Stimvals.wi;
if isempty(or)
or = Expt.Stimvals.or;
end
if isfield(Expt.Trials,'seedoffset')
    id = find([Expt.Trials.ob] > 120);
    seedoffset = median([Expt.Trials.ob]);
end
for j = 1:length(Expt.Trials)
    lens(j) = length(Expt.Trials(j).Seedseq);
end
gid = lens >= prctile(lens,90);
badid = find(lens < prctile(lens,90));
for j = 1:length(badid)
    fprintf('Trials %d only %d seeds\n',badid(j),lens(badid(j)));
end
pid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] == 130 & gid);
nid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] == 130 &gid);

if length(pid) == 0 | length(nid) == 0
    kernel = NaN;
    return;
end
zid = [pid nid];
choices = [ones(size(pid)) ones(size(nid)).*-1];
if isfield(Expt.Trials,'imseed')
    seedoffset = median([Expt.Trials(zid).imseed]);
end

if isempty(ftpwr)
fts = filterim([sf sf/2], [or bw], wi, 'seedoffset', seedoffset, 'pix2deg', px,'nseeds', 1000, 'getft','noplot');
pwr = abs(fts);
else
    pwr = ftpwr;
end

if meanimage
    kernel = fftshift(squeeze(mean(pwr)));
    a = kernel'; %gives vertial blob for vertical ori
    imagesc(fliplr(a));
    aor = AddOriSum(a);
    [a,b] = GetAngle(aor);
    title(sprintf('OR= %.0f->%.0f (%.0f), bw = %.1f',or,a,b,bw));
    details.meanor = aor;
    return;
end

ses = [Expt.Trials(zid).Seedseq];
sevals = unique(ses);

for p = permute:-1:1
    id = randperm(length(zid));
    ipid = zid(id(1:length(pid)));
    inid = zid(id(length(pid)+1:end));
    ichoices = ones(size(zid));
    ichoices(id(length(pid)+1:end)) = -1;
    [ki, a] = CalcPKernel(Expt,ipid, inid, pwr);
    k(p,:,:) = ki;
    a = PredictChoice(squeeze(k(p,:,:)),pwr,ses, ichoices, a);
    pperf(p) = a.predchoice(1);
    pvar(p) = var((ki(:)));
end

if permute
end

[kernel, details] = CalcPKernel(Expt,pid, nid, pwr);

details.kb = CalcPKernel(Expt,pid(2:end), nid, pwr);

if permute
    details.kernelsd = squeeze(std(k,[],1));
    details.permpred = pperf;
    details.permpkvar = pvar;
end


for j = 1:size(ses,2)
alltriggers{j} = [];
end
for j = 1:size(ses,2)
    [ij,ik] = find(ses == sevals(j));
    for k = 1:length(ik)
        alltriggers{ik(k)}= {alltriggers{ik(k)} ij(k)};
    end
end

details = PredictChoice(kernel, pwr, ses, choices, details,pargs{:});
details.pkor = AddOriSum(kernel');

if sprc
    [secounts, sevals] = Counts(ses);
    sxmean = zeros(size(fts,2),size(fts,3));
    for j = 1:length(sevals);
        sxmean = sxmean + squeeze(fts(sevals(j),:,:)) .* secounts(j);
    end
    details.allsxmean = sxmean./sum(secounts);
    details.delays = delays;
    for di = 1:length(delays)
        delay = delays(di);
        spseeds = [];
        nf = size(ses,1);
        k = 1:length(zid);
        for j = 1:length(zid)
            frames = floor((Expt.Trials(zid(j)).Spikes-delay)./166);
            frames = frames(find(frames >0 & frames <= nf));
            spseeds = [spseeds ses(frames,k(j))'];
        end

        spsum = zeros(size(pwr,2),size(pwr,3));
        spsxsum = zeros(size(pwr,2),size(pwr,3));
        [nspk, sv] = Counts(spseeds);
        for j = 1:length(sv)
            se = sv(j);
            spsum = spsum + nspk(j) .* squeeze(pwr(se,:,:));
            spsxsum = spsxsum + nspk(j) .* squeeze(fts(se,:,:));
        end
        spmean = spsum./sum(nspk) - details.allmean;
        details.spmean(di,:,:) = fftshift(spmean);
        details.spvar(di) = var(spmean(:));
        sxmean = spsxsum./sum(nspk) - details.allsxmean;
        details.sxmean(di,:,:) = fftshift(sxmean);
        details.sxvar(di) = var(abs(sxmean(:)));
        nspks(di,sv) = nspk;
    end
    [v,sxt] = max(details.sxvar);
    [v,t] = max(details.spvar);
    nspk = nspks(t,:);
    spk = fftshift(details.spmean(t,:,:)); %%undo fftshift above
    sxspk = fftshift(details.sxmean(t,:,:)); %%undo fftshift above
    for j = 1:length(sevals)
        id = find(sv == sevals(j));
        sid = find(sevals == j);
        if isempty(id) | isempty(sid)
            ps(sevals(j)) = 0;
        else
            ps(sevals(j)) = nspk(sevals(j))./secounts(j);
            pss(sevals(j),:) = nspks(:,sevals(j))./secounts(j);
        end
        fpwr(sevals(j)) = sum(sum(pwr(sevals(j),:,:) .* spk));
        sxpwr(sevals(j)) = abs(sum(sum(fts(sevals(j),:,:) .* sxspk)));
    end
    sevars = var(pss);
    lo = mean(abs(squeeze(pwr(ps<prctile(ps,30),:,:))));
    hi = mean(abs(squeeze(pwr(ps>prctile(ps,70),:,:))));
    for j = 1:length(delays)
        p = pss(:,j);
        lo = mean(abs(squeeze(pwr(p<prctile(p,30),:,:))));
        hi = mean(abs(squeeze(pwr(p>prctile(p,70),:,:))));
        psmean(j,:,:) = fftshift(squeeze(hi-lo));
    end
    [a,b] = sort(fpwr(sevals));
    k = 1;
    nf = 0;
%if some seeds are missing, than length(fwr) = largest seed, not # of
%seeds
    bincount = sum(secounts)./nbins;
    for j = 1:length(sevals) 
        porder(sevals(b(j))) = k;
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(porder == k) < 100;
        porder(porder == k) = k-1;
    end
    [a,b] = sort(ps(sevals));
    k = 1; nf = 0;
    for j = 1:length(sevals) 
        forder(sevals(b(j))) = k;
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(forder == k) < 100;
        forder(forder == k) = k-1;
    end
    [a,b] = sort(sxpwr(sevals));
    k = 1; nf = 0;
    for j = 1:length(sevals) 
        sorder(sevals(b(j))) = k;
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(sorder == k) < 100;
        forder(sorder == k) = k-1;
    end
    subplot(2,2,1);
    hold off;
    plot(fpwr(sevals),ps(sevals),'o');
    pfit = polyfit(fpwr(sevals),ps(sevals),2);
    netspk(sevals) = polyval(pfit,fpwr(sevals));
    hold on;
    [a,b] = sort(fpwr(sevals));
    plot(a,netspk(b),'linewidth',2);
    [x,y] = xysmooth(fpwr(sevals),ps(sevals),100);
    plot(x,y,'r');
    [x,y] = xysmooth(fpwr(sevals),ps(sevals),10);
    plot(x,y,'g');
    latency = 500;
    for j = zid
        Expt.Trials(j).Predcount(1) = sum(netspk(Expt.Trials(j).Seedseq));
        Expt.Trials(j).Predcount(2) = sum(ps(Expt.Trials(j).Seedseq));
        Expt.Trials(j).count = sum(Expt.Trials(j).Spikes > latency & Expt.Trials(j).Spikes < Expt.Trials(j).dur + latency);
        Expt.Trials(j).sepower = porder(Expt.Trials(j).Seedseq);
        Expt.Trials(j).frpower = forder(Expt.Trials(j).Seedseq);
        Expt.Trials(j).sxpower = sorder(Expt.Trials(j).Seedseq);
    end
    cfig = gcf;
    details.nspks = sum([Expt.Trials(zid).count]);
    if isfield(Expt.Header,'frameperiod')
        frameperiod = Expt.Header.frameperiod;
    else
        frameperiod = 166.6;
    end
    details.rc = PlotRevCorAny(Expt,'exp','sepower','Trials',zid,'filltimes',frameperiod,'nmin',100,'box','sdfw',framerate);
    netspk = PlotRC(details.rc,'netspk','sdsmooth',5,'timerange',[200 1000]);
    details.predcount = [Expt.Trials(zid).Predcount];
    for j = 1:max(porder)
        id = find(porder == j);
        netspkx(j) = mean(fpwr(id));
    end
    figure(cfig);
    plot(netspkx,netspk,'ro-');
    for j = zid
        Expt.Trials(j).Predcount(3) = sum(netspk(Expt.Trials(j).sepower));
    end
    subplot(2,2,2);
    pc = cat(1,Expt.Trials.Predcount)';
    hold off;
    plot(pc(1,:),[Expt.Trials(zid).count],'o');
    hold on;
    plot(pc(2,:),[Expt.Trials(zid).count],'go');
    plot(pc(3,:),[Expt.Trials(zid).count],'ro');
    subplot(2,2,3);
    imagesc(fliplr(squeeze(details.spmean(t,:,:))'));
    xc = corrcoef(kernel(:), details.spmean(t,:,:));
    details.psxcorr(1) = xc(2,1);
    [details.spkor, spring] = AddOriSum(squeeze(details.spmean(t,:,:))');
    xc = corrcoef(details.pkor, details.spkor);
    details.psxcorr(2) = xc(2,1);
    [a,b,c] = GetAngle(details.spkor);
    title(sprintf('VarRatio %.2f Pkxc %.3f,%.3f Or %.0f',details.spvar(t)./details.spvar(1),details.psxcorr(1),details.psxcorr(2),a)); 
    details.cp = CalcCP([Expt.Trials(pid).count],[Expt.Trials(nid).count]);
    details.ccp(1) = CalcCP([Expt.Trials(pid).count]-pc(1,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(1,ismember(zid,nid)));
    details.ccp(2) = CalcCP([Expt.Trials(pid).count]-pc(3,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(3,ismember(zid,nid)));
    details.ccp(4) = CalcCP([Expt.Trials(details.fixpid).count],[Expt.Trials(details.fixnid).count]);
    details.ccp(3) = CalcCP([Expt.Trials(pid).count]-pc(2,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(2,ismember(zid,nid)));
    ppid = find(details.pkresp >= details.criterion(1));
    pnid = find(details.pkresp < details.criterion(1));
    details.ccp(5) = CalcCP([Expt.Trials(zid(ppid)).count],[Expt.Trials(zid(pnid)).count]);
    
    details.ccp(6) = CalcCP(pc(1,ppid),pc(1,pnid));
    details.ccp(7) = CalcCP(pc(3,ppid),pc(3,pnid));
    details.ccp(8) = CalcCP(pc(2,ppid),pc(2,pnid));
    if isfield(details,'pkaresp') %used jacknife
        ppid = find(details.pkaresp >= details.criterion(2));
        pnid = find(details.pkaresp < details.criterion(2));
        details.ccp(9) = CalcCP([Expt.Trials(ppid).count],[Expt.Trials(pnid).count]);
    end
    subplot(2,2,4);
else
    spring.r = 0;
end
    



if rotate
    [xi, yi] = meshgrid([1:wi],[1:wi]);
    xr = (xi-128) .* cos(pi/4) + (yi-128) .* sin(pi/4);
    yr = (yi-128) .* cos(pi/20) - (xi-128) .* sin(pi/20);
end
imagesc(fliplr(kernel'));
[details.pkor, b] = AddOriSum(kernel');
if b.r < spring.r
    imagesc(fliplr(kernel'));
    AddOriSum(kernel','maxr',spring.r);
end
[a,b,c] = GetAngle(details.pkor);
if isfield(details,'cp')
title(sprintf('CP %.2f, %.2f, %.2f Or%.0f',details.cp,details.ccp(1),details.ccp(2),a));
end
details.pkvar = var(kernel(:));
details.fts = fts;

function [kernel, details] = CalcPKernel(Expt, pid, nid, pwr)



[a,b] = Counts([Expt.Trials(pid).Seedseq]);
[c,d] = Counts([Expt.Trials(nid).Seedseq]);

psum = zeros(size(pwr,2),size(pwr,3));
nsum = zeros(size(pwr,2),size(pwr,3));
for j = 1:length(a)
    se = b(j);
    psum = psum + squeeze(pwr(se,:,:)) .* a(j);
end
for j = 1:length(c)
    nsum = nsum + squeeze(pwr(d(j),:,:)) .* c(j);
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


function [kernel, details] = PlotKernel(kernel, details, plottype, varargin)


j = 1; 
while j <= length(varargin)
    if strncmpi(varargin{j},'var',3)
    end
    j = j+1;
end


ntrials = length(details.predcount);

if strncmpi(plottype,'RCVar',4)
    subplot(2,2,1);
    hold off;
    [a,t] = max(details.spvar);
    sk = squeeze(details.spmean(t,:,:))';
    imagesc(fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('PWR %d trials. Or %.0f, peak %.0f r%.2f',ntrials,a,b,c));
    subplot(2,2,2);
    hold off;
    plot(details.delays./10,details.spvar./details.spvar(1));
    hold on;
    plot(details.delays./10,details.sxvar./details.sxvar(1),'r');
    legend('pwr','sta');
    subplot(2,2,3);
    hold off;
    [a,t] = max(details.sxvar);
    sk = squeeze(abs(details.sxmean(t,:,:)))';
    imagesc(fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('STA. Or %.0f, peak %.0f r%.2f',a,b,c));
    subplot(2,2,4);
    hold off;
    imagesc(ifft2(squeeze(details.sxmean(t,:,:))));
    
elseif strncmpi(plottype,'BWRCall',4)
    subplot(2,3,1);
    hold off;
    [a,t] = max(details.spvar);
    sk = squeeze(details.spmean(t,:,:))';
    imagesc(fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('PWR %d trials. Or %.0f, peak %.0f r%.2f',ntrials,a,b,c));
    subplot(2,3,4);
    hold off;
    plot(details.delays./10,details.spvar./details.spvar(1));
    hold on;
    plot(details.delays./10,details.sxvar./details.sxvar(1),'r');
    legend('pwr','sta');

    subplot(2,3,5);
    hold off;
    [a,t] = max(details.sxvar);
    sk = squeeze(abs(details.sxmean(t,:,:)))';
    imagesc(fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('STA. Or %.0f, peak %.0f r%.2f',a,b,c));
    
    subplot(2,3,6);
    hold off;
    imagesc(ifft2(squeeze(details.sxmean(t,:,:))));
    subplot(2,3,2);
    hold off;
    imagesc(fliplr(kernel'));
    [a,pb] = AddOriSum(kernel');
    [a,t] = max(details.spvar);
    [a,b,c] = GetAngle(details.pkor);
    title(sprintf('Pk: var %.2f, xc %.3f,%.3f Or%.0f',details.pkvar,details.psxcorr(1),details.psxcorr(2),a));
    details.plotr = max([pb.r sb.r]);
else
    subplot(2,1,1);
    hold off;
    imagesc(fliplr(kernel'));
    [a,pb] = AddOriSum(kernel');
    [a,t] = max(details.spvar);
    subplot(2,1,2);
    hold off;
    [a,t] = max(details.spvar);
    sk = squeeze(details.spmean(t,:,:))';
    imagesc(fliplr(sk));
    [a,sb] = AddOriSum(sk);
    if sb.r < pb.r
        subplot(2,1,2);
        hold off;
        imagesc(fliplr(sk));
        [a,sb] = AddOriSum(sk,'maxr',pb.r);
        details.plotr = pb.r;
    else
        subplot(2,1,1);
        hold off;
        imagesc(fliplr(kernel)');
        [a,pb] = AddOriSum(kernel','maxr',sb.r);
        details.plotr = sb.r;
    end
    subplot(2,1,1);
    [a,b,c] = GetAngle(details.pkor);
    title(sprintf('Pk: var %.2f, xc %.3f,%.3f Or%.0f',details.pkvar,details.psxcorr(1),details.psxcorr(2),a));
    subplot(2,1,2);
    [a,b,c] = GetAngle(details.spkor);
    title(sprintf('VarRatio %.2f, CP %.2f, %.2f, %.2f  Or%.0f,%.0f',details.spvar(t)./details.spvar(1),...
        details.cp,details.ccp(1),details.ccp(2),a,b));
end
    

function [veco, maxo, r] = GetAngle(ors)
%assumes ors is [0:1180] degrees, as in AddOriSum
    [a,maxo] = max(ors(2:end));
    R = ors(2:end) - min(ors);
    %double the angles - its orientation not diretion. 
    sa = sin([pi/90:pi/90:2*pi]);
    ca = cos([pi/90:pi/90:2*pi]);
    r = sum(R .* (ca + i * sa))./sum(R);
    veco = angle(r).* 90/pi; %halve
    r = abs(r);

function [sor, details]  = AddOriSum(kernel, varargin)
coss = cos([0:pi/180:pi]);
sins = sin([0:pi/180:pi]);


fixr = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'maxr',4)
        j = j+1;
        fixr = varargin{j};
    end
    j = j+1;
end
for r = 0.5:45.5
    xi = r .* coss;
    yi = r .* sins;
    sir(ceil(r),:) = interp2(kernel, 129+xi, 129+yi);
end
svar = var(sir');
id = find(svar > max(svar)./50);
if fixr
    r= fixr;
else
    r = id(end);
end
details.r = r;

xi = 129 + [r .* coss; (r-1) .* coss];
yi = 129 + [r .* sins; (r-1) .* sins];
uxi = 129 - [r .* coss; (r-1) .* coss];
uyi = 129 - [r .* sins; (r-1) .* sins];
cr = caxis;
sor = sum(sir);
sor = (range(cr) .* (sor - min(sor))./range(sor)) + cr(1);
set(gca,'xlim',[129-r 129+r],'ylim',[129-r 129+r]);
hold on;
pcolor(xi,yi,fliplr([sor; sor]));
pcolor(uxi,uyi,fliplr([sor; sor]));
[a, b] = max(sor);
[a, c] = min(sor);
plot([uxi(1,b) xi(1,b)],[yi(1,b) uyi(1,b)],'r');
plot([uxi(1,c) xi(1,c)],[yi(1,c) uyi(1,c)],'b');
shading('flat');




function details = PredictChoice(kernel, pwr, seeds, choices, details, varargin)

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
for j = 1:size(pwr,1)
    fpwr(j) = sum(sum(squeeze(pwr(j,:,:)).*sk));
end
resp = sum(fpwr(seeds));
aresp = resp; %copy so resp can be used for jacknife
if jacknife
%   fsums = sum(pwr(seeds,:,:)); Ntrials * 200 x 256 x 256 = too much
%   memory
    for j = 1:length(choices)
        fsum = squeeze(sum((pwr(seeds(:,j),:,:))));
        if choices(j) == 1
        psum = details.psum - fsum;
        pn = details.nframes(1) - nf;
        nsum = details.nsum;
        nn = details.nframes(2);
        else
        nsum = details.nsum - fsum;
        nn = details.nframes(2)-nf;
        psum = details.psum;
        pn = details.nframes(1);
        end
        ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
        resp(j) = sum(sum(fsum.*ski));
    end
else
%calculate projection of each frame onto kernel

end
% set a low criterion to get more + choices
p = sum(choices==-1)./length(choices);
crit = prctile(resp,p.*100);
details.criterion(1) = crit;
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
        fsum = squeeze(sum((pwr(seeds(:,pid(id)),:,:))));
        pid = setdiff(pid,pid(id));
        psum = psum - fsum;
    else
        [a, id] = min(aresp(nid));
        fsum = squeeze(sum((pwr(seeds(:,nid(id)),:,:))));
        nsum = nsum - fsum;
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
