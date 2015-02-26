function res = PlotSpikeC(dat, th, varargin)
%Build Structure with one large matrix of continuous voltages, from
%componend files:
%
% res = PlotSpikeC(files, th, 'probe',[1:24], 'sumv', 'makev', 
%                                  ... ,'submean') removes mean (across all probes) from each.         
%plot spikes from a continuous voltage record.
%
% PlotSpikeC(a, th,......)
% if th < -0, finds all local minima < th
%
%for lem58 use
%a = load('/bgc/data/tests/lemM158.spkblk8.mat');
%b = load('/bgc/data/tests/lemM158A.spkblk39.mat')
%res = PlotSpikeC(a,5,b,'probe',[1:24],'submean','sumv','makev');
%AllVPcs(res.V,'tchan',[3],'nprobepc',[1],'nspk',8000,'spts',-10:20,'tryall
%','smoothsd',2); 

spk = 13;
clplot = 0;
j = 1;
if iscell(dat)
    C = dat;
elseif isstruct(dat)
    C = dat;
else
C{1} = dat;
end
clear dat;
sumreject = 5;
sumv = 0; 
fullcov = 0;
nf = 1;
prct = 0;
rmcluster = 0;
testsplit = 0;
autocut = 0;
rwtmplt = [0.183803 0.227398 0.333533 0.494481 0.644137 0.843963 1.180304 1.697385 2.256425 2.645227 2.743000 2.589913 2.315911 2.089542 1.961463 1.861570 1.736660 1.598173 1.473007 1.368351 1.277644 1.196907 1.125743 1.061659 1.000064 0.938558 0.878896 0.823377 0.774765 0.734389 0.698591 0.662071 0.621361 0.577299 0.534229 0.493684 0.453253 0.413407 0.379653 0.355637 0.339941 0.325690 0.304629 0.275406 0.243394 0.213982 0.190089 0.171977 0.158428 0.146101 0.128928 0.105108 0.081825 0.065891 0.054421 0.040305 0.022889 0.006364 -0.008726 -0.022893 -0.033612];
checkmains = 0;
plotv = 3;
clusterprops = [];
tmains = [];
checkrw = 0;
showrwv = 1;
nspk = 0;  %number of events to find
makeV = 0;
prefix = [];
lastblk = 0;
checkblocks = 1;
useall = 0;
Expt = [];

while j <= length(varargin)
    if isstruct(varargin{j})
       nf = nf+1;
       C{nf} = varargin{j};
    elseif strncmpi(varargin{j},'centile',4)
        j = j+1;
        prct = varargin{j};
    elseif strncmpi(varargin{j},'cluster',4)
        autocut = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j =j+1;
            clusterprops = varargin{j};
        end
    elseif strncmpi(varargin{j},'Expt',4)
        j = j+1;
        Expt = varargin{j};
    elseif strncmpi(varargin{j},'submean',7)
        checkrw = 2;
    elseif strncmpi(varargin{j},'checkrw',7)
        checkrw = 1;
        j = j+1;
        Ev = varargin{j};
        npts = 20000;
    elseif strncmpi(varargin{j},'density',4)
        clplot = 1;
    elseif strncmpi(varargin{j},'fullcov',4)
        fullcov = 1;
    elseif strncmpi(varargin{j},'lastblk',5)
        if strncmpi(varargin{j},'lastblktest',9)
            testsplit = 1;
        end
        j = j+1;
        lastblk = varargin{j};
    elseif strncmpi(varargin{j},'makev',5)
        makeV =1;
    elseif strncmpi(varargin{j},'mains',5)
        checkmains = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j =j+1;
            tmains = varargin{j};
        end
    elseif strncmpi(varargin{j},'nov',3)
        plotv = 0;
    elseif strncmpi(varargin{j},'quickv',6)
        plotv = 2;
    elseif strncmpi(varargin{j},'prefix',4)
        j = j+1;
        prefix = varargin{j};
    elseif strncmpi(varargin{j},'probe',4)
        j = j+1;
        spk = varargin{j};
    elseif strncmpi(varargin{j},'rmcluster',4)
        rmcluster = 1; %redo PCA after removing strongest cluster
    elseif strncmpi(varargin{j},'sumv',4)
        sumv = 1; %redo PCA after removing strongest cluster
    end
    j = j+1;
end
sds(1) = 0;
sds(16) = 1;
samper=1./40000;
avg = [];
maxpts = 200000000; %made v big July 2012
res.errmsg = {};
res.errdata = [];

if length(spk) > 1 && sumv
    if isempty(C)
        return;
    end

    if isstruct(C) %given a list of "probes" structs
    [chspk, lens, starts, offsets, f, X] = CountProbeChans(C, spk);
    else
    [chspk, lens, starts, offsets, f, X] = CountChans(C, spk);
    end
     if isfield(X,'errs')
         res.errmsg = {res.errmsg{:}, X.errs{:}};
     end
     if isfield(X,'errdata')
         res.errdata = {res.errdata X.errdata};
     end
    if useall
            np = max(chspk);
            V(np,npts) = 0; % preallocate
            res.samper = X.samper;
    end
    probelist = unique(chspk);
    id = find(X.nchans < length(probelist));
    if length(id)
        fprintf('Some Blocks are missing channels:  %s\n', sprintf(' %d',id));
%These blocks are NOT included in the FullV
        sumt = 0;
        badid = id;
        badblk = [];
        for j = 1:length(id)
            bid = find(X.blks == id(j));
            msg = sprintf('%s %d Chans at %.1f-%.1f. Excluded Time Range',X.names{X.ix(bid(1))},...
                X.nchans(id(j)),X.blkstart(id(j)),X.blkstart(id(j))+X.blklen(id(j)).*X.samper);
            fprintf('%s\n',msg);
            res.errmsg = {res.errmsg{:} msg};
            sumt = sumt + X.blklen(id(j)).*X.samper;
            badblk = [badblk bid];
        end
        res.error = 0;
        res.badblk = badblk;
        res.nchans = X.nchans;
        res.probelist = probelist;
        res.sumerrtime = sumt;
        if sumt > 10
            res.error = 1;
            fprintf('%.1f Seconds of missing channels - Aborting write\n',sumt);
            return;
        end
        id = setdiff(1:length(X.blks),badblk);
        X.blks = X.blks(id);
        X.ix = X.ix(id);
        X.pid = X.pid(id);
        X.totaloffsets = X.totaloffsets(id);
        f = f(id);
        offsets = offsets(id);
        lens = lens(id);
        chspk = chspk(id);

        id = setdiff(1:length(X.blklen),badid);
        X.blkstart= X.blkstart(id);
        X.blklen= X.blklen(id);
        blks = unique(X.blks);
        tlen = 1;
        for j = 1:length(blks)
            id = find(X.blks == blks(j));
            X.totaloffsets(id) = tlen;
            tlen = tlen+X.blklen(j);
        end
    end
    npts = sum(X.blklen);
    if lastblk > 0
        nhalf = npts/2;
        allnpts(1) = sum(X.blklen);
        lasti = lastblk+1; %lstblk is the last block in previous file
        allnpts(2) = sum(X.blklen(1:lasti-1));
        X.blklen = X.blklen(lasti:end);
        X.blkstart = X.blkstart(lasti:end);
        id = find(X.blks == lasti);
        o = max(X.totaloffsets(id)); %make sure all offsets are positive
        id = find(X.blks > lastblk);
        X.blks = X.blks(id);
        X.ix = X.ix(id);
        X.totaloffsets = X.totaloffsets(id)-o+1;
        X.pid = X.pid(id);
        starts = starts(id);
        chspk = chspk(id);
        lens = lens(id);
        offsets = offsets(id);
        npts = sum(X.blklen);
        res.firstblk = lastblk;
        fprintf('npts %d -> %d(%.0f) %d(%.0f)\n',allnpts(1),allnpts(2),abs(allnpts(2)-nhalf),npts,abs(npts-nhalf));
        fprintf('npt diffs %.0f %.0f\n',abs(allnpts(2)-nhalf),abs(npts-nhalf));
    elseif npts > maxpts
        fprintf('%d Points - splitting %s\n',npts,X.names{1});
        alln = cumsum(X.blklen);
        if npts/2 < maxpts
        [a, lasti] = min(abs(alln-npts/2)); %nearest to 50/50 split
        else
        id = find(cumsum(X.blklen) < maxpts);
        lasti = id(end)-1;
        end
        X.blklen = X.blklen(1:lasti);
        X.blkstart = X.blkstart(1:lasti);
        id = find(X.blks <= lasti);
        X.blks = X.blks(id);
        X.ix = X.ix(id);
        X.totaloffsets = X.totaloffsets(id);
        starts = starts(id);
        chspk = chspk(id);
        lens = lens(id);
        offsets = offsets(id);
        npts = sum(X.blklen)
        res.lastblk = lasti;
    end
    if testsplit
        return;
    end
    np = max(chspk);
    if isempty(np) || npts == 0
        suffs = unique([C.suffix]);
        res = AddError(res,'No Data available for Expt %d',suffs(1));
        res.error = 1; %won't be saved
        return;
    end
    V(np,npts) = 0; % preallocate
    res.samper = X.samper;
    if isstruct(C) && ~isempty(prefix)
        res.name = [prefix '/' X.names{1}];
    end
    if ~isempty(prefix)
    res.prefix = prefix;
    end
    res.name = X.names{1};
    maxl = min(lens-offsets);
    start = max(starts);
%    [a,b] = sort(chspk);
    fid = unique(X.ix);
    for k = 1:length(fid)
        id = find(X.ix == fid(k));
        if isstruct(C)
            name = [prefix '/' X.names{fid(k)} '.mat'];
            F = load(name);
            for j = id
                s = round(X.totaloffsets(j));
                p = chspk(j);
                ch = C(X.pid(j)).var;
                V(p,s:s+lens(j)-1) = F.(f{j}).values(1+offsets(j):1+offsets(j)+lens(j)-1);
            end
        else
            F = C{fid(k)};
            for j = id
                s = round(X.totaloffsets(j));
                p = chspk(j);
                ch = C(X.pid(j)).var;
                V(p,s:s+lens(j)-1) = F.(f{j}).values(1+offsets(j):1+offsets(j)+lens(j)-1);
            end
        end

        clear F;
    end
    first = 1;
    if makeV
        clear C;
        clear dat;
    end
    for j = 1:length(X.blklen)
        last = first+X.blklen(j)-1;
        vt(first:last) = X.blkstart(j)+[1:X.blklen(j)].*samper;
        first = last+1;
    end
    res.blklen = X.blklen;
    res.blkstart = X.blkstart;
    res.chspk = sort(unique(chspk));
    sumv = mean(V);
    if checkrw == 2
        %mx = sumv * sumv';
        g = V * sumv';
%mean(g) is the same as sumv * sumv' because sumv is the mean of V;        
        g = g./mean(g);
        res.meangain = g;
        for j = 1:length(g)
%            gV(V(j,:) -sumv .*g(j);
        end
        for j = 1:size(V,1) %use less memory
            V(j,:) = V(j,:) - sumv .* g(j);
        end
%        V = V - repmat(sumv,size(V,1),1);
    elseif checkrw
        stdV = std(V);
        rid = find(Ev.codes(:,1) == 113 & Ev.times >= vt(1)-10 & Ev.times <= vt(end));
        rwtimes = cat(1,vt(1),Ev.times(rid));
        smv = smooth(sumv,3,'gauss');
        
        id = find(abs(sumv) > 3); %smoothed mean across channels = saturated
        id = id(diff(id)>3);
        for j = 1:length(id)
            if smv(id(j)) > 0
                pre = find(smv(id(j)-1000:id(j)) <= 0);
                post = find(smv(id(j):id(j)+1000) <= 0);
            else
                pre = find(smv(id(j)-1000:id(j)) >= 0);
                post = find(smv(id(j):id(j)+1000) >= 0);
            end
            sumv(id(j)+[pre(end)-1000:post(1)]) = 0;
            smv(id(j)+[pre(end)-1000:post(1)]) = 0;
            V(:,id(j)+[pre(end)-1000:post(1)]) = 0;
        end
                
        
        
        
        a = diff(smv);
        id = find(diff(sign(diff(a))) < 0 & a(1:end-2) > std(a).*3);
        vid = find(diff(sign(a)) < 0)+1; %peaks
        nid = find(diff(sign(diff(a))) > 0 & a(1:end-2) < std(a).* -3);
        nvid = find(diff(sign(a)) > 0)+1; %troughs
        id = union(id,nid);
        s = a(id);
        rV = zeros(length(id),601);
        rawV = zeros(size(rV));
        for j = length(id):-1:1
            if s(j) < 0
                pid = find(nvid > id(j));
                t = nvid(pid(1));
            else
                pid = find(vid > id(j));
                t = vid(pid(1));
            end
            ti(j) = t;
            d = vt(t)-rwtimes;
            rid = find(d > 0);
            if t > 100 & t < length(smv)-500
            trw(j) = d(rid(end));
            rV(j,:) = smv(t-100:t+500);
            rawV(j,:) = sumv(t-100:t+500);
            elseif t >= length(smv)-500
                n = 101+length(smv)-t;
                trw(j) = d(rid(end));
                rV(j,1:n) = smv(t-100:end);
                rawV(j,1:n) = sumv(t-100:end);
            end
        end
        scores = rwtmplt * rawV(:,90:150)';
        v = (sumv(id+2)+sumv(id+3)+sumv(id+4)+sumv(id+5))/4;
        rid = find(trw > 0.6);
        crit = mean(scores(rid))+std(scores(rid)).*3;
        it = find(scores > crit & trw < 0.6);
        v = rawV(:,100);
        s = min(cat(1,stdV(ti),stdV(ti-1),stdV(ti-2)));
        v = v';
        it = find(abs(v) > 4.*s);
        it = union(it,find(v > 3 .* s & trw < 0.4));
        if showrwv
            DATA.V = V;
            DATA.rti = ti;
            [a, isnew] = GetFigure('RWscatter');
            hold off;
            set(a,'UserData',DATA);
            myscatter(s,v,'.','buttonpressfcn',@PlotV);
            hold on;
            myscatter(s(it),v(it),'r.','buttonpressfcn',@PlotV,'ids',it);
        end
        for j = 1:length(it)
            sumv(id(it(j))-2:id(it(j))+50) = 0;
            V(:,id(it(j))-2:id(it(j))+50) = 0;
        end
        
        
        old = 0;
        if old  %actually check time w.r.t reward events
        id = find(Ev.codes(:,1) == 113 & Ev.times < vt(end) & Ev.times > vt(1));
        for j = 1:length(id)
            ti = find(vt >= Ev.times(id(j)));
            ts(j) = ti(1);
            rV(j,:) = sumv(ti(1):ti(1)+npts);
        end
        aV = rV;
        GetFigure('FullV');
        plot(rV');
        dv = diff(rV(:,1:end-100),1,2);
        vcrit = std(dv(:)).*3;
        [a,b] = find(dv > vcrit);
        na = 1;
        scrit = 2;
        for j = 1:length(a)
            s(j) = sum(rV(a(j),b(j)+2:b(j)+20));
            sa(j) = dv(a(j),b(j));
            if s(j) > scrit
                af(na,:) = rV(a(j),b(j):b(j)+100);
                na = na+1;
                sumv(ts(a(j))+b(j)-1:ts(a(j))+b(j)+50) = 0;
            end
        end
        plot(aV');
        return;
        end
    end
    if checkmains
        period = round(0.01662./samper);
        tv = 1:(period+10);
        d = tmains - start;
        id = find(d > 0);
        t = tmains(id(1));
        tid = find(vt >= t);
        ti = tid(1); % first sampale
        x = ti:ti+tv;
        ts(length(vt)) = 0;
        ts(ti+tv-1) = tv;
        ts(1:ti-1) = tv(end-ti+2:end);
        mv = V(:,ti+tv-1);
        mt(1) = ti;
        j = 2;
        while t < vt(end) && ti+2*period+100 < length(vt)
            tid = find(vt(ti+period-5:ti+period+5) >= tmains(id(j)));
            ti = tid(1)+ti+period-5;
            x = ti:ti+tv;
            ts(ti+tv-1) = tv;
            mv = mv+V(:,ti+tv-1);
            mt(j) = ti;
            j = j+1;
        end
        ts(ti+tv(end):end) = 1;
        mv = mv./j;
        mv = mv-(repmat(mean(mv,2),1,size(mv,2)));
       G = sum(mv(:).*mv(:));
        rmv = zeros(size(V));
    for j = 1:length(mt)-1
        v = V(:,mt(j):mt(j+1)-1);
        g = sum(sum(v .* mv(:,1:length(v))))./G;
        rmv(:,mt(j):mt(j+1)-1) = mv(:,1:length(v));
        gains(j) = g;
    end
        GetFigure('FullV');
        hold off;
    plot(V(1,:));
        hold on;
    plot(rmv(1,:),'r');
    V = V-rmv;
    plot(V(1,:),'g');
    
    end
    for j = size(V,1):-1:1
%        smv(j,:) = smooth(V(j,:),80);
        smv = smooth(V(j,:),80);
        V(j,:) = V(j,:)-smv;
    end

    if makeV
        res.V = V;
        res.start = start;
        res = AddTimesToFullV(res);
        res.meanV = sumv;
        return;
    end
    ispk = 1;
sgn = diff(sign(diff(V,1,2)),1,2);
if th(1) < 0
id = find(sgn(ispk,:) > 0 & V(ispk,2:end-1) < th(1));
else
id = find(sgn(ispk,:) < 0 & V(ispk,2:end-1) > th(1));
end
if isempty(id)
    id = find(sgn > 0)+1;
    prc = 100000./length(id); % get 1000 spikes
    id = id(V(ispk,id) < prctile(V(ispk,id),prc));
end
GetFigure('Spikes');
hold off;
plot([0 40],[0 0]);
hold on;
allid = [];
id = id(id > 8 & id < maxl-32);
vtimes = start+id.*samper;
if size(id,1) == 1
    allid = repmat(id,41,1) + repmat([-8:32]',1,length(id));
else
    allid = repmat(id',41,1) + repmat([-8:32]',1,length(id));
end
AllV = reshape(V(:,allid),[size(V,1) size(allid)]);
fullcov = 1;
%for j = 1:length(id)
%    allid = [allid [id(j)-8:id(j)+32]'];
%end


else

    if checkblocks
        [chspk, lens, starts, offsets, f, X] = CountProbeChans(C, spk);
    end

if length(spk) > 1
    for j = 1:length(spk)
        res{j} = PlotSpikeC(C{1}, th, varargin{:}, 'probe', spk(j));
        drawnow;
        k = res{j}.spk;
        if res{j}.sd > 0.1
            amps(res{j}.spk,:) = res{j}.scale;
        end
        sds(k) = res{j}.sd;
        peaks(k) = min(res{j}.avg(j,:));
    end
    GetFigure('Ampls');
    subplot(2,1,1);
    imagesc(amps);
    subplot(2,1,2);
    dsum = zeros(length(amps),1);
    dn = dsum;
    for j = 1:size(amps,1)
        for k = 1:size(amps,2)
            d = abs(j-k)+1;
            if sds(k) > 0.1
                dsum(d) = dsum(d)+amps(j,k);
                dn(d) = dn(d)+1;
            end
        end
    end
    plot(0:length(dsum)-1,dsum./dn,'o-');
    return;
end

ns = 1;
for p = 1:length(C)
f = fields(C{p});
for j = 1:length(f)
    if isfield(C{p}.(f{j}), 'title') &&strncmp(C{p}.(f{j}).title,'Spike',5)
        chspk(ns) = sscanf(C{p}.(f{j}).title,'Spike %d');
        lens(ns) = length(C{p}.(f{j}).values);
        starts(ns) = C{p}.(f{j}).start;
        if chspk(ns) == spk
            start = starts(ns);
            V = C{p}.(f{j}).values;
            samper =  C{p}.(f{j}).interval;
        end
        ns = ns+1;
    end
end
end
ispk = find(chspk == spk);
res.spkchans = chspk;
%start = min(starts); %start has to match reference channel

if checkmains
    vt = start + [1:length(V)].*samper;
    period = round(0.01662./samper);
    tv = 1:(period+10);
    d = tmains - start;
    id = find(d > 0); 
    t = tmains(id(1));
    tid = find(vt >= t);
    ti = tid(1); % first sampale
    x = ti:ti+tv;
    ts(length(vt)) = 0;
    ts(ti+tv-1) = tv;
    ts(1:ti-1) = tv(end-ti+2:end);
    mv = V(ti+tv-1);
    mt(1) = ti;
    j = 2;
    while t < vt(end) && ti+2*period+100 < length(vt)
        tid = find(vt(ti+period-5:ti+period+5) > tmains(id(j)));
        ti = tid(1)+ti+period-5;
        x = ti:ti+tv;
        ts(ti+tv-1) = tv;
        mv = mv+V(ti+tv-1);
        mt(j) = ti;
        j = j+1;
    end
    ts(ti+tv(end):end) = 1;
    mv = mv./j;
    mv = mv - mean(mv);
    G = sum(mv.*mv);
    GetFigure('FullV');
    hold off;
    plot(V);
    hold on; 
    rmv = zeros(size(V));
    for j = 1:length(mt)-1
        v = V(mt(j):mt(j+1)-1);
        g = sum(v .* mv(1:length(v)))./G;
        rmv(mt(j):mt(j+1)-1) = mv(1:length(v)).*g;
    end
    plot(rmv,'r');
    V = V-rmv;
end
if checkrw
    vt = start + [1:length(V)].*samper;
    id = find(Ev.codes(:,1) == 113 & Ev.times < vt(end) & Ev.times > vt(1));
    for j = 1:length(id)
        ti = find(vt >= Ev.times(id(j)));
        ts(j) = ti(1);
        rV(j,:) = V(ti(1):ti(1)+npts);
    end
    aV = rV;
    GetFigure('FullV');
    plot(rV');
    dv = diff(rV(:,1:end-100),1,2);
    vcrit = std(dv(:)).*3;
    [a,b] = find(dv > vcrit);
    na = 1;
    scrit = 2;
    for j = 1:length(a)
        s(j) = sum(rV(a(j),b(j)+2:b(j)+20));
        sa(j) = dv(a(j),b(j));
        if s(j) > scrit
            af(na,:) = rV(a(j),b(j):b(j)+100);
            na = na+1;
            V(ts(a(j))+b(j)-1:ts(a(j))+b(j)+50) = 0;
        end
    end
    plot(aV');
    return;
end
smv = smooth(V,80);
GetFigure('FullV');
hold off;
plot(V);
hold on;
plot(smv,'r');
V = V-smv;
plot(V,'m');
res.sd = std(V);
if prct ~= 0
    th = prctile(V,prct);
elseif th > 50 %percentile
    th = prctile(V,th);
end

maxl = lens(ispk);

sgn = diff(sign(diff(V)));
if th(1) < 0
id = find(sgn > 0 & V(2:end-1) < th(1));
else
id = find(sgn < 0 & V(2:end-1) > th(1));
end
if isempty(id)
    id = find(sgn > 0)+1;
    prc = 100000./length(id); % get 1000 spikes
    id = id(V(id) < prctile(V(id),prc));
end
GetFigure('Spikes');
hold off;
plot([0 40],[0 0]);
hold on;
allid = [];
id = id(id > 8 & id < maxl-32);
vtimes = start+id.*samper;
allid = repmat(id',41,1) + repmat([-8:32]',1,length(id));
%for j = 1:length(id)
%    allid = [allid [id(j)-8:id(j)+32]'];
%end

ns = 1;
for p = 1:length(C)
f = fields(C{p});
for j = 1:length(f)
    if strncmp(C{p}.(f{j}).title,'Spike',5)
        offset = round((C{p}.(f{j}).start - start)./C{p}.(f{j}).interval); 
        id = find(allid(1,:)-offset > 0 & allid(end,:)-offset <= lens(ns));
        avg(ns,:) = mean(C{p}.(f{j}).values(allid(:,id)-offset),2);
        AllV(ns,:,id) = C{p}.(f{j}).values(allid(:,id)-offset);
        ns = ns+1;
    end
end
end

if sumreject > 0
    sumv = squeeze(mean(AllV,1));
    mv = mean(sumv);
    id = find(mv < std(mv) .* sumreject);
    if plotv == 2
        GetFigure('Spikes');
        bid = setdiff(1:size(AllV,3),id);
        plot(squeeze(AllV(ispk,:,bid)),'g');
    end
    AllV = AllV(:,:,id);
end

end

if fullcov
TV = AllV(1,:,:);
for j = 2:length(chspk)
    TV = cat(2,TV,AllV(j,:,:));
end
TV = squeeze(TV)';
elseif ispk == size(AllV,1)
TV = squeeze(cat(2,AllV(ispk-2,:,:),AllV(ispk-1,:,:),AllV(ispk,:,:)))';
elseif ispk == 1
TV = squeeze(cat(2,AllV(ispk,:,:),AllV(ispk+1,:,:),AllV(ispk+2,:,:)))';
else
TV = squeeze(cat(2,AllV(ispk-1,:,:),AllV(ispk,:,:),AllV(ispk+1,:,:)))';
end

GetFigure('Covar');
imagesc(cov(TV));
[pc, E] = eig(cov(TV));
pcs = TV*pc;
res.Evec = pc;
GetFigure('Means');
subplot(2,1,1);
if length(avg)
    plot(avg');
subplot(2,1,2);
imagesc(avg);
end

if autocut
for j = 1:size(pcs,2)
    dip(j) = HartigansDipTest(sort(pcs(:,j)));
end
res.dipvals = dip;
[dipval, best] = max(dip);

if length(clusterprops)
    p = 1+size(pcs,2) - clusterprops(1:2);
    cx = clusterprops(3);
    cy = clusterprops(4);
    rx = clusterprops(5);
    ry = clusterprops(6);
    r = ((pcs(:,p(1))-cx)./rx).^2 + ((pcs(:,p(2))-cy)./ry).^2;
    id = find(r < 1);
    nid = find(r>1);
else
    npeaks = 0;
w = std(pcs(:,best))./10;
while npeaks ~= 2
    [pdf,x] = smhist(pcs(:,best),'sd',w);
    sgn = diff(sign(diff(pdf)));
    id = find(sgn < 0);
    id = id(find(id > length(pdf) *0.05 & id < length(pdf)*0.95));
    npeaks = length(id);
    if npeaks > 2
        w = w*2;
    elseif npeaks < 2
        w = w*0.8;
    end
    
    
end
[minval, dippos] = min(pdf(id(1):id(2)));
dippos = dippos+id(1);
a = WeightedSum(x(1:dippos),pdf(1:dippos));
b = WeightedSum(x(dippos:end),pdf(dippos:end));
if abs(b) > abs(a)
    id = find(pcs(:,best) > x(dippos));
    nid = find(pcs(:,best) < x(dippos));
else
    id = find(pcs(:,best) < x(dippos));
    nid = find(pcs(:,best) > x(dippos));
end
end
for j = 1:size(AllV,1)
    avg(j,:) = mean(AllV(j,:,id),3);
end
GetFigure('Means');
subplot(2,1,1);
plot(avg');
subplot(2,1,2);
imagesc(avg);
else
    nid = 1:size(AllV,3);
    id = [];
end

if plotv
GetFigure('Spikes');
if plotv == 2
    plot(squeeze(AllV(ispk,:,nid)),'b');
    hold on;
    plot(squeeze(AllV(ispk,:,id)),'r');
else
for j = 1:size(AllV,3)
    if ismember(j,id)
        plot(AllV(ispk,:,j),'r');
    else
        plot(AllV(ispk,:,j));
    end
end
end

res.avg = avg;
for j = 1:size(avg,1)
    scale(j) = sum(avg(j,:) .* avg(ispk,:));
end
res.scale = scale./scale(ispk);
res.spk = spk;


title(sprintf('Probe %d, %d spikes',spk,length(id)));
if rmcluster
[pc, E] = eig(cov(TV(nid,:)));
pcs = TV(nid,:)*pc;
id = [];
end
else
id = [];
end
GetFigure('PCs');
subplot(2,4,1);
PlotPCs(pcs,1,2,clplot,id);

subplot(2,4,2);
PlotPCs(pcs,1,3,clplot,id);

subplot(2,4,3);
PlotPCs(pcs,1,4,clplot,id);

subplot(2,4,4);
PlotPCs(pcs,1,5,clplot,id);
subplot(2,4,5);
PlotPCs(pcs,2,3,clplot,id);
subplot(2,4,6);
PlotPCs(pcs,2,4,clplot,id);
subplot(2,4,7);
PlotPCs(pcs,2,5,clplot,id);
subplot(2,4,8);
PlotPCs(pcs,3,5,clplot,id);
res.pcs = pcs;

if length(clusterprops)
    subplot(2,4,5); %need to find right graph
    hold on;
    DrawEllipse(clusterprops(3:6),'r');
    hold off; 
end
res.energy = squeeze(sum(diff(AllV,[],2).^2,2));
res.spkvar = squeeze(var(AllV,[],2));


function DrawEllipse(E,varargin)

a = E(3); %x radius
b = E(4);
sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)]+E(1);
y = [y fliplr(-y) -y fliplr(y)]+E(2);
plot(real(x),real(y),varargin{:});


function PlotPCs(pcs, a,b, type, id)
ptsz = 1;

if type == 0
    hold off;
    plot(pcs(:,end-a+1),pcs(:,end-b+1),'.','markersize',ptsz);
    if length(id)
        hold on;
        plot(pcs(id,end-a+1),pcs(id,end-b+1),'r.','markersize',ptsz);
    end 
else
    DensityPlot(pcs(:,end-a+1),pcs(:,end-b+1),'sd',[2 2]);
end


function [chspk, lens, starts, offsets, chnames, X] = CountChans(C, spk);

samper=1./40000;
ns = 1;
names = {};
MAXOVERLAP=0.0015;
for p = 1:length(C)
    f = fields(C{p});
    if isfield(C{p},'name');
        [a,b,c] = fileparts(C{p}.name);
        id = strfind(C{p}.name,'A.spkblk');
        if isempty(id)
            id = strfind(C{p}.name,'.spkblk');
        end
        b = C{p}.name(1:id(1)-1);
        if isempty(strmatch(b,names,'exact'))
            names = {names{:} b};
        end
    end
    for j = 1:length(f)
        if isfield(C{p}.(f{j}),'title') & strncmp(C{p}.(f{j}).title,'Spike',5)
            chspk(ns) = sscanf(C{p}.(f{j}).title,'Spike %d');
            lens(ns) = length(C{p}.(f{j}).values);
            starts(ns) = C{p}.(f{j}).start;
            chnames{ns} = f{j};
            sampers(ns) = C{p}.(f{j}).interval;
            ix(ns) = p;
            ns = ns+1;
        end
    end
end

X.names = names;
ends = starts + (lens-1).*samper;
good = ones(size(ends));
for j = unique(chspk)
    chid = find(chspk ==j);
    [s, order] = sort(starts(chid));
    id = find(starts(chid(order(2:end))) < ends(chid(order(1:end-1)))-MAXOVERLAP);
    good(chid(order(id))) = 0;
end
bad = find(good ==0);
good = find(good);
chspk = chspk(good);
lens = lens(good);
starts = starts(good);
chnames = chnames(good);
ix = ix(good);
ends = ends(good);
[a,b] = Counts(sampers);
[a,id] = max(a);
X.samper = sampers(id);

laststart = 0;
nblk = 0;
useid = [];
ul = 0;
vs = 1;
id = find(ends > laststart);
while ~isempty(id)
    nblk = nblk+1;
    blkend(nblk) = min(ends(id));
    vs = vs+ul;
    id = find(starts < blkend(nblk));
    start = max(starts(id)); %defines block
    id = find(ends >= blkend(nblk) & starts <= start);
    blkstart(nblk) = start;
    kid = length(useid)+[1:length(id)];
    offsets(kid) = round((start- starts(id))./samper);
    uselen(kid) = min(lens(id)-offsets(kid));
    X.totaloffsets(kid) = vs;
    ul = max(uselen(kid));
    nextstart = min(starts(starts > min(ends(id))));
    useid(kid) = id;
    X.blklen(nblk) = round(mean(uselen(kid)));
    X.blks(kid) = nblk;
    if isempty(nextstart)
        blkend(nblk) = min(ends(id));
        id = [];
    else
    laststart = max(ends(ends < nextstart));
    blkend(nblk) = laststart;
    id = find(ends > laststart);
    end
end
starts = starts(useid); %some files have > 1 chunk. These are listed n times
chspk = chspk(useid);  
chnames = chnames(useid);
X.ix = ix(useid);
X.blkstart = blkstart;
    [start, a] = min(blkstart);
%    offsets  = round((start- starts)./samper);
id = find(uselen > lens(useid));
lens = uselen;
    

function [chspk, lens, starts, offsets, chnames, X] = CountProbeChans(C, spk);



k = 0;


MAXOVERLAP=0.001;

samper=1./40000;
ns = 1;
names = {};
suffs = unique([C.suffix]);
for p = 1:length(C)
    [a,b,c] = fileparts(C(p).file);
    id = strfind(C(p).file,'A.spkblk');
    id = []; %
    if isempty(id)
        id = strfind(C(p).file,'.spkblk');
    end
    b = [C(p).file(1:id(1)-1) '.spkblk' num2str(C(p).fileid)];
    nid = strmatch(b,names,'exact');
    if isempty(nid)
        names = {names{:} b};
        nid = length(names);
    end
    chspk(ns) = C(p).probe;
    if isfield(C,'len')
        lens(ns) - C(p).len;
    else
    lens(ns) = (C(p).end - C(p).start)./samper;
    end
    starts(ns) = C(p).start;
    chnames{ns} = C(p).var;
    sampers(ns) = samper;
    ix(ns) = nid;
    pid(ns) = p;
    ns = ns+1;
end

X.names = names;
ends = starts + (lens-1).*samper;
good = ones(size(ends));
for j = unique(chspk)
    chid = find(chspk ==j);
    [s, order] = sort(starts(chid));

    id = find(starts(chid(order(2:end))) < ends(chid(order(1:end-1)))-MAXOVERLAP);
    bzid = find(starts(chid(order(2:end))) < ends(chid(order(1:end-1))));
    if ~isempty(id)
        X = AddError(X,'CH%d %d blocks end(n) > start(n+1). %d by > %.0fms  (Starts %s) (Ends %s)\n',...
            chspk(j),length(bzid),length(id), MAXOVERLAP,...
            sprintf('%.4f ',starts(chid(order(id+1)))),sprintf('%.4f ',ends(chid(order(id)))));
        
    elseif ~isempty(bzid)
        fprintf('CH%d %d blocks end(n) > start(n+1). But all < 2ms  (Starts %s) (Ends %s)\n',chspk(j),length(bzid), ...
            sprintf('%.4f ',starts(chid(order(bzid+1)))),sprintf('%.4f ',ends(chid(order(bzid)))));
    end
    for k = 1:length(id)
        a = chid(order(id(k)));
        b = chid(order(id(k)+1));
        if (ends(a)-starts(a)) > (ends(b)-starts(b))
            good(b) = 0;
        else
            good(a) = 0;
        end
    end
end
bad = find(good ==0);
good = find(good);
chspk = chspk(good);
lens = lens(good);
starts = starts(good);
chnames = chnames(good);
ix = ix(good);
ends = ends(good);
[a,b] = Counts(sampers);
[a,id] = max(a);
X.samper = sampers(id);

laststart = 0;
nblk = 0;
useid = [];
ul = 0;
vs = 1;
id = find(ends > laststart);
while ~isempty(id)
    nblk = nblk+1;
    blkend(nblk) = min(ends(id));
    vs = vs+ul;
    id = find(starts < blkend(nblk)-MAXOVERLAP);
    start = max(starts(id)); %defines block
    id = find(ends >= blkend(nblk) & starts <= start);
    blkstart(nblk) = start;
    kid = length(useid)+[1:length(id)];
    offsets(kid) = round((start- starts(id))./samper);
    uselen(kid) = min(lens(id)-offsets(kid));
    X.totaloffsets(kid) = vs;
    ul = max(uselen(kid));
    nextstart = min(starts(starts > min(ends(id))-MAXOVERLAP));
    useid(kid) = id;
    X.blklen(nblk) = round(mean(uselen(kid)));
    X.blks(kid) = nblk;
    nchans(nblk) = length(id);
    if isempty(nextstart)
        blkend(nblk) = min(ends(id));
        id = [];
    else
        laststart = max(ends(ends < nextstart+MAXOVERLAP));
        blkend(nblk) = laststart;
        id = find(ends > laststart+MAXOVERLAP);
    end
end
starts = starts(useid); %some files have > 1 chunk. These are listed n times
chspk = chspk(useid);  
chnames = chnames(useid);
X.ix = ix(useid);
X.pid = pid(useid);
X.nchans = nchans;
X.blkstart = blkstart;
 [start, a] = min(blkstart);
%    offsets  = round((start- starts)./samper);
offsets = round(offsets);
id = find(uselen > lens(useid));
lens = round(uselen);


function PlotV(a,b,id)
DATA = GetDataFromFig(a);

GetFigure('RWvoltages');
plot(DATA.V(:,DATA.rti(id)-100:DATA.rti(id)+200)');

fprintf('Sample %d at %d\n',id,DATA.rti(id));
        