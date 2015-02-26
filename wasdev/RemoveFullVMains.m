function [FullV, MainsTrig] = RemoveFullVMains(FullV, mt, MainsTrig, varargin)


calconly = 0;
newcycle = 0;
mcycle = [];
samples = [];
block = 0;
byblock = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'calconly',8)
        calconly = 1;
    elseif strncmpi(varargin{j},'byblock',6)
        byblock = 1;
    elseif strncmpi(varargin{j},'blocks',6)
        for j = 1:length(FullV.blkstart);
            [a,b] = RemoveFullVMains(FullV, mt, [], 'block',j);
            avgim(j,:) = a.mainsavg;
        end
        FullV = avgim;
        imagesc(avgim);
        return;
    elseif strncmpi(varargin{j},'block',5)
        j = j+1;
        block = varargin{j};
    elseif strncmpi(varargin{j},'segments',8)
        j = j+1;
        seglen = varargin{j};
        starts = 1:seglen:length(FullV.V);
        if length(starts) > 30
            nseg = 30;
        else
            nseg = length(starts)-1;
        end
        for j= 1:nseg
            [a,b] = RemoveFullVMains(FullV, mt, [], 'samples', starts(j):starts(j)+seglen);
            avgim(j,:) = a.mainsavg;
        end
        imagesc(avgim);
        return;
    elseif strncmpi(varargin{j},'samples',8)
        j = j+1;
        samples = varargin{j};
    end
    j = j+1;
end

firstblk = 1;
lastblk = length(FullV.blklen);
if block > 0
    lastblk = block;
    firstblk = block;
end

if ~isempty(samples)
    FullV.V = FullV.V(samples);
end
    
mp = prctile(diff(mt), 90)./10000;  %period of mains cycle in sec
ns = ceil(mp./FullV.samper); %# of samples per mains cycle.
mt = mt./10000;  %work in sec
if ~isfield(FullV,'t')
    first = 1;
    vt(sum(length(FullV.blklen))) = 0;
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        vt(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        first = last+1;
    end
end

tst = now;
if isempty(MainsTrig)
    starts = [1 cumsum(FullV.blklen)];
    ends = starts(2:end)-1;
    mcycle(1:length(FullV.V)) = NaN;
    for j = firstblk:lastblk
        ts = vt(starts(j));
        te = vt(ends(j));
        ti = vt(starts(j):ends(j));
        mid = find(mt > ts & mt < te+0.02);
        k = 1;
        tic;
        nextt = mt(mid(k));
        if mid(k) == 1
            lastt = mt(mid(k))-0.0167;
        else
            lastt = mt(mid(k)-1);
        end
        for t = 1:length(ti)
            if ti(t) > nextt
                k = k+1;
                lastt = nextt;
                nextt = mt(mid(k));
            end
            mcycle(starts(j)+t-1) = ti(t)-lastt;
        end
    end
    FullV.mainsbuild(1) = mytoc(tst);
    mcycle = ceil(mcycle.*30000);
    mcycle = mcycle(1:length(FullV.V));
    MainsTrig.mcycle = mcycle;
    newcycle = 1;
end

if calconly
    for j = ns:-1:1
        id = find(mcycle == j);
        mtrig(j) = mean(FullV.V(id));
        mainsfullv(id) = mtrig(j);
    end
    FullV.fullv = mainsfullv;
    FullV.mtrig = mtrig;
    return;
end

adjgain = 1;
if adjgain
mmean = zeros(size(FullV.V));
rawv = FullV.V;
end
if byblock
    starts = [1 cumsum(FullV.blklen)];
    for b = 1:length(FullV.blklen)
        tid = starts(b)+[1:FullV.blklen(b)]-1;
        for j = ns:-1:1
            if newcycle
                MainsTrig.mids{b,j} = find(mcycle(tid) == j);
            end
            id = MainsTrig.mids{b,j};
            mtrig(j) = mean(FullV.V(tid(id)));
            mmean(tid(id)) = mtrig(j);
            FullV.V(tid(id)) = FullV.V(tid(id))-mtrig(j);
        end
    end
    MainsTrig.mmean = mmean;
else
for j = ns:-1:1
    if newcycle
        MainsTrig.mids{j} = find(mcycle == j);
    end
    id = MainsTrig.mids{j};
    mtrig(j) = mean(FullV.V(id));
    mmean(id) = mtrig(j);
    FullV.V(id) = FullV.V(id)-mtrig(j);
end
end
FullV.mainsbuild(2) = mytoc(tst);
FullV.mainsavg = mtrig; 

