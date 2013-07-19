function [pavg, spk] = FixSpkMains(spk, tics, varargin)

%pavg = FixSpkMains(spk, tics, ...)
%calculate mean resp triggered on Mains signal
%
%[pavg, spks]  = FixSpkMains(spk, tics, 'fix')
%removes the avg from the traces

period = median(diff(tics));
if period > 100  %in tics, not secs
npts = ceil(period./(spk.resolution.*10000));
else
npts = ceil(period./spk.resolution);
end
n = 0;
k = 1;
dofix = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fix',3)
        dofix = 1;
    end
    j = j+1;
end
aid = [];
bid = [];
bts = [1:101];
pn = zeros(1,npts);
psum = zeros(size(pn));
pavg = psum;
nb = length(spk.times);
id = find(tics > spk.times(1));
start =  id(1);
if start > 1
    start = start-1;
end
if start == 0
    fprintf('Tics all before Spikes\n');
    return;
end
ns = size(spk.values,2);
allpts = zeros(size(spk.values));
for j = start:length(tics)
    t = spk.times(bts)-(tics(j)+period);
    while(t(end) < 0 & bts(end) < nb)
        if bts(end)+100 > nb
            bts = bts(end):nb;
        else
            bts = bts+100;
        end
        t = spk.times(bts)-(tics(j)+period);
    end
    b = find(t > 0 & t < period);
    for k = 1:length(b)
        pt = 1+round(npts.* t(b(k))./period);
        pts = pt:pt+ns-1;
        pts = 1+mod(pts-1,npts);
        allpts(bts(b(k)),:) = pts;
        pn(pts) = pn(pts)+1;
        psum(pts) = psum(pts)+spk.values(bts(b(k)),:);
    end
end
pavg = psum./pn;
if dofix
    id = find(allpts(:,1) > 0);
    spk.values(id,:) = spk.values(id,:) - pavg(allpts(id,:));
end

