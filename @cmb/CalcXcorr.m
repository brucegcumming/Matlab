function xc = CalcXcorr(DATA, eids, sa, sb)
dt = [-1000:5:1000];
sign = DATA.syncsign;
colors = 'bmr';
if DATA.spikelist(1) >= 0
spikelist = DATA.spikelist;
else
spikelist = [0:4];
end
xc = [];
fprintf('Calculating Cross Correlation %d,%d\n',sa,sb);
tic;
for j = 1:length(eids)
Expt =  DATA.Expts{eids(j)};
trange = [Expt.Trials(1).Start(1) Expt.Trials(end).End(end)];
if isfield(DATA,'AllClusters')
eid = eids(j);
aSpikes =DATA.AllClusters{eid}(sa); 
bSpikes =DATA.AllClusters{eid}(sb); 
else
aSpikes = DATA.AllSpikes{sa};
bSpikes =  DATA.AllSpikes{sb};
end
aid = find(aSpikes.times > trange(1) & aSpikes.times < trange(2));
bid = find(bSpikes.times > trange(1) & bSpikes.times < trange(2));
%        ts = now;
%      xcorrtimes(aSpikes.times(aid),bSpikes.times(bid));
%    mytoc(ts); %not faster
isspk = [0 0];
if sum(aSpikes.codes(aid,2)) == 0 %% no clusters cut
aspikelist = 0;
else
aspikelist = spikelist;
isspk(1) = 1;
end
if sum(bSpikes.codes(bid,2)) == 0 %% no clusters cut
bspikelist = 0;
else
bspikelist = spikelist;
isspk(2) = 1;
end

for t = 1:length(Expt.Trials)
durs(t) = Expt.Trials(t).End(end) - Expt.Trials(t).Start(1);
end
ptimes = 0:10:mean(durs);
for t = length(Expt.Trials):-1:1
trange = [Expt.Trials(t).Start(1)-0.1 Expt.Trials(t).End(end)+0.1];
if sign == 1 %use only positive triggers
aid = find(aSpikes.times > trange(1) & aSpikes.times < trange(2) & ...
aSpikes.values(:,9) > 0 & ismember(aSpikes.codes(:,2),aspikelist));
at = aSpikes.times(aid);
bid = find(bSpikes.times > trange(1) & bSpikes.times < trange(2) & ...
bSpikes.values(:,9) > 0 & ismember(bSpikes.codes(:,2),bspikelist));
bt = bSpikes.times(bid);
elseif sign < 0 %use only negative triggers
aid = find(aSpikes.times > trange(1) & aSpikes.times < trange(2) & ...
aSpikes.values(:,9) < 0 & ismember(aSpikes.codes(:,2),aspikelist));
at = aSpikes.times(aid);
bid = find(bSpikes.times > trange(1) & bSpikes.times < trange(2) & ...
bSpikes.values(:,9) < 0 & ismember(bSpikes.codes(:,2),aspikelist));
bt = bSpikes.times(bid);
else
aid = find(aSpikes.times > trange(1) & aSpikes.times < trange(2) & ismember(aSpikes.codes(:,2),aspikelist));
at = aSpikes.times(aid);
bid = find(bSpikes.times > trange(1) & bSpikes.times < trange(2) & ismember(bSpikes.codes(:,2),aspikelist));
bt = bSpikes.times(bid);
end
dts = [];
if length(at) & length(bt)
for k = length(at):-1:1
dts(k,:) = hist(at(k) - bt,dt);
end
xc(t,:) = sum(dts,1);
psth{1}(t,:) = hist(at-Expt.Trials(t).Start(1),ptimes);
psth{2}(t,:) = hist(bt-Expt.Trials(t).Start(1),ptimes);
end
counts(t,:) = [length(at) length(bt)];
end
end
if length(xc)
xc = sum(xc);
psth{1} = sum(psth{1},1);
psth{2} = sum(psth{2},1);
fprintf('%d,%d spikes, took %.2f (%s)\n',sum(psth{1}),sum(psth{2}),toc,datestr(now));
if sa == sb
xc(find(dt == 0)) = 0;
end
bar(dt(2:end-1)./10,xc(2:end-1),1,colors(sum(isspk)+1));
fgf
end

