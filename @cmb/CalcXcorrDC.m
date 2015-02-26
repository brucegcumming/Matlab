function xc = CalcXcorrDC(DATA, eids, sa, sb)
dt = 2;
colors = 'bmr';
fprintf('Calculating DC Cross Correlation %d,%d\n',sa,sb);
tic;
for j = 1:length(eids)
Expt =  DATA.Expts{eids(j)};
trange = [Expt.Trials(1).Start(1) Expt.Trials(end).End(end)];
if DATA.syncsign > 0
aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
DATA.AllSpikes{sa}.values(:,9) > 0);
bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
DATA.AllSpikes{sb}.values(:,9) > 0);
elseif DATA.syncsign < 0
aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
DATA.AllSpikes{sa}.values(:,9) < 0);
bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
DATA.AllSpikes{sb}.values(:,9) < 0);
else 
aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2));
bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2));
end

[ai,bi] = cmb.FindSync(DATA.AllSpikes{sa}.times(aid),...
DATA.AllSpikes{sb}.times(bid),dt);
dcs(1,:) = mean(DATA.AllSpikes{sa}.values(aid(ai),24:end),2);
dcs(2,:) = mean(DATA.AllSpikes{sb}.values(bid(bi),24:end),2);
end
scatter(dcs(1,:),dcs(2,:),'.');


