function xc = CalcXcorrV(DATA, eids, sa, sb)
dt = 2;
pts = 24:32;
colors = 'bmr';
fprintf('Calculating V Cross Correlation %d,%d\n',sa,sb);
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
dcs(1,:,:) = DATA.AllSpikes{sa}.values(aid(ai),:);
dcs(2,:,:) = DATA.AllSpikes{sb}.values(bid(bi),:);
end
if ~isempty(pts)
a = mean(dcs(1,:,pts),3);
b = mean(dcs(2,:,pts),3);
scatter(a(:),b(:),'.');
else
scatter(dcs(1,:),dcs(2,:),'.');
end



