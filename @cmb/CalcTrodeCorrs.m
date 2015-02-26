function OUT = CalcTrodeCorrs(DATA,eid)
cell = 0;
mu=0;
clabel = {};
spikes = [];
GetFigure('MeanSpike');
colors = mycolors;
hold off;
GetFigure('MeanMU');
subplot(2,1,1);
hold off;

for p = 1:length(DATA.probelist)
DATA = cmb.SetProbe(DATA,DATA.probelist(p));
eid = eid(find(eid <= length(DATA.Expts)));
for e = 1:length(eid)
[DATA, allspks] = SetExptSpikes(DATA,eid(e),'setrange');
for j = 1:length(allspks);
xc = corrcoef(reshape(DATA.AllData.Spikes.values(allspks(j),:),32,4));
xcs(j,:) = [xc(1,2) xc(2,3) xc(3,4) xc(1,3) xc(2,4) xc(1,4)];
end
tx(p,e,:) = prctile(xcs,50);
ds(p,e,:) = [1 1 1 2 2 3];
cid = unique(DATA.AllData.Spikes.codes(allspks,2));
for j = 2:length(cid)
cell = cell+1;
id = find(DATA.AllData.Spikes.codes(allspks,2) ==cid(j));
mspk = mean(DATA.AllData.Spikes.values(allspks(id),:));
spikes(cell).v = mspk;
spikes(cell).probe = p;
spikes(cell).cluster = cid(j);
GetFigure('MeanSpike');
plot(mspk,'color',colors{cell});
hold on;
mspk = reshape(mspk,32,4);
[a,c] = max(var(mspk));
amp = mspk' * mspk(:,c);
amp = amp./amp(c);
cd = abs([1:4] - c);
cds(cell,:) = cd(cd > 0);
ccs(cell,:) = amp(cd > 0);
clabel{cell} = sprintf('P%dC%d',p,cid(j));
end
ts = [9 41 73 105];
v = prctile(DATA.AllData.Spikes.values(allspks,ts),5);
for j = 1:length(ts)
mu = mu+1;
id = find(DATA.AllData.Spikes.values(allspks,ts(j)) < v(j));
mspk = mean(DATA.AllData.Spikes.values(allspks(id),:));
mus(mu).v = mspk;
mus(mu).probe = p;
mus(mu).subprobe = j;
GetFigure('MeanMU');
subplot(2,1,1);
plot(mspk,'color',colors{j});
hold on;
mspk = reshape(mspk,32,4);
amp = mspk' * mspk(:,j);
amp = amp./amp(j);
cd = abs([1:4] - j);
mus(mu).amp = amp(cd > 0);
mus(mu).d = cd(cd > 0);
muds(mu,:) = cd(cd > 0);
mucs(mu,:) = amp(cd > 0);
end
end
end
GetFigure('MeanSpike');
legend(clabel);
GetFigure('MeanMU');
subplot(2,1,2);
hold off;
for j = 1:4
id = find([mus.subprobe] == j);
plot(mean(cat(1,mus(id).v)),'linewidth',2,'color',colors{j});
hold on;
end
GetFigure('TrodeCorr');
hold off;
plot(ds(:),tx(:),'o');
hold on;
if exist('cds','var')
plot(cds(:),ccs(:),'ro');
else 
cds = [];
ccs = [];
end
dval = unique(ds);
plot(muds(:),mucs(:),'go');
for j = 1:length(dval)
id = find(cds == dval(j));
cmean(j) = mean(ccs(id));
id = find(ds == dval(j));
dmean(j) = mean(tx(id));
id = find(muds == dval(j));
mumean(j) = mean(mucs(id));
end
plot(dval,dmean,'-');
plot(dval,cmean,'r-');
plot(dval,mumean,'g-');

DATA.TrodeMean.spikes = spikes;
DATA.TrodeMean.mu = mus;
OUT = DATA;

