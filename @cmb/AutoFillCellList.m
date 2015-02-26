function AutoFillCellList(caller,b, varargin)

DATA = GetDataFromFig(caller);
eid = DATA.exabsid;
trange(1) = DATA.Expts{eid(1)}.Trials(1).Trial;
trange(2) = DATA.Expts{eid(end)}.Trials(end).Trial;
ncells = size(DATA.CellList,1);
for p = 1:length(DATA.probelist)
ts = [];
for j = 1:length(eid)
E = DATA.Expts{eid(j)};
if cmb.iscluster(E.Cluster,1,p) & isfield(E.Cluster{1,p},'quality')
Q(p,j) = E.Cluster{1,p}.quality;
else
Q(p,j) = NaN;
end
ts(j,:) = [E.Trials(1).Trial E.Trials(end).Trial];
end
for c = 1:size(DATA.CellList,1)
id = find(DATA.CellList(c,:) == p);
nc(c) = length(id);
if nc(c)
cc(c) = c;
else
cc(c) = 0;
end
end
[a,b] = max(nc);
id = find(~isnan(Q(p,:)));
qu = mean(Q(p,id));
did = find(diff(Q(p,id)));
if qu > 0
did(length(did)+1) = length(id);
first  = 1;
if a > 0
c = cc(b);
else qu > 0;
ncells = ncells+1;
c = ncells;
end
cellid(p) = c;
for k = 1:length(did)
qu = mean(Q(p,id(first):id(did(k))));
trange = ts(id(first),1):ts(id(did(k)),2);
DATA.CellList(c,trange) = p;
DATA.CellQuality(c,trange) = qu;
DATA.CellListCluster(c,trange) = 1;
first = did(k)+1;
end
else
cellid(p) = 0;
end
end
set(DATA.toplevel,'UserData',DATA);
cmb.Recmb.PlotCellList(caller,b);

