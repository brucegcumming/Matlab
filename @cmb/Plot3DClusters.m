function DATA = Plot3DClusters(a,b);
recalc = 0;
if isfield(a,'AllData');
DATA = a;
recalc = b;
else
DATA = GetDataFromFig(a);
end


[xyzfig, isnew] = GetFigure('Cluster3D');
if isnew
if ~isfield(DATA.plot,'clusterZ')
DATA.plot.clusterZ = 5;
end
z.parentfigtag = DATA.tag.top;
z.parentfig = DATA.toplevel;
set(xyzfig,'UserData',z);
bp = [10 10 100 DATA.plot.ch];
xyfig = xyzfig;
hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
for j = 1:length(DATA.spkvarorder)
k = DATA.spkvarorder(j);
uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 1, 3, k});
end
hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
for j = 1:length(DATA.spkvarorder)
k = DATA.spkvarorder(j);
uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 2,3, k});
end
hm = uimenu(xyfig,'Label','Z','Tag','3DplotMenu');
uimenu(hm,'Label','None','Callback',{@cmb.SetXYCluster, 3, 0});
for j = 1:length(DATA.spkvarorder)
k = DATA.spkvarorder(j);
uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@cmb.SetXYCluster, 3, 3, k});
end


DATA.figs.xyz = xyzfig;
end
Spks = DATA.AllData.Spikes;
ispk = DATA.spklist;
DATA = CheckForPCA(DATA, ispk, 0);
if isempty(DATA.AllData.pcs)
PCs = [];
else
PCs = DATA.AllData.pcs(ispk,:);
end
classify = 1;
[cz, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterZ, classify,PCs);
DATA.Spikes.cz(ispk) = cz;
if recalc
[cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs);
DATA.Spikes.cx(ispk) = cx;
[cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs);
DATA.Spikes.cy(ispk) = cy;
end
cl = unique(Spks.codes(ispk,2));
hold off;
for j = 1:length(cl)
id = find(Spks.codes(ispk,2) == cl(j));
plot3(DATA.Spikes.cx(ispk(id)),DATA.Spikes.cy(ispk(id)),cz(id),'.',...
'color',DATA.spkcolor{cl(j)+1},'markersize',DATA.ptsize);
hold on;
end
xlabel(['X:' DATA.spkvarnames{DATA.plot.clusterX}]);
ylabel(['Y:' DATA.spkvarnames{DATA.plot.clusterY}]);
zlabel(['Z:' DATA.spkvarnames{DATA.plot.clusterZ}]);


