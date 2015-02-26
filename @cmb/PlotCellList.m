function PlotCellList(DATA, varargin)

QUALITY = 1;
CELLNUM=2;
ONECELL=3;
JOINCELL=4;
CLUSTERQUALITY=5;
CLUSTERANDCELL=6;
DPRIME=7;
CELLANDDPRIME = 9;
CQDPRIME=8;
replotshape = 0;

cell = 1;
if isfield(DATA,'CellList') && (ndims(DATA.CellList) == 3 || size(DATA.CellList,1) == length(DATA.Expts))
cmb.PlotNewCellList(DATA, varargin)
return;
end
if ~isfield(DATA.plot,'showspikeshapes')
showspikeshape = 0;
else
showspikeshape = DATA.plot.showspikeshapes;
end
if isfield(DATA.plot,'cellplot')
plottype = DATA.plot.cellplot;
else
plottype = QUALITY;
end

n = length(DATA.probelist);
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'bynumber',4)
plottype = CELLNUM;
elseif strncmpi(varargin{j},'onecell',4)
plottype = ONECELL;
elseif strncmpi(varargin{j},'joincell',4)
plottype = JOINCELL;
elseif strncmpi(varargin{j},'clusterquality',7)
plottype = CLUSTERQUALITY;
elseif strncmpi(varargin{j},'plotshapes',8)
replotshape = 1;
end
j = j+1;
end

if ~isfield(DATA,'CellList') | isempty(DATA.CellList)
DATA.CellList = zeros(1,length(DATA.AllData.Trialids));
DATA.CellQuality = zeros(1,length(DATA.AllData.Trialids));
DATA.CellListCluster = zeros(1,length(DATA.AllData.Trialids));
%    return;
end
it = findobj('Tag','CombinerCellList');
it = findobj(it, 'Tag','CellNumber');
if ~isempty(it)
cell = get(it(1),'value');
end
cfig = gcf;
if showspikeshape
fign = findobj('Tag','SpikeShapes','Type','Figure');
if (isempty(fign) | replotshape) & exist(DATA.meanspkfile,'file');
PlotSpikeShapes(DATA.meanspkfile,'cells',DATA.CellList,DATA.CellListCluster);
end
if isfield(DATA, 'exabsid')
PlotSpikeShapes([DATA.exabsid(1) DATA.exabsid(end)]);
end
figure(cfig);
end

cmb.ClearAnnotations(gca);
Q = zeros(length(DATA.probelist),DATA.Expts{end}.Trials(end).Trial);
im = zeros(length(DATA.probelist)*2,DATA.Expts{end}.Trials(end).Trial);
dp = zeros(size(Q));
for j = 1:size(DATA.CellList,1)
id = find(DATA.CellList(j,:) > 0);
k = DATA.CellList(j,id);
ind = sub2ind(size(im),k*2,id);
aid = find(DATA.CellList(j,:) > 0 & DATA.CellListCluster(j,:) <= 1);
if length(aid)
k = DATA.CellList(j,aid)*2;
indb = sub2ind(size(im),k-1,aid);
else
indb = [];
end
bid = find(DATA.CellList(j,:) > 0 & DATA.CellListCluster(j,:) == 2 & DATA.CellQuality(j,:) > 0);
if length(bid)
k = DATA.CellList(j,bid)*2;
indc = sub2ind(size(im),k,bid);
else
indc = [];
end
indb = [indb indc];

if isfield(DATA,'CellQuality')
im(ind) = DATA.CellQuality(j,id);
im(indb) = DATA.CellQuality(j,[aid bid]);
else
im(ind) = j;
im(indb) = j;
end
end
np = length(DATA.probelist);
xl = get(gca,'xlim');
yl = get(gca,'ylim');
ap = get(gca,'position');

for e = 1:length(DATA.Expts)
id = [DATA.Expts{e}.Trials.Trial];
if isfield(DATA.Expts{e}.Trials,'ed')
eds(id) = [DATA.Expts{e}.Trials.ed];
else
eds(id) = DATA.Expts{e}.Stimvals.ed;
end
if isfield(DATA.Expts{e},'Cluster') && size(DATA.Expts{e}.Cluster,1) > 0
for k = 1:min([size(DATA.Expts{e}.Cluster,2) np])
if isfield(DATA.Expts{e}.Cluster{1,k},'dprime')
dp(k,id) = DATA.Expts{e}.Cluster{1,k}.dprime;
end
if isfield(DATA.Expts{e}.Cluster{1,k},'quality')
Q(k,id) = DATA.Expts{e}.Cluster{1,k}.quality;
end
end
end
t = [DATA.Expts{e}.Trials(1).Trial DATA.Expts{e}.Trials(end).Trial];
if t(1) >= xl(1) && t(1) <= xl(2)
xp = ap(1)+ap(3) * (t(1)-xl(1))/range(xl);
hlines(j) = annotation('line',[xp xp],[ap(2)+ap(4) 1.0]);
set(hlines(j),'color','r');
end
end

if plottype == QUALITY || plottype == CLUSTERANDCELL
if plottype == CELLANDDPRIME
for j = 1:np
im((j*2)-1,:) = dp(j,:);
end
end
if plottype == CLUSTERANDCELL
for j = 1:np
im((j*2)-1,:) = Q(j,:);
end
end
hold off;
imagesc([1:size(im,2)],[0.75 np+0.25],im);
end
if plottype ==CLUSTERQUALITY
if plottype == CLUSTERANDCELL
subplot(2,1,2);
else
subplot(1,1,1);
end
Q = zeros(length(DATA.probelist),DATA.Expts{end}.Trials(end).Trial);
for e = 1:length(DATA.Expts)
id = [DATA.Expts{e}.Trials.Trial];
if isfield(DATA.Expts{e}.Trials,'ed')
eds(id) = [DATA.Expts{e}.Trials.ed];
else
eds(id) = DATA.Expts{e}.Stimvals.ed;
end
if isfield(DATA.Expts{e},'Cluster') && size(DATA.Expts{e}.Cluster,1) > 0
for k = 1:size(DATA.Expts{e}.Cluster,2)
if isfield(DATA.Expts{e}.Cluster{1,k},'quality')
Q(k,id) = DATA.Expts{e}.Cluster{1,k}.quality;
end
end
end
end
hold off;
imagesc(Q);
caxis([0 8]);
hold on;
cl = (DATA.CellListCluster-1.5)/2;

plot((DATA.CellList')+cl','linewidth',2);
legend(num2str([1:size(DATA.CellList,1)]'),'location','NorthWest');
set(gca,'YaxisLocation','right');
[a,b] = cmb.TrialRange(DATA);
plot([a a],[1 length(DATA.probelist)],'w:');
plot([b b],[1 length(DATA.probelist)],'w:');
elseif plottype == CELLANDDPRIME
for j = 1:np
im((j*2)-1,:) = dp(j,:);
end
hold off;
imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
caxis([0 7]);
elseif plottype == CQDPRIME
for j = 1:np
im((j*2)-1,:) = dp(j,:);
im((j*2),:) = Q(j,:);
end
hold off;
imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
caxis([0 7]);
elseif plottype ==DPRIME
for j = 1:np
im((j*2)-1,:) = dp(j,:);
im((j*2),:) = dp(j,:);
end
hold off;
imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
caxis([0 7]);
elseif plottype == CELLNUM
DATA.CellList(find(DATA.CellList == 0)) = NaN;
hold off;
plot(DATA.CellList','linewidth',2);
set(gca,'ylim',[min(DATA.CellList(:))-1 max(DATA.CellList(:))+1]);
legend(num2str([1:size(DATA.CellList,1)]'),'location','SouthWest')
elseif plottype == ONECELL
trials = find(DATA.CellList(cell,:) > 0);
hold off;
plot(trials,DATA.CellList(cell,trials),'o-');
elseif plottype == JOINCELL
hold off;
colors = mycolors;
mid = (size(DATA.CellList,1)-1)/2;
for j = 1:size(DATA.CellList,1)
trials = find(DATA.CellList(j,:) > 0);
plot(trials,DATA.CellList(j,trials)+0.05*(j-mid),'o-','color',colors{j});
hold on;
end
[a,b] = cmb.TrialRange(DATA);
plot([a a],[1 length(DATA.probelist)],'k:');
plot([b b],[1 length(DATA.probelist)],'k:');

legend(num2str([1:size(DATA.CellList,1)]'),'location','SouthWest')
hold off;
set(gca,'ylim',[min(DATA.CellList(:))-1 max(DATA.CellList(:))+1]);
end

if ~ismember(plottype,[CELLNUM])
hold on;
DATA.CellList(find(DATA.CellList == 0)) = NaN;
cl = (DATA.CellListCluster-0.5)/2;

h = plot((DATA.CellList')+cl','linewidth',2);
for j = 1:length(h)
colors{j} = get(h(j),'color');
end
legend(num2str([1:size(DATA.CellList,1)]'),'location',[0 0.5 0.08 0.5])
set(gca,'YaxisLocation','right');
[a,b] = cmb.TrialRange(DATA);
plot([a a],[1 length(DATA.probelist)],'w:');
plot([b b],[1 length(DATA.probelist)],'w:');
id = find(eds > 0);
%plot depth with sign inverte, so that it matches predicted movement
%of cells on probe. Ingreasing depth = cells move to lowe probe #s
if range(eds(id)) > 8 .* 0.15;
aid = find(eds(a:b) >0) + a -1;
a = mean(eds(aid));
if a-min(eds(id)) < 0.65
a = min(eds(id)+0.65);
end
%add 0.6 because y axis starts at 0.5;
plot(id, 0.6+((a-eds(id))./0.15),'w');
else
plot(id, (max(eds(id))-eds(id))./0.15,'w');
end
for j = 1:size(DATA.CellList,1)
id = find(DATA.CellList(j,:) > 0);
if length(id)
text(id(1),DATA.CellList(j,id(1))+0.1,num2str(j),'color',colors{j});
did = find(abs(diff(DATA.CellList(j,id)))>0)+1;
for k = 1:length(did)
text(id(did(1)),DATA.CellList(j,id(did(1)))+0.1,num2str(j),'color',colors{j});
end
end
end
end

