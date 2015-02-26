function PlotNewCellList(DATA, varargin)
plotmahal = 0;
force = 0;
reload = 0;
offset = 0;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'showfig',5) %Force figure creation
force = 1;
elseif strncmpi(varargin{j},'reload',5) %reload from disk
reload = 1;
end
j = j+1;
end


if force
f = cmb.SetFigure(DATA, DATA.tag.celllist);
else
f = findobj('Tag',DATA.tag.celllist,'type','figure');
if isempty(f)
return;
else
figure(f);
end
end
if reload
%        DATA = cmb.LoadCellFile(DATA);
end

if isempty(DATA.CellList)
return;
end
colors = mycolors;
subplot(1,1,1);
set(f,'UserData',DATA.toplevel);
hold off;
for j = 1:size(DATA.CellList,1)
for k = 1:size(DATA.CellList,2)
id = find(DATA.CellList(j,k,:) > 0);
nclusters(j,k) = length(id);
end
end

nc = max(DATA.nclusters); %max # clusters in any one expt, for each probe

for j = 1:size(DATA.CellList,2)*2
CellIm(:,j) = DATA.CellList(:,ceil(j/2),1+offset);
end
id = find(nc >1);
for j = 1:length(id)
tid = find(DATA.nclusters(:,id(j)) > 1);
CellIm(tid,id(j)*2) = DATA.CellList(tid,id(j),offset+2);
end

imagesc([1 size(DATA.CellList,2)],[1 size(DATA.CellList,1)],CellIm,'buttondownfcn',{@cmb.HitImage, 3});
hold on;
cells = unique(DATA.CellList(:,:,1+offset));
cells = cells(cells > 0);

for j = 1:length(cells)
[x,y] = find(DATA.CellList(:,:,1+offset) == cells(j));
for k = 1:length(y)
if nclusters(x(k),y(k)) > 1
plot([y(k)-0.25 y(k)-0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@cmb.HitImage, 3});
else
plot([y(k) y(k)], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@cmb.HitImage, 3});
end
end
[a,b] = min(x);
h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
set(h,'color',colors{j});
end
if size(DATA.CellList,3) > 1+offset
cells = unique(DATA.CellList(:,:,2+offset));
cells = cells(cells > 0);
for j = 1:length(cells)
[x,y] = find(DATA.CellList(:,:,2+offset) == cells(j));
for k = 1:length(y)
plot([y(k)+0.25 y(k)+0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@cmb.HitImage, 3});
end
[a,b] = min(x);
h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
set(h,'color',colors{j});
end
end
iscellim = sum(DATA.CellList,3) > 0;
p = GetProbe(DATA, DATA.currentexpt(1), DATA.probe)
DATA.currentpoint = [DATA.currentexpt(1) p];
DATA = cmb.MarkCurrentCell(DATA);
set(DATA.toplevel,'UserData',DATA);
%set(gcf, 'KeyPressFcn',{@cmb.KeyPressed,3});

