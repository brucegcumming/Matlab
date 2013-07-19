function PlotCellIm(CellList, CellDetails, nclusters, varargin)
%PlotCellIm(CellList, CellDetails, nclusters, varargin)
%plot cell list image for like PlotClusters
%CellList is a ExPxC matrix
%CellDetails is a structre with exluded trials etc
%nclusters is and ExP matrix recording the max cluster # on that probe;
%          if empty, it will be calculated here


plotmahal = 0;
colors = mycolors;
callback = @HitCellIm;
markmat = {};
nmark = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'callback',6)
        j=j+1;
        callback = varargin{j};
    elseif strncmpi(varargin{j},'colors',6)
        j=j+1;
        colors = varargin{j};
    elseif strncmpi(varargin{j},'mark',4)
        j=j+1;
        nmark = nmark+1;
        markmat{nmark} = varargin{j};
        j = j+1;
        markcolor = varargin{j};
    end
    j=j+1;
end


if strcmp(CellList,'mark')
    [x,y] = find(CellDetails);
    for j = 1:length(x)
        DrawBox(x(j),y(j),callback,'color',nclusters);
    end
    return;
elseif strcmp(CellList,'marklist')
    for j = 1:size(CellDetails,2)
        DrawBox(CellDetails(1,j),CellDetails(2,j),callback,'color',nclusters);
    end
    return;
end
if plotmahal
im = zeros(size(CellList));
id = find(sum(CellList,3) > 0);
X = squeeze(Mahal(:,:,1));
im(id) = X(id);
end
if isempty('nclusters')
    for j = 1:size(CellList,1)
    for k = 1:size(CellList,2)
        id = find(CellList(j,k,:) > 0);
        nclusters(j,k) = length(id);
    end
    end
end
offset = 0;
nc = max(nclusters); %max # clusters in any one expt, for each probe
if plotmahal
    imagesc(im,'buttondownfcn',{@HitImage, 1});
    caxis([0 5]);
else
    for j = 1:size(CellList,2)*2
        CellIm(:,j) = CellList(:,ceil(j/2),1+offset);
    end
    id = find(nc >1);
    for j = 1:length(id)
        tid = find(nclusters(:,id(j)) > 1);
        CellIm(tid,id(j)*2) = CellList(tid,id(j),offset+2);
    end
        
%show locations of duplicate cells with a -1
CellIm(CellIm <0) = -1;
imagesc([0.75 size(CellList,2)+0.25],[0.75 size(CellList,1)+0.25],CellIm,'buttondownfcn',callback);
end
hold on;
cells = unique(CellList(:,:,1+offset));
cells = cells(cells > 0);

for j = 1:length(cells)
    [x,y] = find(CellList(:,:,1+offset) == cells(j));
    for k = 1:length(y)
        if nclusters(x(k),y(k)) > 1
            plot([y(k)-0.25 y(k)-0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',callback);
        else
            plot([y(k) y(k)], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',callback);
        end
    end
    [a,b] = min(x);
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
    set(h,'color',colors{j});
end
cells = unique(CellList(:,:,2+offset));
cells = cells(cells > 0);
for j = 1:length(cells)
    [x,y] = find(CellList(:,:,2+offset) == cells(j));
    for k = 1:length(y)
            plot([y(k)+0.25 y(k)+0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
    end
    [a,b] = min(x);
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
    set(h,'color',colors{j});
end
iscellim = sum(CellList,3) > 0;


function h = DrawBox(ex, p, callback, varargin)

    lx = -0.5;
    j = 1;
    passon = {};
    while j <= length(varargin)
        if strcmp(varargin{j},'box2')
            lx = 0;
        else
            passon = {passon{:} varargin{j}};
        end
        j = j+1;
    end
    
    if length(ex) > 1 || length(p) > 1
        square = [min(p)+lx max(p)+0.5 max(p)+0.5 min(p)+lx min(p)+lx; min(ex)-0.5 min(ex)-0.5 max(ex)+0.5 max(ex)+0.5 min(ex)-0.5];
    else
        square = [lx 0.5 0.5 lx lx; -0.5 -0.5 0.5 0.5 -0.5];
        square = [square(1,:)+p; square(2,:)+ex];
    end
    h = plot(square(1,:),square(2,:),'k-','buttondownfcn',callback,'linewidth',2);
    if length(passon)
        set(h,passon{:});
    end

function HitCellIm(a,b)

ax = get(a,'Parent');
xy = get(ax,'currentpoint');
ex = round(xy(1,2));
p = round(xy(1,1));
if xy(1,1) > p
    setcl = 2;
else
    setcl = 1;
end
