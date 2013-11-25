function [X,D] = ListClusterBackup(prefix,varargin)
%ListClusterBackup(prefix)
%Finds cluster backup files and compares them with current definitions
j = 1;
exptlist = [];
plotargs = {};
D = [];

while j <= length(varargin)
    if strncmpi(varargin{j},'expt',4)
        j = j+1;
        exptlist = varargin{j};
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end
if isfield(prefix,'savetime') %plot previous result
    X = PlotBackup(prefix,plotargs{:});
    return;
end

if isempty(exptlist)
    d = mydir([prefix '/*ClusterTimes*.mat']);
else
    d = mydir(sprintf('%s/Expt%dClusterTimes*.mat',prefix,exptlist(1)));
end
path = fileparts(prefix);
[dates, did] = sort([d.datenum],'descend');
nb = 0;
for j = 1:length(d)
    clear Clusters;
    load(d(did(j)).name);
    if exist('Clusters','var')
        nb = nb+1;
        saves = CellToMat(Clusters,'savetime');
        exs = CellToMat(Clusters,'exptno');
        X(nb).nspks = CellToMat(Clusters,'nspks');
        X(nb).space = CellToMat(Clusters,'space');
        X(nb).ncut = CellToMat(Clusters,'ncut');
        X(nb).marked = CellToMat(Clusters,'marked');
        X(nb).mahal = CellToMat(Clusters,'mahal');
        X(nb).dropi = CellToMat(Clusters,'dropi');
        X(nb).fitdprime = CellToMat(Clusters,'fitdprime');
        X(nb).recluster = CellToMat(Clusters,'recluster');
        X(nb).manual = CellToMat(Clusters,'manual');
        [a,b] = max(saves(:,end));
        X(nb).savetime =a;
        X(nb).savetimes = saves(:,end);
        X(nb).p = Clusters{b}.probe(1);
        if isfield(Clusters{b},'user')
            X(nb).user = Clusters{b}.user;
        else
            X(nb).user = '?';
        end
        X(nb).savecl = b;
        X(nb).expt = median(exs);
        X(nb).name = d(did(j)).name;
    end
%    fprintf('%s %d %d %s\n',d(did(j)).filename,b,Clusters{b}.manual,datestr(a));
end

[a,b] = sort([X.savetime]);
X = X(b);
PlotBackup(X,plotargs{:});

function D = PlotBackup(X, varargin)
j = 1;
plottype = 'saves';

if ishandle(X)
    F = GetFigure(X);
    X = getappdata(F,'ClusterData');
end
while j <= length(varargin)
    if sum(strncmpi(varargin{j},{'recluster' 'savetimes' 'fitdprime'},6))
        plottype = varargin{j};
    elseif sum(strncmpi(varargin{j},{'marked' 'nspks' 'saves' 'ncut' 'space' 'mahal' 'dropi'},4))
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'maxage',6)
        j = j+1;
        id = find(now - [X.savetime] < varargin{j});
        X = X(id);
    end
    j = j+1;
end

[F, isnew] = GetFigure('ClusterBackups');
if isnew
   hm = uimenu(F,'label','Plots');
   uimenu(hm,'label','mahal','callback',{@PlotBackup, 'mahal'});
   uimenu(hm,'label','Fit dprime','callback',{@PlotBackup, 'fitdprime'});
   uimenu(hm,'label','Drop Index','callback',{@PlotBackup, 'dropi'});
   uimenu(hm,'label','Save Times','callback',{@PlotBackup, 'savetimes'});
   uimenu(hm,'label','Recluster type','callback',{@PlotBackup, 'recluster'});
   uimenu(hm,'label','N events','callback',{@PlotBackup, 'nspks'});
   uimenu(hm,'label','N Spikes Cut','callback',{@PlotBackup, 'ncut'});
   uimenu(hm,'label','Marks','callback',{@PlotBackup, 'marked'});
   uimenu(hm,'label','date/probe','callback',{@PlotBackup, 'saves'});
end
   setappdata(F,'ClusterData',X);
D.savetime = [X.savetime];
imh = 0;
hold off;
if sum(strcmp(plottype,{'marked' 'savetimes' 'nspks' 'ncut' 'recluster'}))
    f = plottype;
    for j = 1:length(X)
        marked(j,1:length(X(j).(f))) = X(j).(f);
    end
    imh = imagesc(marked);
    D.(f) = marked;
elseif sum(strcmp(plottype,{'space' 'fitdprime' 'mahal' 'dropi'}))
    f = plottype;
    ix = 1;
    if strcmp(plottype,'mahal')
        ix = 4;
    elseif strcmp(plottype,'dropi')
        ix = 3;
    end
    for j = 1:length(X)
        marked(j,1:length(X(j).(f))) = X(j).(f)(:,ix);
    end
    if strcmp(plottype,'fitdprime')
        marked = -marked;
    end
    imh = imagesc(marked);
    D.(f) = marked;    
else
    for j = 1:length(X)
        plot([X(j).savetime],[X(j).p],'o','buttondownfcn',{@HitPoint,j,X(j).p});
        hold on;
    end
    datetick;
end
if imh
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','spool','Callback',{@ImageMenu,  'spool'});
    uimenu(cmenu,'label','->FullV','Callback',{@ImageMenu,  'getfullv'});
    set(imh,'buttondownfcn',{@HitImage},'uicontextmenu',cmenu);
end
setappdata(gcf,'ClusterData',X);

function HitPoint(a,b,id,p)
X = getappdata(GetFigure(a),'ClusterData');
PrintString(X(id),p);

function HitImage(a,b)
ax = gca;
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
id = round(xy(1,2));
p = round(xy(1,1));
X = getappdata(GetFigure(a),'ClusterData');
PrintString(X(id),p);

function PrintString(X,p)

fprintf('%s on %s %d/%d spks. Drop %.1f Mahal 1D %.2f 2D%.1f. Fit %.1f,Space %s\n',...
    X.user,datestr(X.savetimes(p)),X.ncut(p),X.nspks(p),...
    X.dropi(p,3), X.mahal(p,4), X.mahal(p,1), X.fitdprime(p,1),...
    sprintf('%d ',X.space(p,:)));

function ImageMenu(a,b, fcn)
ax = gca;
xy = get(ax,'currentpoint');
id = round(xy(1,2));
p = round(xy(1,1));
X = getappdata(GetFigure(a),'ClusterData');

if strcmp('getfullv',fcn)
    args = {};
    load(X(id).name);
    fullv = GetFullVName(strrep(X(id).name,'backup',''));
    args = {args{:} 'tchan' p};
    AllVPcs(fullv,args{:},'reapply',Clusters{p},'nocheck');
end 
