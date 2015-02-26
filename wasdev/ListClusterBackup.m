function [X,D] = ListClusterBackup(prefix,varargin)
%ListClusterBackup(datadir)
%Finds cluster backup files and compares them with current definitions
%Default plot is a summary for allexpts
%ListClusterBackup(datadir/backup,'expts', e) makes a plot showing all
%backup files for a given expt. Without setting e, it will attemp to 
% ALL files, which can take a while.

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
backupdir = [prefix '/backup'];
D.prefix = prefix;
D.exptlist = exptlist;
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
        for c = 1:length(Clusters{b}.next)
            X(nb).xcluster(c) = isfield(Clusters{b}.next{c},'space');
        end
        X(nb).ncl =  max(find(X(nb).xcluster));
        X(nb).savecl = b;
        X(nb).expt = median(exs);
        X(nb).name = d(did(j)).name;
        X(nb).Clusters = Clusters;
        if ~isempty(strfind('Auto',X(nb).name));
            X(nb).auto = 1;
        else
            X(nb).auto = 0;
        end
    end
%    fprintf('%s %d %d %s\n',d(did(j)).filename,b,Clusters{b}.manual,datestr(a));
end

if isdir(backupdir)
    expts = unique([X.expt]);
    probes = 1:max([X.p]);
    for e = 1:length(expts)
        id = find([X.expt] == expts(e) & [X.auto] == 0);
        if isempty(id)
            id = find([X.expt] == expts(e));
        end
        d = mydir([backupdir '/Expt' num2str(expts(e)) 'ClusterTimes*.mat']);
        [a,b] = max([d.datenum]);
        load(d(b).name);
        Expts(expts(e)) = X(id(1));
        if exist('Clusters')
            Backups(expts(e)).Clusters = Clusters;
        end
    end
    PlotBackupExpt(Expts, Backups);
    D.toplevel = gcf;
    D.Expts = Expts;
    D.Backups = Backups;
    set(gcf,'UserData',D);
    return;
end
[a,b] = sort([X.savetime]);
X = X(b);
PlotBackup(X,D,plotargs{:});

function PlotBackupExpt(E, B,varargin);

if ishandle(E)
    F = GetFigure(E);
    E = getappdata(F,'Expts');
    B = getappdata(F,'Backups');
end
plottype = 'shape';

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'savetime',8)
        plottype = 'savetime';
    end
    j = j+1;
end
for j = 1:length(E)
    e = E(j).expt;
    for c = 1:length(B(j).Clusters);
        eC = E(j).Clusters{c};
        bC = B(j).Clusters{c};
        nspks(e,c,1) = eC.nspks;
        nspks(e,c,2) = bC.nspks;
        ncut(e,c,1) = eC.ncut;
        ncut(e,c,2) = bC.ncut;
        xc = corrcoef(eC.MeanSpike.ms,bC.MeanSpike.ms);
        xcorrs(e,c) = xc(1,2);
        savetimes(e,c*2) = eC.savetime(end);
        savetimes(e,c*2-1) = bC.savetime(end);
    end
end

[F, isnew] = GetFigure('ShapeCorr');
if isnew
   hm = uimenu(F,'label','Plots');
   uimenu(hm,'label','Default','callback',{@PlotBackupExpt, 'mahal'});
   uimenu(hm,'label','SaveTime','callback',{@PlotBackupExpt, 'savetime'});
   uimenu(hm,'label','SaveTime','callback',{@PlotBackupExpt, 'nclusters'});
end

if strcmp(plottype,'savetime')
    subplot(1,1,1);
    h(1) = imagesc(1:size(savetimes,1),[1 size(savetimes,2)/2],savetimes);
    title(sprintf('Dates %s - %s',datestr(min(savetimes(:))),datestr(max(savetimes(:)))));
else
    subplot(2,1,1);
    h(1) = imagesc(xcorrs);
    set(h,'ButtonDownFcn',@HitExImage);
    colorbar;
    title('Ratio of Event Counts Backup/Current');
subplot(2,1,2);
hold off;
Im = squeeze(ncut(:,:,1))./squeeze(ncut(:,:,2));
h(2) = imagesc(Im);
set(h,'ButtonDownFcn',@HitExImage);
colorbar;
title('Ratio of N spikes classified in Cluster 1')
end
xlabel('Probe');
ylabel('Expt');
cmenu = uicontextmenu;
hm = uimenu(cmenu,'label','AllVPcs Current','callback',{@CallAllVPcs, 'current'});
hm = uimenu(cmenu,'label','AllVPcs Previous','callback',{@CallAllVPcs, 'prev'});
hm = uimenu(cmenu,'label','See all For this Expt','callback',{@CallAllVPcs, 'fullexpt'});
set(h,'uicontextmenu',cmenu);
setappdata(F,'Expts',E);
setappdata(F,'Backups',B);


function HitExImage(a,b)

pos = get(gca,'currentpoint');
p = round(pos(1,1));
e = round(pos(1,2));
smode = get(gcf,'SelectionType');

D = GetDataFromFig(a);
fprintf('E%dP%d %d/%d,%d/%d spikes\n',e,p,D.Expts(e).ncut(p),D.Expts(e).nspks(p),...
    D.Backups(e).Clusters{p}.ncut,D.Backups(e).Clusters{p}.nspks);
Ca = D.Expts(e).Clusters{p};
Cb = D.Backups(e).Clusters{p};
if strcmp(smode,'normal');
    str = sprintf('E%dP%d',e,p);
    CompareClusters(Ca,Cb,'quiet','labels',{str, 'Current' 'Backup'});
end


function CallAllVPcs(a,b, fcn)
pos = get(gca,'currentpoint');
p = round(pos(1,1));
e = round(pos(1,2));
D = GetDataFromFig(a);
name = [D.prefix '/Expt' num2str(e) 'Spikes.mat'];
if strcmp(fcn,'fullexpt')
    ListClusterBackup([D.prefix '/backup'],'expts', e);
    return;
end
if strcmp(fcn,'current')
    C = D.Expts(e).Clusters{p};
elseif strcmp(fcn, 'previous')
    C = D.Backups(e).Clusters{p};
end
AllVPcs(name,'tchan',p,'reapply',C);

function D = PlotBackup(X, varargin)
j = 1;
plottype = 'saves';

if ishandle(X)
    F = GetFigure(X);
    X = getappdata(F,'ClusterData');
end
while j <= length(varargin)
    if isfield(varargin{j},'prefix')
        D = varargin{j};
    elseif sum(strncmpi(varargin{j},{'recluster' 'savetimes' 'fitdprime'},6))
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
   uimenu(hm,'label','SaveTime scatter','callback',{@PlotBackup, 'savescatter'});
   uimenu(hm,'label','Recluster type','callback',{@PlotBackup, 'recluster'});
   uimenu(hm,'label','N events','callback',{@PlotBackup, 'nspks'});
   uimenu(hm,'label','N Spikes Cut','callback',{@PlotBackup, 'ncut'});
   uimenu(hm,'label','Marks','callback',{@PlotBackup, 'marked'});
   uimenu(hm,'label','date/probe','callback',{@PlotBackup, 'saves'});
else
    D = get(F,'UserData');
    D.savetime = [X.savetime];
end
   setappdata(F,'ClusterData',X);
imh = 0;
hold off;
xlb = [];
ylb = [];
tstr = [];
if sum(strcmp(plottype,{'marked' 'savetimes' 'nspks' 'ncut' 'recluster'}))
    f = plottype;
    for j = 1:length(X)
        marked(j,1:length(X(j).(f))) = X(j).(f);
    end
    imh = imagesc(marked);
    D.(f) = marked;
    tstr = plottype;
    xlb = 'Probe';
    ylb = plottype;
    ylb = 'Backup File';
elseif sum(strcmp(plottype,{'space' 'fitdprime' 'mahal' 'dropi'}))
    f = plottype;
    ix = 1;
    tstr = plottype;
    xlb = 'Probe';
    ylb = 'Backup File';
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
    tstr = 'Time of save';
    ylb = 'Probe';
    xlb = 'Save Time';
    for j = 1:length(X)
        plot([X(j).savetime],[X(j).p],'o','buttondownfcn',{@HitPoint,j,X(j).p});
        hold on;
    end
    datetick;
end
title(sprintf('Ex%d %s',D.exptlist,tstr));
if ~isempty(ylb)
    ylabel(ylb);
end
if ~isempty(xlb)
    xlabel(xlb);
end

if imh
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','spool','Callback',{@ImageMenu,  'spool'});
    uimenu(cmenu,'label','->FullV','Callback',{@ImageMenu,  'getfullv'});
    set(imh,'buttondownfcn',{@HitImage},'uicontextmenu',cmenu);
end
setappdata(gcf,'ClusterData',X);
set(gcf,'UserData',D);

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

if isfield(X,'Clusters')
if isfield(X.Clusters{p},'hostname')
    hst = [X.Clusters{p}.hostname ' '];
else
    hst = [];
end
end
fprintf('%s %son %s file %s \n%d/%d spks. Drop %.1f Mahal 1D %.2f 2D%.1f. Fit %.1f,Space %s\n',...
    X.user,hst,datestr(X.savetimes(p)),X.name,X.ncut(p),X.nspks(p),...
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
