function [X,D] = ListClusterBackup(prefix,varargin)
%LISTCLUSTERBACKUP ListClusterBackup(prefix)
%Finds cluster backup files and compares them with current definitions
%Default plot is a summary for allexpts
%ListClusterBackup(datadir/backup,'expts', e) makes a plot showing all
%backup files for a given expt. Without setting e, it will attemp to 
% ALL files, which can take a while.

j = 1;
exptlist = [];
missing = [];
plotargs = {};
D = [];

while j <= length(varargin)
    if strncmpi(varargin{j},'clusters',7)
        plotargs = {plotargs{:} varargin{j} varargin{j+1}};
        j = j+1;
        EClusters = varargin{j};
    elseif strncmpi(varargin{j},'expt',4)
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
            X(nb).ncl =  max(find(X(nb).xcluster));
        end       
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
    tic;
    expts = unique([X.expt]);
    probes = 1:max([X.p]);
    for e = 1:length(expts)
        Backups(expts(e)).Clusters{1} = X(e).Clusters;
        id = find([X.expt] == expts(e) & [X.auto] == 0);
        if isempty(id)
            id = find([X.expt] == expts(e));
        end
        d = mydir([backupdir '/Expt' num2str(expts(e)) 'ClusterTimes*.mat']);
        [a,b] = sort([d.datenum],'descend');
        if length(expts) > 1 %load  last backup backups
            d = d(b(1)); %most recent
        end
        for j = 1:length(d)
            x = load(d(b(j)).name);
            if isfield(x,'Clusters')
                Backups(expts(e)).Clusters{end+1} = x.Clusters;
                nb = length(Backups(expts(e)).Clusters);
                Backups(expts(e)).names{nb} = d(b(j)).name;
            end
        end
        Expts(expts(e)) = X(id(1));
    end
    if length(expts) == 1
        D.oneexpt = expts;
    else
        D.oneexpt = 0;
    end
    D.Expts = Expts;
    D.Backups = Backups;
    F = PlotBackupExpt(Expts, Backups, D);
    D = get(F,'UserData');
    toc;
    return;
end

[a,b] = sort([X.savetime]);
X = X(b);
PlotBackup(X,D,plotargs{:});

function F = PlotBackupExpt(E, B,varargin);

if ishandle(E)
    F = GetFigure(E);
    E = getappdata(F,'Expts');
    B = getappdata(F,'Backups');
end
plottype = 'shape';

j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'oneexpt')
        DATA = varargin{j};
    elseif strncmpi(varargin{j},'checkbreaks',8)
        plottype = 'checkbreaks';
        CheckBreaks(E,B);
    elseif strncmpi(varargin{j},'savetime',8)
        plottype = 'savetime';
    end
    j = j+1;
end


for j = 1:length(E)
    e = E(j).expt;
    if e > 0
        exptid(j) = e;
    else
        exptid(j) = 0;
    end
    for b = 1:length(B(j).Clusters)
    for c = 1:length(B(j).Clusters{b});
        eC = E(j).Clusters{c};
        bC = B(j).Clusters{b}{c};
        nspks(e,c,1) = eC.nspks;
        nspks(e,c,1+b) = bC.nspks;
        ncut(e,c,1) = eC.ncut;
        ncut(e,c,1+b) = bC.ncut;
        xc = corrcoef(eC.MeanSpike.ms,bC.MeanSpike.ms);
        xcorrs(e,c,b) = xc(1,2);
        savetimes(e,c,1) = eC.savetime(end);
        savetimes(e,c,1+b) = bC.savetime(end);
    end
    end
end
expts = exptid(exptid > 0);

[F, isnew] = GetFigure('ShapeCorr');
if isnew
   hm = uimenu(F,'label','Plots');
   uimenu(hm,'label','Default','callback',{@PlotBackupExpt, 'mahal'});
   uimenu(hm,'label','SaveTime','callback',{@PlotBackupExpt, 'savetime'});
   uimenu(hm,'label','# of clusters','callback',{@PlotBackupExpt, 'nclusters'});
   uimenu(hm,'label','Check Breaks','callback',{@PlotBackupExpt, 'checkbreaks'});
   D = DATA;
   D.toplevel = F;
else
    D = GetDataFromFig(F);
end

if length(expts) == 1
    D.oneexpt = e;
    str = [GetName(D.prefix) ' Expt' num2str(e)];
else
    D.onexpt = 0;
    str = GetName(D.prefix);
end
SetData(D);

if strcmp(plottype,'savetime')
    subplot(1,1,1);
    if length(expts) == 1
        h(1) = imagesc(1:size(savetimes,3),[1 size(savetimes,2)],squeeze(savetimes(e,:,:)));
        ylabel('probe');
        xlabel('File #');
    else
        h(1) = imagesc(1:size(savetimes,1),[1 size(savetimes,2)/2],savetimes);
    end
    title(sprintf('Dates %s - %s',datestr(min(savetimes(:))),datestr(max(savetimes(:)))));
else
    subplot(2,1,1);
    if length(expts) == 1
        h(1) = imagesc(squeeze(xcorrs(e,:,:)));
        ylabel('probe');
        xlabel('Backup#');
    else
        h(1) = imagesc(xcorrs);
    end
    set(h,'ButtonDownFcn',@HitExImage);
    colorbar;
    title([str ' Ratio of Event Counts Backup/Current']);
    subplot(2,1,2);
    hold off;
    if length(expts) == 1
        for j = 1:size(ncut,3)
            Im(:,j) = squeeze(ncut(e,:,j))./squeeze(ncut(e,:,end));
        end
        h(2) = imagesc(Im);
        xlabel('Backup #');
        ylabel('probe');
    else
        Im = squeeze(ncut(:,:,1))./squeeze(ncut(:,:,2));
        h(2) = imagesc(Im);
    end
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


function CheckExptBreaks(E, D)

nprobes = length(E(1).Clusters);
nexp = length(E);
for p = 1:nprobes;
for j = 1:length(E)
        if j < nexp
            eC = E(j+1).Clusters{p};
        else
            eC = D.Clusters{p};
        end
        bC = E(j).Clusters{p};
        ncut(j,p) = eC.ncut-bC.ncut;
        xc = corrcoef(eC.MeanSpike.ms,bC.MeanSpike.ms);
        xcorrs(j,p) = xc(1,2);
        effic(j,p) = 1; %to make plot easy to read
        if xcorrs(j,p) < 0.9
            [a,b] = xcorrtimes(eC.times,bC.times);
            effic(j,p) = max(b.efficacy); 
            if max(b.efficacy) < 0.7
                
            end
        elseif isnan(xcorrs(j,p))
            [a,b] = xcorrtimes(eC.times,bC.times);
            effic(j,p) = max(b.efficacy);             
        end
    end
end
hold off; 
imh = imagesc(xcorrs);
    cmenu = uicontextmenu;
    uimenu(cmenu,'label','spool','Callback',{@ImageMenu,  'spool'});
    uimenu(cmenu,'label','->FullV','Callback',{@ImageMenu,  'getfullv'});
    set(imh,'buttondownfcn',{@HitImage, 'sequence'},'uicontextmenu',cmenu);

function CheckBreaks(E,B)
for j = 1:length(E)
    e = E(j).expt;
    nc = length(B(j).Clusters);
    for c = 1:length(nc);
        if c < nc
            eC = B(j+1).Clusters{c};
        else
            eC = B(j).Clusters{c};
        end
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



function HitExImage(a,b)

pos = get(gca,'currentpoint');
p = round(pos(1,1));
e = round(pos(1,2));
smode = get(gcf,'SelectionType');

D = GetDataFromFig(a);
if D.oneexpt > 0
    eid = p;
    p = e;
    e = D.oneexpt;
    nbackup = length(D.Backups(e).Clusters);
    eid = min([eid nbackup]);
    fprintf('E%dP%d Backup%d %d/%d,%d/%d spikes %s\n',e,p,eid,D.Expts(e).ncut(p),D.Expts(e).nspks(p),...
        D.Backups(e).Clusters{eid}{p}.ncut,D.Backups(e).Clusters{eid}{p}.nspks,...
        D.Backups(e).names{eid});
    Ca = D.Expts(e).Clusters{p};
    Cb = D.Backups(e).Clusters{eid}{p};
else
    fprintf('E%dP%d %d/%d,%d/%d spikes\n',e,p,D.Expts(e).ncut(p),D.Expts(e).nspks(p),...
        D.Backups(e).Clusters{p}.ncut,D.Backups(e).Clusters{p}.nspks);
    Ca = D.Expts(e).Clusters{p};
    Cb = D.Backups(e).Clusters{p};
end
if strcmp(smode,'normal');
    str = sprintf('E%dP%d',e,p);
    CompareClusters(Ca,Cb,'quiet','labels',{str, 'Current' 'Backup'});
end


function CallAllVPcs(a,b, fcn)
pos = get(gca,'currentpoint');
p = round(pos(1,1));
e = round(pos(1,2));
D = GetDataFromFig(a);
if D.oneexpt > 0
    nb = p;
    p = e;
    e = D.oneexpt;
    name = [D.prefix '/Expt' num2str(e) 'FullV.mat'];
else
    name = [D.prefix '/Expt' num2str(e) 'Spikes.mat'];
end

if strcmp(fcn,'fullexpt')
    ListClusterBackup([D.prefix '/backup'],'expts', e, 'Clusters', D.Expts(e).Clusters);
    return;
end
if strcmp(fcn,'current')
    C = D.Expts(e).Clusters{p};
elseif strcmp(fcn, 'prev')
    if D.oneexpt > 0
        C = D.Backups(e).Clusters{nb}{p};
    else
        C = D.Backups(e).Clusters{p};
    end
end
res = AllVPcs(name,'tchan',p,'reapply',C);
B = res.cluster;
CompareClusters(C,res.cluster);
ids = MatchTimes(C.ncut,B.ncut, 0.001);
fprintf('Reapply gets %d/%d events, %d/%d spikes %d match\n',res.cluster.nspks,C.nspks,res.cluster.ncut,C.ncut,length(ids));

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
    elseif strcmpi(varargin{j}, 'clusters')
        j = j+1;
        D.Clusters = varargin{j};
    elseif sum(strncmpi(varargin{j},{'checkbreak' 'recluster' 'savetimes' 'fitdprime'},6))
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
   uimenu(hm,'label','Check Breaks','callback',{@PlotBackup, 'checkbreak'});
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
if strcmp(plottype,'checkbreak')
    CheckExptBreaks(X, D);
    return;
end

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

function HitImage(a,b, varargin)

comparelast = 0;
for j = 1:length(varargin)
    if strcmp(varargin{j},'sequence')
        comparelast = 1;
    end
    j = j+1;
end

ax = gca;
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
id = round(xy(1,2));
p = round(xy(1,1));
X = getappdata(GetFigure(a),'ClusterData');
if comparelast
    if id > 1
    CompareClusters(X(id).Clusters{p},X(id-1).Clusters{p},'print')
    end
end
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
