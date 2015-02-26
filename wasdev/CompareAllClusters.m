function result = CompareAllClusters(A,B,varargin)
plottype = 'overlap';
result.plotxyseq = 0;
result.showtypes = [1 2 3 4];
result.showcells = [0 1 2 3];

j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'cells',4) %only plot identified cells
        result.showtypes = [0 1 2 3 4];
        result.showcells = [1 2 3];
    elseif strncmp(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmp(varargin{j},'xyseq',5)
        result.plotxyseq = 1;
    end
    j = j+1;
end
if ischar(A)
    if strcmp(A,'setexp')
        it = findobj('type','figure','Tag','CompareClusters');
        DATA = get(it,'UserData');
        DATA.currentpoint = B;
        set(DATA.toplevel,'UserData',DATA);
        HitPoint(DATA.toplevel,[],DATA.currentpoint);
        return;
    end
end
if isfield(A,'dpxc')
    PlotClusterComp(A, plottype);
    return;
end

if isfield(A{1}{1},'prefix')
    cellfile = [A{1}{1}.prefix '/CellList.mat'];
    result.prefix{1} = A{1}{1}.prefix;
    load(cellfile);
    if exist('CellListB','var')
    result.Acellid = CellList+CellListB;
    else
    result.Acellid = CellList;
    end
else
    result.prefix{1} = '1';
end
if isfield(B{1}{1},'prefix')
    result.prefix{2} = B{1}{1}.prefix;
    cellfile = [B{1}{1}.prefix '/CellList.mat'];
    load(cellfile);
    if exist('CellListB','var')
    result.Bcellid = CellList+CellListB;
    else
    result.Bcellid = CellList;
    end
else
    result.prefix{2} = '1';
end
for j = 1:length(A)
    ea(j) = A{j}{1}.exptno;
    for k = 1:length(A{j});
            A{j}{k}.cellid = result.Acellid(j,k);
    end
end
for j = 1:length(B)
    eb(j) = B{j}{1}.exptno;
    for k = 1:length(A{j});
            B{j}{k}.cellid = result.Bcellid(j,k);
    end
end
    
for j = 1:length(A)
    bj = find(eb == ea(j));
    for k = 1:length(A{j})
        res(j,k) = CompareClusterPair(A{j}{k},B{bj}{k});
    end
end
f = fields(res);
for j = 1:length(f)
    for k = 1:size(res,1)
    result.(f{j})(k,:) = cat(2,res(k,:).(f{j}));
    end
end

F = PlotClusterComp(result,plottype);
setappdata(F,'AClusters',A);
setappdata(F,'BClusters',B);

function F = PlotClusterComp(DATA, plottype, varargin)
F = GetFigure('CompareClusters');
hold off; 
colors = mycolors;
if strcmp(plottype,'xcorr')
    x = DATA.spikexc;
    y = DATA.dpxc;
else
x = DATA.overlap./DATA.nspka;
y = DATA.overlap./DATA.nspkb;
end
for j = 1:size(DATA.overlap,1)
for k = 1:size(DATA.overlap,2)
    celltype = [DATA.Acellid(j,k) > 0  DATA.Bcellid(j,k) > 0];
    celltype = celltype(1) + 2 * celltype(2);
    ptype = DATA.auto(j,k)+1;

    if ismember(DATA.auto(j,k),DATA.showtypes) && ismember(celltype,DATA.showcells)
    plot(x(j,k),y(j,k),'o','buttondownfcn',{@HitPoint, [j k]},'color', colors{ptype});
    hold on;
    end
end
end
xlabel(sprintf('common/%s',DATA.prefix{1}));
ylabel(sprintf('common/%s',DATA.prefix{2}));
DATA.toplevel = gcf;
DATA.currentpoint = [1 1];
DATA.nrows = size(DATA.overlap,1);
DATA.nprobes = size(DATA.overlap,2);
set(DATA.toplevel,'UserData',DATA);
set(gcf, 'KeyPressFcn',{@KeyPressed,3});

function HitPoint(a,b,pos);
    
    showxyseq = 0;  %takes time
DATA = GetDataFromFig(a);
e = pos(1);
p = pos(2);
fprintf('E%.0fP%.0f  xc %.2f dpxc %.2f N %d %d %d\n',e,p,DATA.dpxc(e,p),DATA.spikexc(e,p),DATA.overlap(e,p),DATA.nspka(e,p),DATA.nspkb(e,p));
fprintf('%s %s\n',datestr(DATA.ctimea(e,p)),datestr(DATA.ctimeb(e,p)));
A = getappdata(DATA.toplevel,'AClusters');
B = getappdata(DATA.toplevel,'BClusters');
GetFigure('MeanSpike1');
PlotMeanSpike(A{e}{p});
if isfield(A{1}{1},'prefix');
set(gcf,'name',A{1}{1}.prefix);
end
GetFigure('MeanSpike2');
PlotMeanSpike(B{e}{p});
if isfield(B{1}{1},'prefix');
set(gcf,'name',B{1}{1}.prefix);
end
GetFigure('XY2');
PlotClusterXY(B{e}{p})
if isfield(B{1}{1},'prefix');
set(gcf,'name',B{1}{1}.prefix);
end
GetFigure('XY1');
PlotClusterXY(A{e}{p})
if isfield(A{1}{1},'prefix');
set(gcf,'name',A{1}{1}.prefix);
end
DATA.currentpoint(1) = (e);
DATA.currentpoint(2) = (p);
if DATA.plotxyseq
GetFigure('XYCompare');
hold off;
aid = find(A{e}{p}.clst >  1);
bid = find(B{e}{p}.clst > 1);
at = A{e}{p}.times(aid);
bt = B{e}{p}.times(bid);

plot(A{e}{p}.times(aid),A{e}{p}.xy(aid,1),'.');
hold on;
plot(B{e}{p}.times(bid),B{e}{p}.xy(bid,1),'r.');
[a, ida, idb] = intersect(round(at.*1000),round(bt.*1000));
plot(A{e}{p}.times(aid(ida)),A{e}{p}.xy(aid(ida),1),'g.');
end

set(DATA.toplevel,'UserData',DATA);


 function KeyPressed(src, ks, fcn)

DATA = GetDataFromFig(src);

if strmatch(ks.Key,'rightarrow')
    if DATA.currentpoint(2) < DATA.nprobes
        DATA.currentpoint(2) = DATA.currentpoint(2)+1;
    elseif DATA.currentpoint(1) < DATA.nrows
        DATA.currentpoint(2) = 1;
        DATA.currentpoint(1) = DATA.currentpoint(1)+1;
    end
    HitPoint(src,ks,DATA.currentpoint);
elseif strmatch(ks.Key,'leftarrow')
    if DATA.currentpoint(2) > 1
        DATA.currentpoint(2) = DATA.currentpoint(2)-1;
    elseif DATA.currentpoint(1) > 1
        DATA.currentpoint(2) = DATA.nprobes;
        DATA.currentpoint(1) = DATA.currentpoint(1)-11;
    end
    HitPoint(src,ks,DATA.currentpoint);
elseif strmatch(ks.Key,'downarrow')
    if DATA.currentpoint(1) < DATA.nrows
        DATA.currentpoint(1) = DATA.currentpoint(1)+1;
    end
    HitPoint(src,ks,DATA.currentpoint);
elseif strmatch(ks.Key,'uparrow')
    if DATA.currentpoint(1) > 1
        DATA.currentpoint(1) = DATA.currentpoint(1)-1;
    end
    HitPoint(src,ks,DATA.currentpoint);
elseif strmatch(ks.Key,'add')
elseif strmatch(ks.Key,'subtract')
elseif strmatch(ks.Key,'space')
end


function PlotClusterPoints(C, uid, varargin)
        plotgmcid = 0;
        clnum = 1;
        j = 1;
        while j <= length(varargin)
            if strncmpi(varargin{j},'plotgmcid',8)
                plotgmcid = 1;
            end
            j = j+1;
        end
        if isempty(uid)
            uid = 1:size(C.xy,1);
        end
        plot(C.xy(uid,1),C.xy(uid,2),'.','markersize',1);
        hold on;
        if plotgmcid && isfield(C,'gmfit2d')
            cid = cluster(C.gmfit2d,C.xy);
            id = find(cid == clnum+1);
        elseif isfield(C,'clst')
            id = find(C.clst(uid) == clnum+1);
        elseif C.sign < 0
            id = find(C.xy(uid,1) <  C.crit(1));
        else
            id = find(C.xy(uid,1) >  C.crit(1));
        end
        plot(C.xy(uid(id),1),C.xy(uid(id),2),'r.','markersize',1);
        if isfield(C,'clst') & max(C.clst) > 2
            id = find(C.clst(uid) == 3);
            plot(C.xy(uid(id),1),C.xy(uid(id),2),'g.','markersize',1);
        end

function PlotClusterXY(C, varargin)
    titlemode = 0;
    plot.density = 0;
    clnum = 1;
    xyargs = {};
    SpaceTypeLabels = {'PCs',  'ADC', 'Template', 'Var-E', 'undef', 'N-D'};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plotgmcid',8)
            xyargs = {xyargs{:} varargin{j}};
        elseif strncmpi(varargin{j},'shorttitle',8)
            titlemode = 1;
        end
        j = j+1;
    end
    lincolor = 'r';
    
    hold off;
    p = C.probe(1);
    e = C.exptno;
uid = 1:length(C.xy);
    if plot.density
        lincolor = 'w';
        DensityPlot(C.xy(uid,1),C.xy(uid,2),'ynormal');
        hold on;
    else
        PlotClusterPoints(C,uid,xyargs{:});
    end
    if isfield(C,'space')
        spstr = [SpaceTypeLabels{C.space(1)} sprintf(' %d',C.space(2:end))];
    else
        spstr = '';
    end

    axis('tight');
    if isfield(C,'bestspace') && C.space(1) == 6
        if titlemode == 1
           h= title(sprintf('E%dP%d %.2f  %s %.2f',C.exptno,p,C.mahal(1),spstr,C.mahal(4)));
        else
           h= title(sprintf('P%d Ex %.0f Gm %.2f  %s (%.2f,%.2f for space %.0f)%.0f',p,C.exptno,C.mahal(1),spstr,C.mahal(4),C.bestspace(1),C.bestspace(2),C.sign));
        end
    else
        if titlemode == 1
            h=title(sprintf('E%dP%d %.2f  %s %.2f',C.exptno,p,C.mahal(1),spstr,C.mahal(4)));
        else
            h =title(sprintf('P%d Ex %.0f Gm %.2f  %s (%.2f)%.0f',p,C.exptno,C.mahal(1),spstr,C.mahal(4),C.sign));
        end
    end
    if C.shape == 1
        line([C.crit(1) C.crit(1)],get(gca,'ylim'),'color',lincolor);
    elseif C.shape == 2
        line([C.crit(1) C.crit(1)],get(gca,'ylim'),'color',lincolor);
    elseif C.shape == 0
        C.down = 0;
        DrawEllipse(C,'color',lincolor);
    end
    if isfield(C,'next') 
        if C.next.shape == 0
            DrawEllipse(C.next,'color','g');
        end
    end
if C.auto == 1
        set(h,'color','r');
    elseif plot.density
        set(h,'color','w');
end
if C.cellid > 0
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    h = text(mean(xl),yl(2)-diff(yl)/10,sprintf('Cell%d',C.cellid),'color','b','fontweight','bold');
end

    
function h= DrawEllipse(E,varargin)

if E.shape(1) == 1
    h = DrawLine(E,varargin{:});
    return;
end
if ~isfield(E,'down')
    E.down = 0;
end
if ~isfield(E,'xyr')  || E.down 
E.xyr(1) = mean(E.pos([1 3]));
E.xyr(2) = mean(E.pos([2 4]));
E.xyr(3) = abs(diff(E.pos([1 3])))/2;
E.xyr(4) = abs(diff(E.pos([2 4])))/2;
end
a = E.xyr(3); %x radius
b = E.xyr(4);
sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);

sn = sin(E.angle);
cn = cos(E.angle);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+mean(E.xyr(1));
y = yr+mean(E.xyr(2));

if isfield(E,'h') && ishandle(E.h)
    set(E.h,'Xdata',real(x),'Ydata',real(y));
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end

