function res = PlotSpikeShapes(ms, varargin)

%plots mean spike shape for each Expt/Probe in psuedocolor.
%PlotSpikeShapes(name)
% plots mean spike shape for each Expt/Probe in psuedocolor.
% name can be a MeanSpike structure, or a 
% prefix, like 'lemM073', in which case it works out the full path
% 
%PlotSpikeShapes(name,'waves',cells) plots the mean waveform for each cell#
%in the vector cells. The mean of all waveforms for each cell at any given
%probe is shown, with each probe in a different linestyle
%
%PlotSpikeShapes(name,'labelexp') shows expt names

cells = [];
ccl = [];
ExptList = [];
tag = 'SpikeShapes';
remakelist = 0;
showdprimeline = 1;

if ischar(ms)
    if strfind(ms,'.mat') %% named the file
        name = ms;
    cname = regexprep(name,'[so][ps][kp][sk].mat','.cells.mat');
    cellfile = strrep(cname,'cells.mat','cellexp.mat');
        
    elseif strncmp(ms,'lem',3)
        
    name = ['/bgc/data/lem/'  ms(4:end) '/' ms 'spks.mat'];
    oname = ['/bgc/data/lem/'  ms(4:end) '/' ms 'ospk.mat'];
    cname = ['/bgc/data/lem/'  ms(4:end) '/' ms '.cells.mat'];
    idxfile = ['/bgc/data/lem/'  ms(4:end) '/' ms 'idx.mat'];
    cellfile = ['/bgc/data/lem/'  ms(4:end) '/' ms '.cellexp.mat'];
    end
    if exist(cname,'file');
        load(cname);
        cells = CellList;
        ccl = CellListCluster;
    end
    if exist(name,'file');
    load(name);
    ms = MeanSpike;
    elseif exist(oname,'file');
    load(oname);
    ms = MeanSpike;
    elseif length(cells)
        x.CellList = cells;
        x.CellListCluster = ccl;
        x.CellQuality = CellQuality;
        x.plot.cellplot = 1;
        x.probelist = 1:8;
        PlotCellList(x,'quality');
        return;
    else
        return;
    end
elseif isnumeric(ms)
    fign = findobj('Tag',tag,'Type','Figure');
    if length(fign)
        figure(fign(1));
    else
        return;
    end
    DATA = get(gcf,'UserData');
    if isfield(DATA,'hlines')
        for j = 1:length(DATA.hlines)
            if ishandle(DATA.hlines(j))
                delete(DATA.hlines(j));
            end
        end
        DATA.hlines = [];
    end
    yl = get(gca,'Ylim');
    ap = get(gca,'position');
    for j = 1:length(ms)
        yp = ap(2)+ap(4) - (ap(4)) * (ms(j)-yl(1))/range(yl);
        if yp > 0 && yp < 1
        DATA.hlines(j) = annotation('line',[0.1 0.9],[yp yp]);
        set(DATA.hlines(j),'color','r');
        end
    end
    set(gcf,'UserData',DATA);
    cellfile = [];
    return;
else
    cellfile = [];
end



np = size(ms.v,2);
eid = 1:size(ms.v,1);
showauto = 0;
mindp = 0;
sv = ms.Cluster;
cr = [];
dprange = [1 10];
showwaves = [];
edrange = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'auto',4)
        showauto = 1;
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        eid = varargin{j};
    elseif strncmpi(varargin{j},'crange',2)
        j = j+1;
        cr = varargin{j};
    elseif strncmpi(varargin{j},'cells',5)
        j = j+1;
        cells = varargin{j};
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            ccl = varargin{j};
        end
    elseif strncmpi(varargin{j},'cldp',4)
        sv = (ms.Cluster-ms.NotCluster)./ms.ClusterSD;
        sv(find(ms.ClusterSD==0)) = 0;
    elseif strncmpi(varargin{j},'cldiff',5)
        sv = (ms.Cluster-ms.NotCluster);
        sv(find(ms.ClusterSD==0)) = 0;
    elseif strncmpi(varargin{j},'dprange',5)
        j = j+1;
        dprange = varargin{j};
    elseif strncmpi(varargin{j},'dprime',5)
        j = j+1;
        mindp = varargin{j};
    elseif strncmpi(varargin{j},'labelexp',8)
        idxfile = strrep(cname,'.cells.mat','idx.mat');
        if exist(idxfile,'file')
            load(idxfile);
        end
    elseif strncmpi(varargin{j},'edrange',3)
        j = j+1;
        edrange = varargin{j};
    elseif strncmpi(varargin{j},'remake',5)
        remakelist = 1;
    elseif strncmpi(varargin{j},'vmean',5)
        j = j+1;
        sv = ms.v;
    elseif strncmpi(varargin{j},'waves',5)
        j = j+1;
        showwaves = varargin{j};
    end
    j = j+1;
end

[sfig, isnew] = GetFigure(tag);
if isnew
    DATA.cellid = 0;
    DATA.toplevel = sfig;
    DATA.MeanSpike = ms;
    DATA.ExptList = ExptList;
    DATA.cellexpts{1} = [];
    DATA.cellexptprobes{1} = [];
    DATA.cellfile = cellfile;
    if exist(DATA.cellfile,'file') && ~remakelist
        load(DATA.cellfile);
        DATA.cellexpts = CellExpts;
        DATA.cellexptprobes = CellExptProbes;
    end
    bp = [0 0 0.05 0.05];
    uicontrol(sfig,'style','pop','string',num2str([0:24]'),'Units','normalized','Position',bp,'Tag','Cellid',...
        'Callback',@Update);
    set(DATA.toplevel,'UserData',DATA);
else
    DATA = get(sfig,'UserData');
end
set(gcf, 'WindowButtonDownFcn',@ButtonPressed);
set(gcf, 'WindowButtonMotionFcn',@ButtonDragged);
set(gcf, 'WindowButtonUpFcn',@ButtonReleased);



if isempty(cr)
    cr = [min(sv(:)) max(sv(:))];
end
res.cr = cr;
[nr,nc] = Nsubplots(np);
for j = 1:np
    subplot(1,np+1,j);
    hold off;
    if ~showauto
    id = find(ms.autocut(j,:) >= 1);
    sv(id,j,:) = NaN;
    end
    id = find(ms.dprimes(j,:) < mindp);
    sv(id,j,:) = 0;
    ms.dprimes(j,id) = 0;
    imagesc(squeeze(sv(eid,j,:)));
    hold on;
    id = find(ms.autocut(j,:) == 0);
    if length(id)
        plot([0 10],[id; id;],'w-');
    end
    hold off;
    if j > 1
        set(gca,'ytick',[]);
    end
    caxis(cr);
    title(num2str(j));
    sub.probe = j;
    set(gca,'UserData',sub);
end

colors = mycolors;
for j = 1:length(DATA.cellexpts)
    probes = unique(DATA.cellexptprobes{j});
    for p = 1:length(probes)
        subplot(1,np+1,probes(p));
        id = find(DATA.cellexptprobes{j} == probes(p));
        ShowCellRect(j, probes(p), min(DATA.cellexpts{j}(id)), max(DATA.cellexpts{j}(id)), colors{j});
    end
end

if ~isempty(ExptList)
    DATA.ExptList = ExptList;
end
if ~isempty(DATA.ExptList)
    ExptList = DATA.ExptList;
    subplot(1,np+1,1);
    for j = 1:length(ExptList)
        if j ==1 | ~strcmp(ExptList(j).expname,ExptList(j-1).expname)
        DATA.extext(j) = text(-10,j,sprintf('%d: %s',j,ExptList(j).expname),'HorizontalAlignment','Right');
        end
        
    end
end


if isfield(ms,'Trials')
    nex = size(ms.Trials,1);
else
    nex = 0;
end    
if ~isempty(cells) & isfield(ms,'Trials')
    colors = mycolors;
    for j = 1:nex
        et(ms.Trials(j,1):ms.Trials(j,2)) = j:1./diff(ms.Trials(j,:)):j+1;
        if j < nex
        et(ms.Trials(j,2):ms.Trials(j+1,1)) = j+1;
        end
    end
    res.shapes = ones([size(cells,1) np 32]) .*NaN;
    for j = 1:size(cells,1)
        for p = 1:np
            cid = find(cells(j,:) == p);
            if isempty(ccl)
                cl = 1;
            else
                cl = median(ccl(j,cid));
            end
            if length(cid) %this cell defined for this expt
                a = [cid(1) cid(end)];
                ei = unique(floor(et(cid)));
                if length(ei) > 1
                ei = ei(1:end-1);
                end
                cid = cid(find(et(cid) < ei(end)+0.5));
                if length(cid)
                subplot(1,np+1,p);
                hold on;
                plot(31+ones(size(cid)),et(cid),'o','color',colors{j});
                exs = unique(floor(et(cid)));
                text(32+cl*3,exs(1),num2str(j));
                text(32+cl*3,exs(end)-1,num2str(j));
                ei = ei(~isnan(squeeze(ms.Cluster(ei,p,1))));
                res.shapes(j,p,:) = mean(ms.Cluster(ei,p,:),1);
                res.ncshapes(j,p,:) = mean(ms.NotCluster(ei,p,:),1);
                res.eds{j,p} = ms.eds(ei);
                res.eis{j,p} = ei;
                end
            end
        end
    end
end
subplot(1,np+1,np+1);
hold off;
imagesc(ms.dprimes(:,eid)');
title(sprintf('%.1f-%.1f',min(min(ms.dprimes(:,eid))),max(max(ms.dprimes(:,eid)))));
caxis(dprange);

if isfield(ms,'Header') & isfield(ms.Header,'probesep')
   probesep = ms.Header.probesep;
else
    probesep = 150;
end
    

%subplot(1,np+1,1);
if isfield(ms,'eds')
    mdepth = median(ms.eds);
    if isempty(edrange)
        edrange = minmax(ms.eds);
    end
    depths = (edrange(2)-ms.eds)./(probesep/1000); %depth in units of probes
    erange = edrange./(probesep/1000);
    depths = np .* depths./diff(erange);
    depths = 0.95 .* (depths+1.0);
    text(np+1,length(ms.eds)+1,sprintf('%.2fmm',range(ms.eds)));
    hold on;
    
    plot(depths,[1:length(depths)],'w-','linewidth',2);
    set(gca,'xticklabel',{});

    text(np+1,1,sprintf('%.2f',ms.eds(1)));
    if nex
    h = text(1,nex+2,sprintf('%.2f',edrange(1)));
    set(h,'HorizontalAlignment','center');
    h = text(np,nex+2,sprintf('%.2f',edrange(2)));
    set(h,'HorizontalAlignment','center');
    end
    
id = find(abs(diff(ms.eds))>0);
for j = 1:length(id)
    text(np+1,id(j)+1,sprintf('%.2f',ms.eds(id(j)+1)));
end
end

linestyles{1} = '-';
linestyles{2} = '--';
linestyles{3} = '-.';
linestyles{4} = ':';


shownc = 0;
wavetype = 1;
if length(showwaves) && length(cells)
    GetFigure('WaveForms');
    if wavetype == 0
    hold off;
    colors = mycolors;
    for j = 1:length(showwaves)
        id = find(~isnan(mean(res.shapes(showwaves(j),:,:),3)));
        for k = 1:length(id)
        h(j) = plot(squeeze(res.shapes(showwaves(j),id(k),:))','color',colors{j},'linewidth',2,'linestyle',linestyles{k});
        hold on;
        if shownc
        h(j) = plot(squeeze(res.ncshapes(showwaves(j),id(k),:))',':','color',colors{j});
        end
        end
    end
    legend(h,num2str(showwaves'));
    elseif wavetype == 1
    hold off;
    colors = mycolors;
    x = ms.probesep.*[1:size(res.shapes,3)]./(size(res.shapes,3) * 1000);
    dx = x(2);
    xo = 0;
    for j = 1:length(showwaves)
        id = find(~isnan(mean(res.shapes(showwaves(j),:,:),3)));
        for k = 1:length(id)
            ed = mean(res.eds{showwaves(j),id(k)});
            ei = mean(res.eis{showwaves(j),id(k)});
            xo  = (id(k).*ms.probesep/1000) + mean(res.eds{showwaves(j),id(k)});
            yo = 5 * ei/length(ms.eds);  % all expts = 5V range
            wave = squeeze(res.shapes(showwaves(j),id(k),:))';
        h(j) = plot(x+xo,wave+yo,'color',colors{j},'linewidth',2,'linestyle',linestyles{k});
        hold on;
         [a,b] = min(wave);
        t = text(x(b)+xo+dx,yo,sprintf('%d:%.2f',id(k),ed) ,'color',colors{j});
        hn(j) = plot(x+xo,yo+squeeze(res.ncshapes(showwaves(j),id(k),:))',':','color',colors{j});
        end
        if isempty(id)
            good(j) = 0;
        else
            good(j) = 1;
        end
    end
    legend(h(find(good)),num2str(showwaves(find(good))'));
        
    end
elseif length(showwaves)
    GetFigure('WaveForms');
    hold off;
    step = 0.25;
    subplot(1,1,1);
    np = size(ms.Cluster,3);
    for ip = 1:length(showwaves);
        p = showwaves(ip);
        for j = 1:size(ms.Cluster,1)
            if ~isnan(sum(ms.Cluster(j,p,:))) && (ms.autocut(p,j) == 0)
                plot([1:np]+ip*np,squeeze(ms.Cluster(j,p,:))-j*step);
                hold on;
            end
        end
        yl(ip,:) = get(gca,'ylim');
        h = text(ip*np+np/3,1,num2str(p));
        set(h,'color','r');
    end
end     
if showdprimeline
GetFigure('Dprimes')
ng = ceil(np/8);
for j = 1:ng
    subplot(ng,1,j);
    ps = [1:8] + 8 * (j-1);
    plot(ms.dprimes(ps,:)','linewidth',2);
    legend(num2str(ps'));
end
end


function PlotCellList(DATA, varargin)
    
QUALITY = 1;
CELLNUM=2;
ONECELL=3;
JOINCELL=4;
CLUSTERQUALITY=5;
cell = 1;

if isfield(DATA.plot,'cellplot')
    plottype = DATA.plot.cellplot;
else
plottype = QUALITY;
end

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
    end
    j = j+1;
end

if ~isfield(DATA,'CellList') | isempty(DATA.CellList)
    return;
end
it = findobj('Tag','CombinerCellList');
it = findobj(it, 'Tag','CellNumber');
if ~isempty(it)
    cell = get(it(1),'value');
end

if plottype == QUALITY
    if(size(DATA.CellQuality,2) < size(DATA.CellList,2))
        DATA.CellQuality(end,size(DATA.CellList,2)) = 0;
    end
    im = zeros(length(DATA.probelist)*2,size(DATA.CellList,2));
    for j = 1:size(DATA.CellList,1)
        id = find(DATA.CellList(j,:) > 0);
        k = DATA.CellList(j,id)*2;
        ind = sub2ind(size(im),k,id);
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
        
        if plottype == QUALITY & isfield(DATA,'CellQuality')
            im(ind) = DATA.CellQuality(j,id);
            im(indb) = DATA.CellQuality(j,[aid bid]);
        else
        im(ind) = j;
        im(indb) = j;
        end
    end
    hold off;
    imagesc([1:size(DATA.CellList,2)],[1:8],im);
    hold on;
    DATA.CellList(find(DATA.CellList == 0)) = NaN;
    cl = (DATA.CellListCluster-1.5)/2;

    plot((DATA.CellList')+cl','linewidth',2);
    legend(num2str([1:size(DATA.CellList,1)]'),'location','NorthWest')
elseif plottype ==CLUSTERQUALITY
    Q = zeros(length(DATA.probelist),DATA.Expts{end}.Trials(end).Trial);
    for e = 1:length(DATA.Expts)
        id = [DATA.Expts{e}.Trials.Trial];
        if isfield(DATA.Expts{e}.Trials,'ed')
            eds(id) = [DATA.Expts{e}.Trials.ed];
        else
            eds(id) = DATA.Expts{e}.Stimvals.ed;
        end
        if isfield(DATA.Expts{e},'Cluster')
            for k = 1:size(DATA.Expts{e}.Cluster,2)
                if isfield(DATA.Expts{e}.Cluster{1,k},'quality')
                    Q(k,id) = DATA.Expts{e}.Cluster{1,k}.quality;
                end
            end
        end
    end
    hold off;
    imagesc(Q);
    hold on;
    cl = (DATA.CellListCluster-1.5)/2;

    plot((DATA.CellList')+cl','linewidth',2);
    legend(num2str([1:size(DATA.CellList,1)]'),'location','NorthWest')
    [a,b] = TrialRange(DATA);
    plot([a a],[1 length(DATA.probelist)],'w:');
    plot([b b],[1 length(DATA.probelist)],'w:');
    id = find(eds > 0);
    plot(id, (max(eds(id))-eds(id))./0.15,'w');
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
    [a,b] = TrialRange(DATA);
    plot([a a],[1 length(DATA.probelist)],'k:');
    plot([b b],[1 length(DATA.probelist)],'k:');
        
    legend(num2str([1:size(DATA.CellList,1)]'),'location','SouthWest')
    hold off;
    set(gca,'ylim',[min(DATA.CellList(:))-1 max(DATA.CellList(:))+1]);
end

        
function ButtonPressed(src, data)

DATA = get(gcf,'UserData');
pt = get(gca,'CurrentPoint');
DATA.mouse.start = pt;
set(gcf,'UserData',DATA);

function ButtonReleased(src, data)

DATA = get(gcf,'UserData');
pt = get(gca,'CurrentPoint');
DATA.mouse.end = pt;
set(gcf,'UserData',DATA);

a = get(gca,'UserData');
probe = a.probe;
first = ceil(DATA.mouse.start(1,2)-0.5);
last = floor(DATA.mouse.end(1,2)+0.5);
x = [1 first-0.25 30 0.5+last-first];
hold on;
colors = mycolors;
if DATA.cellid > 0
    cellid = DATA.cellid;
    ShowCellRect(cellid, probe, first, last, colors{cellid});
    h = rectangle('Position',x);
    set(h,'Edgecolor',colors{DATA.cellid});
    h = text(34,(last+first)/2,num2str(DATA.cellid));
    set(h,'color', colors{DATA.cellid},'FontWeight','bold');
    if cellid >length(DATA.cellexpts)
        DATA.cellexpts{cellid} = [];
        DATA.cellexptprobes{cellid} = [];
    end
    DATA.cellexpts{cellid} = [DATA.cellexpts{cellid} first:last];
    DATA.cellexptprobes{cellid} = [DATA.cellexptprobes{cellid} sign(first:last) * probe];
    set(gcf,'UserData',DATA);
    CellExpts = DATA.cellexpts;
    CellExptProbes = DATA.cellexptprobes;
    ExptList = DATA.ExptList;
    save(DATA.cellfile,'CellExpts','CellExptProbes','ExptList');
end


function ShowCellRect(cellid, probe, first, last, color)

x = [1 first-0.25 30 0.5+last-first];
if cellid > 0
    h = rectangle('Position',x);
    set(h,'Edgecolor',color);
    h = text(34,(last+first)/2,num2str(cellid));
    set(h,'color', color,'FontWeight','bold');
end

function ButtonDragged(src, data)

DATA = get(gcf,'UserData');
pt = get(gca,'CurrentPoint');
set(gcf,'UserData',DATA);


function Update(src,data)

DATA = get(get(src,'parent'),'UserData');
DATA.cellid = get(findobj(DATA.toplevel,'Tag','Cellid'),'value')-1;
set(get(src,'parent'),'UserData',DATA);
