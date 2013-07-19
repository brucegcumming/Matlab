function res = PlotSpikeShapes(ms, varargin)



cells = [];
ccl = [];
tag = 'SpikeShapes';

if ischar(ms)
    if strncmp(ms,'lem',3)
        
    name = ['/bgc/data/lem/'  ms(4:end) '/' ms 'spks.mat'];
    oname = ['/bgc/data/lem/'  ms(4:end) '/' ms 'ospk.mat'];
    cname = ['/bgc/data/lem/'  ms(4:end) '/' ms '.cells.mat'];
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
    elseif strncmpi(varargin{j},'edrange',3)
        j = j+1;
        edrange = varargin{j};
    elseif strncmpi(varargin{j},'vmean',5)
        j = j+1;
        sv = ms.v;
    elseif strncmpi(varargin{j},'waves',5)
        j = j+1;
        showwaves = varargin{j};
    end
    j = j+1;
end

GetFigure(tag);

if isempty(cr)
    cr = [min(sv(:)) max(sv(:))];
end
res.cr = cr;
[nr,nc] = Nsubplots(np);
for j = 1:np
    subplot(1,np+1,j);
    hold off;
    if ~showauto
    id = find(ms.autocut(j,:) == 1);
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
end

if isfield(ms,'Trials')
    nex = size(ms.Trials,1);
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
                ei = ei(1:end-1);
                subplot(1,np+1,p);
                hold on;
                plot(31+ones(size(cid)),et(cid),'o','color',colors{j});
                exs = unique(floor(et(cid)));
                text(32+cl*3,exs(1),num2str(j));
                text(32+cl*3,exs(end)-1,num2str(j));
                res.shapes(j,p,:) = mean(ms.Cluster(ei,p,:),1);
                res.ncshapes(j,p,:) = mean(ms.NotCluster(ei,p,:),1);
            end
        end
    end
end
subplot(1,np+1,np+1);
imagesc(ms.dprimes(:,eid)');
title(sprintf('%.1f-%.1f',min(min(ms.dprimes(:,eid))),max(max(ms.dprimes(:,eid)))));
caxis(dprange);

if isfield(ms,'Header')
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
    text(9,length(ms.eds)+1,sprintf('%.2fmm',range(ms.eds)));
    hold on;
    
    plot(depths,[1:length(depths)],'w-','linewidth',2);
    set(gca,'xticklabel',{});

    text(9,1,sprintf('%.2f',ms.eds(1)));
    h = text(1,nex+2,sprintf('%.2f',edrange(1)));
    set(h,'HorizontalAlignment','center');
    h = text(np,nex+2,sprintf('%.2f',edrange(2)));
    set(h,'HorizontalAlignment','center');
    
id = find(abs(diff(ms.eds))>0);
for j = 1:length(id)
    text(9,id(j)+1,sprintf('%.2f',ms.eds(id(j)+1)));
end
end


if length(showwaves) && length(cells)
    GetFigure('WaveForms');
    hold off;
    colors = mycolors;
    for j = 1:length(showwaves)
        id = find(~isnan(mean(res.shapes(showwaves(j),:,:),3)));
        for k = 1:length(id)
        h(j) = plot(squeeze(res.shapes(showwaves(j),id(k),:))','color',colors{j});
        hold on;
        h(j) = plot(squeeze(res.ncshapes(showwaves(j),id(k),:))',':','color',colors{j});
        end
    end
    legend(h,num2str(showwaves'));
elseif length(showwaves)
    fprintf('Need to give cell list to plot waveforms');
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

        
