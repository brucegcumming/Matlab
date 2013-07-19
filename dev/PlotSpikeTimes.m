function result = PlotSpikeTimes(C, varargin)
%
%Plot Spike Times from Cluster lists made by AllVPCs
%PlotSpikeTimes(Clusters,'scatter')
%PlotSpikeTimes(Clusters,'xcorr')
%PlotSpikeTimes(Clusters,'xcorrstep')
colors = mycolors;
plottype = 'dips';
pstep = 0;
result.starttime = now;
Details = [];
matrixplot = 0;

hitfcn = @PlotSpikeTimes;
j = 1;
while j <= length(varargin)
    if iscell(varargin{j}) 
        for k = 1:length(varargin{j})
            empty(k) = ~isempty(varargin{j}{k});
        end
        k = find(empty);
        if ishandle(C(1)) %called from graph
        C = varargin{j};
        elseif isfield(varargin{j}{k(1)},'xy')
            Details = varargin{j};
        end
    elseif isempty(varargin{j})
        
    elseif strmatch(varargin{j},{'dips' 'probequality' 'quality' 'mahaldip' 'bestspace' 'probesd'})
        plottype = varargin{j};
    elseif strcmpi(varargin{j},'hitfcn')
        j = j+1;
        hitfcn = varargin{j};
    elseif strcmpi(varargin{j},'hit')
        j = j+1;
        hitpt = varargin{j};
        fprintf('C%d\n',hitpt);
        return;
    elseif strcmpi(varargin{j},'matrix')
        matrixplot = 1;
    elseif strcmpi(varargin{j},'scatter')
        plottype = 'scatter';
    elseif strcmpi(varargin{j},'xcorrstep')
        plottype = varargin{j};
        j = j+1;
        pstep = varargin{j};
    elseif strncmpi(varargin{j},'xcorr',4)
        plottype = varargin{j};
    end
    j = j+1;
end

for j = 1:length(C)
    if ~isempty(C{j})
        if ~isfield(C{j},'mahal')
            C{j}.mahal = 0;
        end
        if isfield(C{j},'clst') %limit spike times to cluster
            id = find(C{j}.clst == 2);
            C{j}.times = C{j}.times(id);
        end
    end
end

if strcmp(plottype,'times')
    k = 0;
    hold off;
    for j = 1:length(C)
        if isfield(C{j},'times')
            k = k+1;
            plot(C{j}.times,k,'.','color',colors{k})
            hold on;
        end
    end
elseif strncmp(plottype,'scatter',4)  && length(Details)
    if length(Details) == 24
        nc= 6;
        nr = 4;
    else
    [nr,nc] = Nsubplots(length(Details));
    end
    subplot(1,1,1);
    for j = 1:length(Details);
        mysubplot(nr,nc,j);
        hold off; 
        if isfield(Details{j},'xy')
        plot(Details{j}.xy(:,1),Details{j}.xy(:,2),'.','markersize',1);
        if C{j}.sign < 0
            id = find(Details{j}.xy(:,1) < C{j}.crit);
        else
            id = find(Details{j}.xy(:,1) > C{j}.crit);
        end
        hold on;
        plot(Details{j}.xy(id,1),Details{j}.xy(id,2),'r.','markersize',1);
        axis('tight');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        text(xl(1),yl(2),sprintf('P%d,%.1f,%.1f',j,C{j}.hdip*100,C{j}.bmc*10),'horizontalalignment','left','verticalalignment','top','fontsize',8);
        end
    end
elseif strncmp(plottype,'probesd',6)
    for j = 1:length(C)
        sds(j,:) = std(C{j}.MeanSpike.ms,0,2);
    end
    imagesc(sds);
    result.sds = sds;
elseif strncmp(plottype,'dips',4) || strncmp(plottype,'quality',4)
    k = 0;
    hold off;
    dpscale = 3;
    [dips, scales, probes] = ClusterQualities(C); 
    for k = 1:size(dips,2)
        j = probes(k);
            plot(dips(1,k),dips(2,k).*scales(2),'o','buttondownfcn',{hitfcn, C,'hit', j});
            hold on;
            plot(dips(1,k),dips(3,k).*scales(3),'ko','buttondownfcn',{hitfcn, C,'hit', j});
            plot(dips(1,k),dips(4,k).*scales(4),'go','buttondownfcn',{hitfcn, C,'hit', j});
            plot(dips(1,k),dips(5,k).*scales(5),'co','buttondownfcn',{hitfcn, C,'hit', j});
            plot(dips(1,k),dips(6,k).*scales(6),'ro','buttondownfcn',{hitfcn, C,'hit', j});
    end
    xlabel('Hartigan');
    legend('mydip','sumdp','dprime','bmc','mahal');
    result.dips = dips;
elseif strncmp(plottype,'mahaldip',7)
    k = 0;
    hold off;
    dpscale = 3;
    [dips, scales, probes] = ClusterQualities(C); 
    for k = 1:size(dips,2)
        j = probes(k);
            plot(dips(1,k),dips(6,k),'o','buttondownfcn',{hitfcn, C,'hit', j});
            hold on;
    end
    xlabel('Hartigan');
    ylabel('mahal');
    result.dips = dips;
elseif strncmp(plottype,'bestspace',7)
%compare mahal distance for best N-D space from gmfit to final mahal
%distance
    k = 0;
    hold off;
    dpscale = 3;
    [dips, scales, probes] = ClusterQualities(C); 
    for k = 1:size(dips,2)
        j = probes(k);
            plot(dips(6,k),dips(7,k),'o','buttondownfcn',{hitfcn, C,'hit', j});
            hold on;
    end
    refline(1);
    xlabel('mahal 2D');
    ylabel('mahal ND');
    result.dips = dips;
elseif strncmp(plottype,'probequality',4)
    k = 0;
    [dips, scales] = ClusterQualities(C); 
    n= size(dips,1);
    hold off;
    dpscale = 3;
    probes = 1:size(dips,2);
    plot(probes,dips(1,:),'mo');
    hold on;
    plot(probes,dips(2,:).*scales(2),'o');
    plot(probes,dips(3,:).*scales(3),'ro');
    plot(probes,dips(4,:).*scales(4),'go');
    plot(probes,dips(5,:).*scales(5),'co');
    plot(probes,dips(6,:).*scales(6),'co');
    legend('Hartigan','mydip','dp','dprime','bmc','mahal');
    

    
elseif strncmp(plottype,'xcorr',4)
    ClearPlot;
    minbmc = 0.2;
    k = 0;
    if matrixplot
        k = length(C);
        ids = 1:k;
    else
        for j = 1:length(C)
            if isfield(C{j},'times') && GoodCluster(C{j})
                k = k+1;
                ids(k) = j;
            end
        end
    end
    result.celllist = ids;
    ncl = k;
    np = ceil(sqrt(ncl.*(ncl-1)/2));
    if strcmp(plottype,'xcorrstep')
        [nr,nc] = Nsubplots(ncl-pstep);
        for j = 1:length(ids)-pstep
            xc = xcorrtimes(C{ids(j)}.times,C{ids(j+pstep)}.times,'clip');
            k = ids(j+pstep);
            result.synci(ids(j),ids(k),:) = SyncIndex(xc);
            result.synci(ids(k),ids(j),:) = SyncIndex(xc);
            result.ccf{ids(j),ids(k)} = xc;
            mysubplot(nr,nc,j);
            plot(xc);
                set(gca,'xticklabel',[],'yticklabel',[]);
                axis('tight');
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                h = text(mean(xl),yl(2),sprintf('%d,%d,%.1f',ids(j),ids(j+pstep),result.synci(ids(j),ids(k),1)),'verticalalignment','top');
                if GoodCluster(C{ids(j)}) > 1
                    set(h,'color','r');
                end
                drawnow;
        end
    else
        n = 0;
        nc = length(ids);
        for j = 2:length(ids)
            for k = 1:j-1
                n = n+1;
                result.times(n) = now;
                tic;
                [xc, details] = xcorrtimes(C{ids(j)}.times,C{ids(k)}.times,'clip');
                result.xctime(n) = toc;
                result.ns(n,1) = length(C{ids(j)}.times);
                result.ns(n,2) = length(C{ids(k)}.times);
                result.synci(ids(j),ids(k),:) = SyncIndex(xc); 
                result.synci(ids(k),ids(j),:) = SyncIndex(xc);
                result.ccf{ids(j),ids(k)} = xc;
                result.cells(n,1) = j;
                result.cells(n,2) = k;
                result.probes(n,1) = j;
                result.probes(n,2) = k;
                result.exid(n) = 1;
                if nc <= 12 || matrixplot
                mysubplot(nc-1,nc-1,k+(j-2)*(nc-1));
                else
                mysubplot(np,np,n);
                end
                h = plot(xc);
                if GoodCluster(C{ids(j)}) > 1
                    set(h,'color','r');
                end
                axis('tight');
                yl = get(gca,'ylim');
                xl = get(gca,'xlim');
                set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0 yl(2)],'buttondownfcn',{@HitXcorrAll,0,[j k], 1});
                text(mean(xl),yl(2),sprintf('%d,%d,%.1f',ids(j),ids(k),result.synci(ids(j),ids(k),1)),'verticalalignment','top');
                drawnow;
            end
        end
    end
end
result.endtime = now;
result.xctimes = details.xpts;
setappdata(gcf,'xcorrs',result);

function [dips, scales, probes] = ClusterQualities(C)
k = 0;
for j = 1:length(C)
if isfield(C{j},'hdip')
    k = k+1;
    dips(1,k) = C{j}.hdip(1);
    if isfield(C{j},'dipsize')
        dips(2,k) = C{j}.dipsize(1);
    end
    dips(3,k) = C{j}.dpsum(1);
    dips(4,k) = abs(C{j}.dprime);
    dips(5,k) = C{j}.bmc;
    dips(6,k) = C{j}.mahal(1);
    if isfield(C{j},'bestd')
        dips(7,k) = max(C{j}.bestd);
    else
        dips(7,k) = 0;
    end
    if isfield(C{j},'bestd')
    else
        dips(8,k) = 0;
    end        
    probes(k) = j;
end
end
scales(1) = 1;
scales = range(dips(1,:))./range(dips,2);
            
            
function synci = SyncIndex(xc)
midpt = round(length(xc)+1)/2;
postmid = [midpt+2:midpt+20];
premid = [midpt-20:midpt-2];
allpts = [1:midpt-2  midpt+2:length(xc)];
refpts = [midpt-1 midpt+1];
synci(1) = xc(midpt)./mean(xc([premid postmid]));
synci(2) = xc(midpt)./mean(xc(allpts));
synci(3) = mean(xc(refpts))./mean(xc(allpts));
synci(4) = mean(xc([midpt-20:midpt+20]))./mean(xc(allpts));        


     
function PlotXcorrs(DATA, xcorrs, expts, bycell)
    ClearPlot;
    if bycell
        xcorrs = xcorrs([xcorrs.valid] == 1);
        cellids = cat(1,xcorrs.cells);
        cellids(isnan(cellids)) = 0;
    else
        cellids = cat(1,xcorrs.probes);
    end
    exids = cat(1,xcorrs.eid);
    weights = prod(cat(1,xcorrs.n)');
    cells = unique(cellids);
    cells = setdiff(cells, find(DATA.plot.xcorrexclude));
    cells = cells(cells > 0);
    probes = cat(1,xcorrs.probes);
    np = length(cells);
    ns = 0;
    for j = 1:length(cells)
        ida = find(cellids(:,1) == cells(j) & ismember(exids,expts));
        idb = find(cellids(:,2) == cells(j) & ismember(exids,expts));
        cellpos(j) = (sum(probes(ida,1))+sum(probes(idb,2)))./(length(ida)+length(idb));
        for k = 1:j
            ida = find(cellids(:,1) == cells(j) & cellids(:,2) == cells(k) & ismember(exids,expts));
            idb = find(cellids(:,2) == cells(j) & cellids(:,1) == cells(k) & ismember(exids,expts));
            if length(ida)+length(idb) > 0
                separation(j,k) =( sum([xcorrs(ida).separation]) - sum([xcorrs(idb).separation]))./(length(ida)+length(idb));
            else
                separation(j,k) = cellpos(j)-cellpos(k);
            end
            separation(k,j) = -separation(j,k);
        end
    end
    order = sum(separation < 0);
    [a,b] = sort(order);
    icells = cells; %unsorted
    cells = cells(b);
    for j = 1:length(cells)
        for k = 1:j
            id = find((cellids(:,1) == cells(j) & cellids(:,2) == cells(k)) | (cellids(:,2) == cells(j) & cellids(:,1) == cells(k)));
            id = id(ismember(exids(id),expts));
            if length(id)
                if length(id) > 1
                    xc = WeightedSum(cat(1,xcorrs(id).xc),weights(id));
                else
                    xc = xcorrs(id).xc;
                end
                mysubplot(np,np,k+(j-1)*np,'leftmargin',0.02);
                h = plot(-200:200,xc,'k-','linewidth',2);
                synci = SyncIndices(xc);
                if synci(2) < DATA.crit.synci
                    set(h,'color','r');
                end
                if ~isnan(synci(2))
                    ns = ns+1;
                    syncis(ns,1:length(synci)) = synci;
                    syncis(ns,3) = cells(j);
                    syncis(ns,4) = cells(k);
                end
                axis('tight');
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k) ns], expts});
                set(h,'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k)],expts});
                if k == j
                    [a,b] = find(DATA.CellList == cells(k));
                    p = 1+mod(b-1,DATA.nprobes);
                    if j == 1 || DATA.plot.xcorrlabeltop == 1
                        title(sprintf('Cell%d',cells(k)));
                    else
                        ii = find(icells == cells(k-1));
                        ij  = find(icells == cells(j));
                        if bycell
                            h = text(xl(1),yl(2),sprintf('C%d at %s (%.1f)',cells(k),ProbeLabel(p, DATA),separation(ii,ij)));
                            set(h,'HorizontalAlignment','left','verticalalignment','bottom');
                        else
                            title(sprintf('P%s',ProbeLabel(k, DATA)));
                        end
                    end
                end
                if k == 1
                    ylabel(sprintf('Cell%d',cells(j)));
                end
            end
        end
    end
    if DATA.plot.synci
    mysubplot(2,2,2);
    myscatter(syncis(:,1),syncis(:,2),'o','ids',syncis(:,3:4));
    end
    
function HitXcorrAll(a,b, type, cells, expts)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');
xcorrs = getappdata(get(a,'parent'),'xcorrs');
exids = cat(1,xcorrs.exid);
if type == 2
%    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    e = cells(1);
    xc = xcorrs(cells(3));
    mysubplot(2,2,2);
    xpts = DATA.xcorrval.times;
    plot(xpts,xc.xc,'k');
    if DATA.plot.xcmax < max(xpts) *1000
        set(gca,'xlim',[-DATA.plot.xcmax DATA.plot.xcmax]./1000);
    end
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    [a,b] = max(xc.xc);
    text(xl(2),yl(2),sprintf('Cell %d->%d P%d->%d %.0fms Exp%s',xc.cells(1),...
        xc.cells(2),cells(2),cells(3),xpts(b).*1000,DATA.expnames{e}),...
        'horizontalalignment','right','verticalalignment','top');
    line([0 0],yl,'linestyle','--');
    SetFigure(DATA,DATA.tag.spkmean);
    subplot(1,2,1);
    h = PlotMeanSpike(Clusters{e}{cells(2)},cells(2),cells(4),'imageonly');
    set(h,'ButtonDownFcn',{@HitXYPlot, e , cells(2)},'UserData',DATA.toplevel);
    subplot(1,2,2);
    h = PlotMeanSpike(Clusters{e}{cells(3)},cells(3),cells(5),'imageonly');
    set(h,'ButtonDownFcn',{@HitXYPlot, e , cells(3)},'UserData',DATA.toplevel);

    SetFigure(DATA,DATA.tag.spikes);
    if DATA.plot.synctmax > 0
        PlotSyncSpikes(DATA, e, cells([2 3]), cells([4 5]));
    end
    DATA.selectprobe(e,cells([2 3])) = 1;
    DATA.xcid = cells(6);
    
    if strcmp(DATA.plotexpttype,'means')
        Expts = getappdata(DATA.toplevel,'Expts');
        SetFigure(DATA,DATA.tag.expt);
        subplot(1,2,1);
        PlotExpt(Expts{e,cells(2)});
        subplot(1,2,2);
        PlotExpt(Expts{e,cells(3)});
    end
    set(DATA.toplevel,'UserData',DATA);
%    SpoolAllProbes(DATA, e, AllSpikes(e,cells([2 3])), Clusters{e});
    return;
end
cellids = cat(1,DATA.xcorrs.cells);
probes = cat(1,DATA.xcorrs.probes);
if type == 0 %all probes, not just cells
 ids = find(probes(:,1) == cells(1) & probes(:,2) == cells(2));
 pa = cells(1);
 pb = cells(2);
 np = length(unique(probes));
else
 cells([1 2]) = sort(cells([1 2]),'descend');
 ida = find(cellids(:,1) == cells(1) & cellids(:,2) == cells(2));
 idb = find(cellids(:,2) == cells(1) & cellids(:,1) == cells(2));
 ids = cat(1,ida, idb);
 ids = ids(ismember(exids(ids),expts));
 eid = ismember(exids,expts);
 icells = unique([DATA.xcorrs(eid).cells]);
 icells = icells(icells > 0);
 np = length(icells);
 pa = mean(cat(1,probes(ida,1), probes(idb,2))); 
 pb = mean(cat(1,probes(idb,1), probes(ida,2))); 
end

col = 1;
row = 1;
for j = 1:length(ids)
    p = ids(j);
    ca = xcorrs.cells(p,1);
    cb = xcorrs.cells(p,2);
    xc = xcorrs.ccf{ca,cb};
    eid = xcorrs.exid(p);
    col = col+1;
    if row > floor(np/2)
        if col > np
            if row > floor(np/2)
                row = row+1;
                col =row+1;
            end
        end            
    else
        if col > floor(np/2)
            row = row+1;
            if row == floor(np/2)
                row = row+1;
            end
            col = row+1;
        end
    end
    k = (row-1)*np + col;
    mysubplot(2,2,2,'leftmargin',0.02);
    h = plot(xcorrs.xctimes,xc,'r-');
    si = SyncIndices(xc);
    if isfield(DATA.plot,'synic') && DATA.plot.synci
        fprintf('E%dP%d,%d synci %s\n',e,ca,cb,sprintf('%.2f ',si));
    end
%    set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll, 2, [e X.probes X.clnum ids(j)],expts});
%    set(h,'buttondownfcn',{@HitXcorrAll, 2, [X.eid X.probes X.clnum ids(j)], expts});
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    text(xl(1),yl(2),sprintf('E%dP%d,%d',eid,ca,cb),'horizontalalignment','left','verticalalignment','top');
end

 mysubplot(2,2,2,'leftmargin',0.02);
 if length(ids) > 1
     weights = prod(cat(1,DATA.xcorrs(ida).n)');
     na = sum(weights);
     if length(ida) > 1
         xca = WeightedSum(cat(1,DATA.xcorrs(ida).xc),weights);
     else
         ida = DATA.xcorrs(ida).xc;
     end
     weights = prod(cat(1,DATA.xcorrs(idb).n)');
     nb = sum(weights);
     if length(idb) == 0
         xcb = [];
         nb = 0;
         xc = xca;
     else
         if length(idb) > 1
             xcb = WeightedSum(cat(1,DATA.xcorrs(idb).xc),weights);
         else
             xcb = DATA.xcorrs(idb).xc;
         end
         xc = WeightedSum(cat(1,xca,fliplr(xcb)),[na nb]);
     end
 else
     xc = xcorrs.ccf{ca,cb};
 end
 xpts = xcorrs.xctimes;
 plot(xpts,xc,'k-');
 line([0 0],get(gca,'ylim'),'linestyle','--','color','r');
yl = get(gca,'ylim');
xl = get(gca,'xlim');
[a,b] = max(xc);
text(xl(2),yl(2),sprintf('C%d->%d %.0fms (%.0f+-%.1f) P%.1f,%.1f',cells(1),cells(2),xpts(b).*1000,prctile(xc,50),std(xc),pa,pb),'horizontalalignment','right','verticalalignment','top');

function synci = SyncIndices(xc)
midpt = ceil(length(xc)./2);
synci(1) = mean([xc(midpt-1) xc(midpt+1)])./xc(midpt);
synci(2) = mean([xc(1:midpt-10) xc(midpt+10:end)])./xc(midpt);
