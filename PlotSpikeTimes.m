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
    subplot(1,1,1);
    minbmc = 0.2;
    k = 0;
    for j = 1:length(C) 
        if isfield(C{j},'times') && GoodCluster(C{j})
            k = k+1;
            ids(k) = j;
        end
    end
    result.cells = ids;
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
                xc = xcorrtimes(C{ids(j)}.times,C{ids(k)}.times,'clip');
                result.xctime(n) = toc;
                result.ns(n,1) = length(C{ids(j)}.times);
                result.ns(n,2) = length(C{ids(k)}.times);
                result.synci(ids(j),ids(k),:) = SyncIndex(xc); 
                result.synci(ids(k),ids(j),:) = SyncIndex(xc);
                result.ccf{ids(j),ids(k)} = xc;
                if nc <= 12
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
                set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0 yl(2)]);
                text(mean(xl),yl(2),sprintf('%d,%d,%.1f',ids(j),ids(k),result.synci(ids(j),ids(k),1)),'verticalalignment','top');
                drawnow;
            end
        end
    end
end
result.endtime = now;


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