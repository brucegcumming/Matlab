function DATA = DrawXYPlot(DATA, expspks, varargin)

ho = ishold;
showhist = 0;
if ~isfield(DATA,'nclusters')
    DATA.nclusters = 0;
end
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
DATA.ptsize = cmb.CheckPtSize(DATA,length(expspks));
delete(get(gca,'children'));
expid = DATA.currentexpt(1);
if isfield(DATA,'AllSpikes')
    expspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
    if ~isfield(DATA,'nclusters')
        DATA.nclusters = 0;
    end
    if length(expspks)
        cmb.PlotSpikeXYClusters(DATA, expspks, 1:DATA.nclusters+1,varargin{:});
    end
else
    if length(expspks) < 10
        expspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
    end
    if isempty(expspks)
       return;
    end
    if DATA.densityplot
        cmb.PlotXYDensity(DATA.Spikes.cx(expspks),DATA.Spikes.cy(expspks));
    else
        nc = cmb.MaxSpikeCode(DATA, expspks);
        cmb.PlotSpikeXYClusters(DATA,expspks,1:nc+1,varargin{:});
    end
end
hold on;
DATA = cmb.DrawClusters(DATA, DATA.cluster, 0);
if isfield(DATA.Expts{expid},'OnlineCluster') & DATA.state.showonlineclusters
    DATA = cmb.DrawClusters(DATA,DATA.Expts{expid}.OnlineCluster,0);
end

cmb.CalcClusterdp(DATA,1);
if showhist
    xl = get(gca,'xlim');
    [a,b] = smhist(DATA.Spikes.cx(expspks));
    plot(b,(a-xl(1)).* diff(xl)./max(a));
end
if ismember(DATA.plot.autoscale,[2 3]) && length(DATA.Expts{DATA.currentexpt(1)}.gui.spks)
    expspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
    DATA.spklist = expspks;
    [xr, yr] = cmb.ClusterRange(DATA.cluster,DATA.probe);
    cx = DATA.Spikes.cx(expspks);
    cy = DATA.Spikes.cy(expspks);
    if DATA.plot.autoscale == 2
        y(2) = prctile(cy,99.5).*2;
        x(2) = prctile(cx,99.5).*2;
    else
        y(2) = prctile(cy,95).*2;
        x(2) = prctile(cx,95).*2;
    end
    y(2) = min([y(2) max(cy)]);
    y(2) = max([y(2) yr(2)]); %make sure cluster ellipse is visible
    x(2) = min([x(2) max(cx)]);
    x(2) = max([x(2) xr(2)]);
    if min(cy) > 0 & min(cy) < y(2)/10
        y(1) = 0;
    else
        y(1) = min(cy);
    end
    if min(cx) > 0 & min(cx) < x(2)/10
        x(1) = 0;
    else
        x(1) = min(cx);
    end
    if x(2) <= x(1)
        x(2) = x(1) + abs(x(1));
    end
    if y(2) <= y(1)
        y(2) = y(1) + abs(y(1));
    end
    DATA.plot.clusterYrange = y;
    DATA.plot.clusterXrange = x;
end
DATA = cmb.FinishXYPlot(DATA);
if ~ho
    hold off;
end

