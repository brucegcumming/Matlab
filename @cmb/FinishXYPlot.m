function DATA = FinishXYPlot(DATA)
if ismember(DATA.plot.autoscale,[0 2 3 4]) % 2 means I autoscale
    DATA = cmb.SetXYRanges(DATA);
    if diff(DATA.plot.clusterXrange) > 0
        set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange);
    end
    
end




xlabel(DATA.spkvarnames{DATA.plot.clusterX});
ylabel(DATA.spkvarnames{DATA.plot.clusterY});
it = findobj(DATA.xyfig, 'Tag','Clusterid');
if isempty(it)
    c = 1;
else
    c = get(it,'value');
end
ei = DATA.currentexpt(1);
if ~isempty(DATA.explabels{DATA.currentexpt(1)})
    expname = DATA.explabels{DATA.currentexpt(1)};
    expname = DATA.Expts{DATA.currentexpt(1)}.Header.expname;
else
    expname = DATA.Expts{DATA.currentexpt(1)}.Header.expname;
end
if DATA.firsttrial > 1
    expname = [expname sprintf('from %d',DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.firsttrial).Trial)];
end
p = DATA.probe;
str = '';
if isfield(DATA,'AllClusters')
    ctype = 1;
else
    ctype = 2;
end

if isfield(DATA,'cluster') & cmb.iscluster(DATA.cluster,c,p)
    if isfield(DATA.cluster{c,p},'dprime')
        str = sprintf('d=%.2f',DATA.cluster{c,p}.dprime);
    end
    if isfield(DATA,'AllSpikes')
        id = find(DATA.AllSpikes{DATA.probe}.codes(DATA.spklist,2) == c);
    else
        id = find(DATA.AllData.Spikes.codes(DATA.spklist,ctype) == c);
    end
    if isfield(DATA.Expts{ei}.Stimvals,'du')
        rate = length(id)./(length(DATA.Expts{ei}.Trials) .* DATA.Expts{ei}.Stimvals.du);
    else
        rate = length(id)./length(DATA.Expts{ei}.Trials);
    end
    
    title(sprintf('%s: C%d %sX%s %s P%d %.1fHz(%d/%d)',expname,c,...
        DATA.spkvarnames{DATA.cluster{c,p}.params(1)},...
        DATA.spkvarnames{DATA.cluster{c,p}.params(2)},str,DATA.probe,...
        rate,length(id),length(DATA.spklist)));
elseif isfield(DATA,'AllClusters') && iscell(DATA.AllClusters)
    probe = GetProbe(DATA, DATA.currentexpt(1), DATA.probe);
    cluster = DATA.AllClusters{DATA.currentexpt(1)}(probe);
    if ~isfield(cluster,'suffix')
        cluster.suffix = NaN;
    end
    title(sprintf('%s: Cell%d (P%d) Sufffix%d',expname,DATA.probe,probe,cluster.suffix));
else
    title(sprintf('%s: C%d',expname,c));
end
if DATA.plot.clusterZ > 0
    cmb.Plot3DClusters(DATA,0);
end
