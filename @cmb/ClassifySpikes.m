function ClassifySpikes(mousept,src,varargin)

DATA = GetDataFromFig(src);
cl = mousept.cluster;
p = get(gca,'UserData');
if isempty(p)
    p = DATA.probe;
end

DATA.currentcluster = cl;
C.x = [mousept.c(1) mousept.r(1) mousept.xrange];
if C.x(2) == 0 & isfield(DATA.cluster{cl,p},'x')
    C.x(2) = DATA.cluster{cl,p}.x(2);
end
C.y = [mousept.c(2) mousept.r(2) mousept.yrange];
if C.y(2) == 0 & isfield(DATA.cluster{cl,p},'y')
    C.y(2) = DATA.cluster{cl,p}.y(2);
end
C.angle = -mousept.angle;
C.h = mousept.lasth;
C.params = [DATA.plot.clusterX DATA.plot.clusterY];
C.Arange = DATA.clusterArange;
C.Brange = DATA.clusterBrange;
C.Erange = DATA.clusterErange;
C.deleted = 0;
if DATA.forceclusterid > 0 & DATA.forceclusterid ~= cl
    C.forceid = DATA.forceclusterid;
end
if isfield(DATA,'AllClusters')
    if iscell(DATA.AllClusters)
        expspks = 1:length(DATA.AllClusters{DATA.currentexpt(1)}(p).times);
    else
        expspks = DATA.AllClusters(p).spklist;
    end
elseif isfield(DATA,'AllSpikes') & isfield(DATA.AllSpikes{p},'spklist')
    if DATA.plot.synccluster > 0
        pa = DATA.syncprobes(2);
        C.firstspk = DATA.AllSpikes{pa}.spklist(1);
        C.lastspk = DATA.AllSpikes{pa}.spklist(end);
        DATA.cluster{cl,pa}.Cluster = C;
        [DATA,dprime,nc] = SetSpkCodes(DATA,DATA.AllSpikes{pa}.spklist,pa,1);
        expspks = DATA.AllSpikes{p}.spklist;
        expspks = DATA.syncspikes(:,1);
    else
        expspks = DATA.AllSpikes{p}.spklist;
    end
elseif isfield(DATA,'spklist') && ~isempty(DATA.spklist)
    expspks = DATA.spklist;
else
    expspks = DATA.spkrange(1):DATA.spkrange(2);
end
C.firstspk = expspks(1);
C.lastspk = expspks(end);
if DATA.firsttrial > 1
    DATA.cluster{cl,p}.Cluster = C;
    DATA.cluster{cl,p}.lastspk = C.firstspk-1;
else
    DATA.cluster{cl,p} = C;
end
DATA.cluster{cl,p}.touched = 1;
DATA.newcluster(mousept.cluster) = 1;
%mousept.lasth
colors = mycolors;

%mode 5 is a touch inside the ellipse. Want to classiy if it has moved
%though (moe goes to 11)
if mousept.mode ~= 5   || strcmp(get(src,'tag'),'AllProbeSpikes')
    [DATA, dprime, nc] = SetSpkCodes(DATA,expspks,p,2);
    if length(nc) > 1
        fprintf('%s: %d Clusters for Probe %d\n',DATA.explabels{DATA.currentexpt(1)},length(nc),p);
    end
end
DATA.state.recut = 1;
cmb.SetGui(DATA);
if isfield(DATA,'AllClusters')
    if iscell(DATA.AllClusters)
        id = find(DATA.AllClusters{DATA.currentexpt(1)}(p).codes(expspks) == 0);
    else
        id = find(DATA.AllClusters(p).codes(expspks) == 0);
    end
elseif isfield(DATA,'AllSpikes')
    id = find(DATA.AllSpikes{p}.codes(expspks,2) == 0);
else
    id = find(DATA.AllData.Spikes.codes(expspks,2) == 0);
end
if ~DATA.densityplot
    % done by SetSpkCodes
    %    cmb.PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{1});
end

%set Expts.cluster before countspikes, because that can copy Expts.cluster
%to DATA.cluster
p = get(gca,'UserData');
if isempty(p)
    DATA.Expts{DATA.currentexpt(1)}.Cluster = DATA.cluster;
end

if DATA.state.autoreplotgraph
    DATA = cmb.CountSpikes(DATA, DATA.currentexpt(1),'replot');
else
    DATA = cmb.CountSpikes(DATA, DATA.currentexpt(1));
end
DATA.Expts{DATA.currentexpt(1)}.gui.classified = 2;
if DATA.plot.showdprime
    GetFigure('DprimeCalc');
    cmb.PlotDprime(DATA);
end
cmb.plotISI(DATA);
set(0,'currentfigure',src);
set(DATA.toplevel,'UserData',DATA);
if DATA.state.verbose
    DATA.cluster{DATA.currentcluster,DATA.probe}
end
if DATA.plot.quickspks
    cmb.PlaySpikes(DATA,DATA.currentexpt(1),'quickspks');
end


