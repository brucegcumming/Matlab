function DATA = ReplotPCs(DATA,E, varargin)        setids = [];meanpos = [];autospace = 0;forcetofront = 0;args = {};j = 1;while j <= length(varargin)    if strncmpi(varargin{j},'autospace',6)        autospace = 1;    elseif strncmpi(varargin{j},'setid',5)        j = j+1;        setids = varargin{j};    elseif strncmpi(varargin{j},'showpos',7)        j = j+1;        meanpos = varargin{j};        j = j+1;        dimpos = varargin{j};    elseif strncmpi(varargin{j},'tofront',5)        forcetofront = 1;    end    j = j+1;endif DATA.interactive < 0    return;endif autospace    DATA.plottype = AllV.WhichPlotType(DATA.cluster, DATA.currentcluster);endif isempty(E)    cluster = DATA.cluster;else    cluster = E;endif forcetofront    figure(DATA.toplevel);else    AllV.SetFigure(DATA.tag.top, DATA);endclusterplot = [];if isfield(DATA,'dipvals')    dipvals = DATA.dipvals;else    dipvals = zeros(1,8);endif DATA.usebmi == 0dipvals = zeros(1,8);endif DATA.plottype == 2;    vpts = DATA.vpts;elseif DATA.plottype == 8;    vpts = DATA.dvpts;    AllVoltages=getappdata(DATA.toplevel,'AllVoltages');    DATA.dV= diff(AllVoltages,1,2);elseif DATA.plottype == 9    vpts = [1 1 2 2; 1 1 3 3];    vpts = repmat(vpts,4,1);end%need to get boundary anyway in case there are other clustersif isempty(E) && isfield(DATA,'cluster') && isfield(DATA.cluster,'space')    if DATA.cluster.space(1) == DATA.plottype %cluster is in this space        E = AllV.BoundaryFromCluster(E,DATA.cluster, DATA.currentcluster);    else        E = AllV.BoundaryFromCluster(E,DATA.cluster,DATA.currentcluster);    endendplots = DATA.pcplots(1:8,:);if DATA.plottype == 3 || DATA.plottype ==4    if DATA.usestdtemplates    plots = DATA.tmplots(17:24,:);    elseif DATA.plottype == 4        plots = DATA.tmplots(9:16,:);    else        plots = DATA.tmplots(1:8,:);    end    if isfield(DATA,'tmpdips')    dipvals = DATA.tmpdips;    endelseif DATA.plottype == 12 %Try new scores    tpt = find(DATA.spts == 0);    p = DATA.probe(1);    AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');    tvals = squeeze(AllVoltages(p,tpt,:));    pcs(:,1) = squeeze(AllVoltages(p,tpt+6,:))-tvals;    pcs(:,2) = squeeze(AllVoltages(p,tpt+19,:))-tvals;    pcs(:,3) = squeeze(AllVoltages(p,tpt-6,:))-tvals;    if p < size(AllVoltages,1)        pcs(:,4) = squeeze(AllVoltages(p+1,tpt,:));    elseif p > 1        pcs(:,4) = squeeze(AllVoltages(p-1,tpt+8,:));    end    if p > 1        pcs(:,5) = squeeze(AllVoltages(p-1,tpt,:));    elseif p< size(AllVoltages,1)        pcs(:,5) = squeeze(AllVoltages(p+1,tpt+8,:));    else        pcs(:,5) = squeeze(AllVoltages(p,tpt+8,:));    end    plots = DATA.pcplots;elseif DATA.plottype == 11 %show spaces needed for all cells    pcs = [];    for j = 1:length(DATA.xy)        pcs = cat(2,pcs,DATA.xy{j}(:,1));        pcs = cat(2,pcs,DATA.xy{j}(:,2));    end    if size(pcs,2) == 4        plots = [1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 1; 1 1];     elseif size(pcs,2) == 3        plots = [1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 1; 1 1];     elseif size(pcs,2) > 8        for j = 1:2:size(pcs,2)            plots(:,ceil(j/2)) = [j j+1];        end        k = 3;        for j = (ceil(j/2)+1):8            plots(:,j) = [1 k];            k = k+1;        end    else        plots = DATA.pcplots;    endelseif DATA.plottype == 7    plots = DATA.tmplots(17:24,:);elseif DATA.plottype == 5    plots = DATA.tmplots(9:16,:);elseif DATA.plottype == 1    Labels = AllV.PCLabels(DATA);    elseif DATA.plottype == 13    plots = [5 6; 1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 5;];end   if DATA.usegmcid & length(DATA.gmcid)    if isfield(DATA,'gmcid')        clid = DATA.gmcid;    elseif isfield(E,'bestcl')        clid = E.bestcl;    else    clid = DATA.clid;    DATA.usegmcid = 0; %don't use it if not defined    endelseif length(setids)    clid = setids;    args = {args{:} 'ptsz' [8 1]};else    args = {args{:} 'ptsz' DATA.ptsz};    if isfield(DATA,'clst')         clid = DATA.clst;    else        clid = DATA.clid;    endendp = DATA.probelist(DATA.probe(1));type = double(DATA.clplot);if DATA.clplot ==1    if DATA.plot.scaledensity        type = 2;    else        type(2) = DATA.plot.densityscale;    endendC = AllV.GetClusterDef(cluster, DATA.currentcluster);args = {args{:} 'clusterdef' 0};for j = 1:size(plots,1)    mysubplot(2,4,j);    clusterplot(:,j) = AllV.GetClusterPlots(DATA,cluster,plots,j);    if size(clusterplot,1)  < DATA.currentcluster%currentcluster no space defined yet        args{end} = 0;    elseif clusterplot(DATA.currentcluster,j) == j %this is the plot defining the cluster        args{end} = 1;    else        args{end} = 0;    end    t = 0;    if ismember(DATA.plottype, [3 4 7])        AllV.PlotPCs(DATA.TemplateScores,plots(j,1),plots(j,2),type,clid,DATA.colors,C,'fixrange',args{:});        tstr = sprintf('%s vs %s',DATA.TemplateLabels{plots(j,1)},DATA.TemplateLabels{plots(j,2)});        t =  title(tstr);        if isfield(E,'gmfit') && E.space(1) == 6 && E.space(2) ==4             [a,xi] = ismember(plots(j,1),DATA.tmplspace(1,:));            [b,yi] = ismember(plots(j,2),DATA.tmplspace(1,:));            if a && b && size(E.gmfit.mu,2) > max([xi yi])                hold on;                plot(E.gmfit.mu(1,xi),E.gmfit.mu(1,yi),'g+','markersize',10,'linewidth',2);                plot(E.gmfit.mu(2,xi),E.gmfit.mu(2,yi),'g+','markersize',10,'linewidth',2);            end                    end        set(gca,'UserData',plots(j,:));    elseif ismember(DATA.plottype,[6])        AllV.PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),type,clid,DATA.colors);        set(gca,'UserData',vpts(j,:));        if DATA.usegmcid && sum(ismember([vpts(1,3) vpts(3,4)],DATA.vspace) == 2)        end    elseif ismember(DATA.plottype,[2 6 8 9])        AllV.PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),type,clid,DATA.colors,'fixrange');        set(gca,'UserData',vpts(j,:));        if DATA.usegmcid && sum(ismember([vpts(j,2) vpts(j,4)],DATA.vspace)) == 2            k = find(DATA.vspace == vpts(j,2));            m = find(DATA.vspace == vpts(j,4));            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)        end    else        if DATA.usegmcid && sum(ismember(DATA.pcplots(j,:),[1:4])) == 2 && isfield(DATA.cluster,'gmfit')            showclustermeans = 1;        else            showclustermeans = 0;        end        if ismember(DATA.plottype,[11 12]) %tests            AllV.PlotPCs(pcs,plots(j,1),plots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});        elseif ismember(DATA.plottype,[13]) %james autocut            AllV.PlotFeatures(DATA,plots(j,1),plots(j,2),type,clid,DATA.colors,C, 'fixrange');            showclustermeans = 0;        elseif ismember(DATA.plottype,[10]) %ICA            AllV.PlotPCs(DATA.icas,DATA.pcplots(j,1),DATA.pcplots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});        else            AllV.PlotPCs(DATA.pcs,DATA.pcplots(j,1),DATA.pcplots(j,2),type,clid,DATA.colors, C,'fixrange',args{:});            t = title(Labels{j});        end        set(gca,'UserData',plots(j,:));        if showclustermeans            k = find(DATA.pcspace == DATA.pcplots(j,1));            m = find(DATA.pcspace == DATA.pcplots(j,2));            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)        end        if j == 1            text(0,1.2,sprintf('E%dP%d',DATA.exptno,p),'units','normalized','fontsize',DATA.gui.fontsize(1));        end    end    if ~ismember(j,[1 4])        set(gca,'ytick',[]);    end    if ~ismember(j,[5:8])        set(gca,'xtick',[]);    end    if isempty(setids)    axis('tight');              end    xl = get(gca,'Xlim');    yl = get(gca,'Ylim');    if t && ishandle(t)        set(t,'VerticalAlignment','top','position',[mean(xl) yl(2)]);    end        if DATA.showdipvals    text(mean(xl),yl(2),sprintf('%.2f',dipvals(j)),'VerticalAlignment','top','color','k');    end    %would like to used cluster, not E in future, so it doesn't depend on    %which is current cluster%    clusterplot(:,j) = AllV.GetClusterPlots(DATA,DATA.cluster,plots,j);endif isfield(DATA.Expt,'Header') && isfield(DATA.Expt.Header,'expname')       exname = DATA.Expt.Header.expname;   else       exname = [];   endif DATA.plottype == 1    p = DATA.probelist(DATA.probe);    if isfield(DATA,'alldips')        if size(DATA.alldips,2) > 1            b = max(DATA.alldips');        else            b = DATA.alldips;        end        tstr = sprintf('%s %.2f, dt%.2f, csd%.2f',num2str(p),b(1).*100,b(2).*100,b(3).*100);    elseif DATA.dvdt        tstr = sprintf('E%dP%d: dvdt %s',DATA.exptno,p(1),exname);    elseif DATA.csd        tstr = sprintf('E%dP%d: csd %s',DATA.exptno,p(1),exname);    else        tstr = sprintf('E%dP%d %s',DATA.exptno,p(1),exname);    endelse   tstr = sprintf('E%dP%d%s %s',DATA.exptno,p(1),DATA.probelabel,exname);endmysubplot(2,4,1);th = text(0.5,0.0,tstr,'units','normalized','VerticalAlignment','Top','fontsize',DATA.gui.fontsize(1));if isfield(DATA,'cluster') && isfield(DATA.cluster,'auto') && DATA.cluster.auto == 0    set(th,'color','r');endset(th,'Tag','PCTitleString');DATA.maintitle = th;DATA.clustericon = AllV.SetClusterIcon(DATA);    Cs = cluster;for j = 1:size(clusterplot,2)for k = 1:size(clusterplot,1)if clusterplot(k,j)    mysubplot(2,4,clusterplot(k,j)); %need to find right graph    hold on;    if k > 1        Cs.next{k-1}.color = DATA.colors{k+1};        DATA.elmousept.h = AllV.DrawEllipse(Cs.next{k-1});    else    Cs.color = DATA.colors{2};    DATA.elmousept.h = AllV.DrawEllipse(Cs);    end    DATA.elmousept.handles(k) = DATA.elmousept.h;    if Cs.shape == 1        tmp = Cs;        tmp.pos(1) = mean(E.pos([1 3])) + diff(E.pos([2 4]))/2;        tmp.pos(3) = mean(E.pos([1 3])) - diff(E.pos([2 4]))/2;        tmp.pos(2) = mean(E.pos([2 4])) - diff(E.pos([1 3]))/2;        tmp.pos(4) = mean(E.pos([2 4])) + diff(E.pos([1 3]))/2;        h = AllV.DrawEllipse(tmp,'r');        set(h,'linestyle',':');    end    hold off; endendendif length(DATA.elmousept.handles) >= DATA.currentcluster;    DATA.elmousept.h = DATA.elmousept.handles(DATA.currentcluster);end[iscell, cellid] =  AllV.isacell(DATA, DATA.exptno, AllV.ProbeNumber(DATA));if iscell    mysubplot(2,4,4);    for j = 1:length(cellid)        if cellid(j) > 0        text(1,0.1*j,sprintf('Cell %d',cellid(j)),'units','normalized','color',DATA.colors{j+1},...            'HorizontalAlignment','Right','fontsize',DATA.gui.fontsize(1));        end    endendif isfield(DATA,'energy')  && DATA.plot.vare    AllV.SetFigure(DATA.tag.vare, DATA);    subplot(1,1,1);    AllV.PlotVarE(DATA);end