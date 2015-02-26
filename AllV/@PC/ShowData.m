  function DATA  = ShowData(DATA, ex,p, varargin)    showall = 1;  j = 1;  while j <= length(varargin)      if strncmpi(varargin{j},'oneprobe',7)          showall = 0;      end      j = j+1;  endif nargin == 1    ex = DATA.currentpoint(1);    p = DATA.currentpoint(2);end[Clusters, DATA] = PC.CheckClusterLoaded(DATA, ex);Expts = getappdata(DATA.toplevel,'Expts');    oldf = gcf;Expt = [];eid = DATA.exptid(ex);if length(Expts) >= eid    DATA.Expt  = Expts{eid};endif showall == 0    it = findobj(DATA.toplevel,'Tag','ProbeList');    set(it,'value',p)    it = findobj(DATA.toplevel,'Tag','ExptList');    set(it,'value',ex)endif DATA.datatype == 2    C = Clusters{ex}{p}.cluster{DATA.templatesrc};elseif DATA.useautoclusters    AutoClusters = getappdata(DATA.toplevel,'AutoClusters');    C = AutoClusters{ex}{p};elseC = Clusters{ex}{p};endautoC = [];if DATA.plot.autocompare && isappdata(DATA.toplevel,'AutoClusters')    AutoClusters = getappdata(DATA.toplevel,'AutoClusters');    if length(AutoClusters) >= ex && length(AutoClusters{1}) >= p        autoC = AutoClusters{ex}{p};    endend    if ~C.probe(1)    C.probe(1)=p;endif C.exptid <= 0    fprintf('!!!Ex %d P %d has 0 exptid\n',ex,p);    C.exptid = ex;endif isfield(C,'Evec')    eveci = C.Evec.Eval(1)./sum(C.Evec.Eval);else    eveci = NaN;endif isfield(C.MeanSpike,'muxc')    muxc = C.MeanSpike.muxc;else    muxc = 0;endif isfield(C,'bestspace')    bestspace = C.bestspace(2);else    bestspace = NaN;    C.bestspace = 0;endDATA.cellplotmenu = PC.AddContextMenu(DATA,'cellplot');DATA = PC.MarkCurrentCluster(DATA);PC.PrintComments(DATA,ex,p);if isfield(C,'User')    user = C.user;else    user = '??';endfprintf('Mahal ND %.2f(%d), 2D %.2f, 1D %.2f. Dropi %.2f (T%.2f). Made %s by %s\n',C.bestspace(1),bestspace,C.mahal(1),C.mahal(4),C.dropi(3),C.Trigger,datestr(C.ctime),user)spkr = C.ncut./C.nspks;if isfield(DATA,'GaussFitdp') && size(DATA.GaussFitdp,1) >= ex    fprintf('Fit %.1f %.1f spkratio %.3f\n',DATA.GaussFitdp(ex,p,1),DATA.GaussFitdp(ex,p,2),spkr);endstr = [];if strmatch(DATA.plot.alltype,{'followcorr' 'templatescore'})    if DATA.datatype == 1        xc = PC.CalcTemplateXcorr(Clusters{DATA.templateid(2)}{DATA.templateid(1)},Clusters{ex}{p});    else        xc = PC.CalcTemplateXcorr(DATA.Templates{DATA.templatesrc},Clusters{ex}{p}.cluster{DATA.templatesrc});    end        str = sprintf(' xc %.2f',xc);elseif strmatch(DATA.plot.alltype,{'BuildTimes'})        str = sprintf(' took %.2f',(DATA.ctimes(ex,p,3)-DATA.ctimes(ex,p,1))*24);elseif strmatch(DATA.plot.alltype,{'PcGms'})    [a,b] = max(C.pcgms);        str = ['GM' sprintf(' %.2f',C.pcgms)];elseif strmatch(DATA.plot.alltype,{'Tagged'}) & DATA.tagged(ex,p)    DATA.tagstrings = {'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear'};    fprintf('Tagged %s\n',DATA.tagstrings{DATA.tagged(ex,p)});endif isfield(C,'first') & isfield(C.first,'needmore')    str = [ str sprintf('Nmore %d',C.first.needmore)];endif isfield(C,'duration')    rates(1) = C.ncut./C.duration;    rates(2) = C.nspks./C.duration;elserates(1) = 0;rates(2) = 0;endfprintf('Evi %.2f, muxc %.2f %s. Rates %.1f (su),%.1f(mu)\n',eveci,muxc,str,rates);exptno = DATA.exptlist(DATA.currentpoint(1));str = sprintf('Ex %.1f: %s',exptno,DATA.expnames{exptno});fprintf('%s\n',str);PC.SetFigure(DATA, DATA.tag.all);title(str);if DATA.plot.allxy && showall    PC.PlotAllProbeXY(DATA);endif DATA.plot.xyseq    PC.PlotXYSequence(DATA, [ex p]);endif DATA.plot.isihist    PC.PlotISIHist(DATA, [ex p]);endif DATA.plot.quickspks    PC.SetFigure(DATA,DATA.tag.spikes);    h = PC.QuickSpikes(DATA, [ex p]);    PC.AddLineContextMenu(DATA, h, ex, p);    drawnow;endif DATA.show.allvpcs && showall == 0    PC.CallAllVPcs(DATA,ex,p);endDATA = PC.ConditionalPlotXY(DATA, C, 0);if DATA.plot.autocompare && ~isempty(autoC)    DATA = PC.ConditionalPlotXY(DATA, autoC, 0);endDATA.Expt = PC.PlotExptCounts(DATA, ex,p, DATA.currentcluster);if DATA.plot.hist || DATA.plot.refitgm    PC.SetFigure(DATA,DATA.tag.hist,'front');    x = PC.PlotClusterHistogram(DATA, C, DATA.plot.refitgm,'cluster', DATA.currentcluster);endif DATA.plot.refitgm    if C.shape == 0        C.r = PC.CalcRadius(C, C.xy);        [a,b] = GMDip(C.r,0,'idlist',C.clst);    else        [a,b] = GMDip(C.xy,0,'idlist',C.clst);    end    scale = max(x.nsp)./max(b.gxy(:,2));    if isfield(DATA,'GaussFitdp')        DATA.GaussFitdp(ex,p,1) = b.gmdprime;        [dp, d, c] = PC.Fit2Gauss(C);        DATA.GaussFitdp(ex,p,2) = dp;        DATA.GaussFitdp(ex,p,3) = now;        DATA.gmfitpos(ex,p,:)  = c.fitpos;    end    plot(b.gxy(:,1),b.gxy(:,2).*scale,'g');    plot(b.gxy(:,1),b.gxy(:,3).*scale,'g');    plot(b.gxy(:,1),sum(b.gxy(:,2:3),2).*scale,'r');    title(sprintf('%d/%d events M%.2f->%.2f Fit %.2f%s',C.ncut,size(C.xy,1),C.mahal(4),b.gmdprime,dp,x.qstr));    [a,b,c] = GMfit(C.xy,2+length(C.next),1);    id = cluster(a,C.xy);    sgn = (mean(C.xy(id==1))-mean(C.xy(id==2))) .* C.sign;    if  sgn > 0        mu =2; su = 1;    else        mu=1; su=2;    end    PC.SetFigure(DATA,DATA.tag.xyplot,'front');    hold off;    if sum(id==mu) < 200        sz =  5;    else        sz = 1;    end    plot(C.xy(id==mu,1),C.xy(id==mu,2),'.','markersize',sz);    hold on;    if sum(id==su) < 200        sz =  5;    else        sz = 1;    end    plot(C.xy(id==su,1),C.xy(id==su,2),'r.','markersize',sz);    for j = 1:length(C.next)        plot(C.xy(id==2+j,1),C.xy(id==2+j,2),'.','markersize',sz,'color','g');    end    title(sprintf('P%d Ex %.0f Gm %.2f(2D%.2f) (%.2f,%.2f for space %.0f)%.0f',p,C.exptno,C.mahal(4),C.mahal(1),C.mahal(2),C.bestspace(1),C.bestspace(2),C.sign));endif DATA.plot.spkmean    PC.SetFigure(DATA, DATA.tag.spkmean,'front');    hold off;    PC.PlotMeanSpike(C,p,0,'addtitle',str,DATA.plotmeantype, DATA);endif DATA.spoolspikes && showall == 1    xs = '';    if rem(C.exptno,1) > 0.001        xs = 'a';    end    a = PC.SpoolSpikeFile(DATA,DATA.currentpoint(1),DATA.currentpoint(2));endif DATA.plot.autocompare%    F = GetFigure('CompareAuto');%    PC.CompareClusters(C,autoC, DATA);end    if DATA.steptype == 2    if ishandle(DATA.markcc)        delete(DATA.markcc);    end    PC.SetFigure(DATA, DATA.tag.clusters);    DATA.markcc = PC.DrawBox(ex, p,1);endDATA.currentpoint = [ex, p];[a,cells,clid] = PC.isacell(DATA,ex,p);for c = 1:length(clid)    if clid(c) > 1 && (length(C.next) < clid(c)-1 || isempty(C.next{clid(c)-1}))        s = sprintf('Expt%d(Row%d)P%d Cluster %d is Cell %d but its empty\n',DATA.exptid(ex),ex,p,clid(c),cells(c));        acknowledge(s,DATA.toplevel,'print');    endendset(DATA.toplevel,'UserData',DATA);DATA = PC.PlotCellList(DATA,'showfig');drawnow;figure(oldf);