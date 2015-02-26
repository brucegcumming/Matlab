function [cl, cluster, xy] = ClassifySpikes(DATA, E, varargin)    quick = 0;    forcesign = 0;    replot = 1;    recluster = 0;    j = 1;%    if quickmode.quick is nonzero, then avoid slow steps that are not%    always neeeded.  other flags in quickmode detemine what fits are used.    quickmode.quickest = 0;    quickmode.quick = 0;    quickmode.fit1cut = 0;    quickmode.fit2gauss = 0;    quickmode.dropi = 0;    quickmode.mean = 0;    quickmode.dips = 0;    quickmode.triggerhist = 1;    quickmode.fitallcut = 0;        while j <= length(varargin)        if isfield(varargin{j}, 'quickest')            quickmode = varargin{j};            quickmode.quick = 1;        elseif strncmpi(varargin{j},'noplot',4)            replot = 0;        elseif strncmpi(varargin{j},'notrigger',6)            quickmode.triggerhist = 0;        elseif strncmpi(varargin{j},'quick',4)            quickmode.quick = 1;        elseif strncmpi(varargin{j},'recluster',4) %Not chnaging            recluster = 1;            if isfield(E,'sign')                forcesign = E.sign;            end        elseif strncmpi(varargin{j},'sign',4)            j = j+1;            forcesign = varargin{j};        end        j = j+1;    end        if isfield(E,'automode') && strcmp(E.automode,'james')%        DATA = AllV.JamesAutoCut(DATA, 'reapply', E);%        return;    end    cl = [];    cluster = [];    xy = [];    if ~isfield(E,'next') && isfield(DATA.cluster,'next')        E.next = DATA.cluster.next;    end        if ~isfield(E,'pcplot') && isfield(E,'space')        cluster = E;        if cluster.space(1) == 6            if cluster.space(2) == 3                DATA.plottype = 2;            elseif cluster.space(2) == 4                DATA.plottype = 3;            else                DATA.plottype = cluster.space(2);            end        elseif cluster.space(1) == 5  %Var/E plot is another plot            DATA.plottype = 1;        else            DATA.plottype = cluster.space(1);        end        if DATA.currentcluster > 1 & isempty(cluster.next{DATA.currentcluster-1})            return;        end        E = AllV.BoundaryFromCluster([],cluster, DATA.currentcluster);    else        E.boundarytype = 1;    end    p = E.pcplot;    cx = E.xyr(1);    cy = E.xyr(2);    if length(E.xyr) > 3    rx = E.xyr(3);    ry = E.xyr(4);    end    if forcesign        cluster.sign = forcesign;        E.sign = forcesign;    end            %need to maek this exclsion not incluaion    f = {'bestspace' 'bestd' 'auto' 'autotook' 'firstspace' 'firstbmi' 'bestll' 'gmdip' 'gmdipres' 'gmfit1d' 'gmfit' 'gmfits' 'aspectratio'...        'next' 'pos' 'TemplateUsed' 'mumeanUsed' 'MeanSpike' 'eventrate' 'exptreadmethod'...        'mydip' 'mydipsize' 'manual' 'automode' 'trigparams' 'jamescluster' 'neednewtemplate'};    f = setdiff(fields(E),{'newscores' 'pcgms' 'ctime' 'space' 'pcplot' 'h'});    cluster = CopyFields(cluster,E,f);    ispk = DATA.probe(1);    cluster.nspks = DATA.nevents;    cluster.minenergy = DATA.minenergy;    cluster.minvar = DATA.minvar;    cluster.Trigger = DATA.Trigger;    cluster.spts = DATA.spts;    cluster.dvdt = DATA.dvdt;    cluster.csd = DATA.csd;    cluster.ctime = now;    cluster.eveci = DATA.Evec.Eval(1)./sum(DATA.Evec.Eval);    cluster.pcgms = DATA.dipvals;    cluster.probe = AllV.ProbeNumber(DATA);    cluster.ctype = 0;    cluster.neednewtemplate = 0;         if length(DATA.excludetrialids)        cluster.excludetrialids = DATA.excludetrialids;    end    clnum = DATA.currentcluster+1;        if isfield(DATA.Expt,'exptno')        cluster.exptno = DATA.Expt.exptno;    end    if ~isempty(DATA.lastcut)        cluster.first = DATA.lastcut;    end    if isfield(DATA,'clst')        cl.clst = DATA.clst;    else        cl.clst = ones(DATA.nevents,1);    end    if isempty(cl.clst)        cl.clst = ones(DATA.nevents,1);    end            if isfield(E,'cluster')        cluster.cluster = E.cluster;    elseif DATA.hidecluster        cluster.cluster = DATA.hidecluster+1;    else        cluster.cluster = 1;    end    cluster.cluster = clnum-1;    if clnum == 2        C = cluster;    elseif length(cluster.next) > clnum-3        C = cluster.next{clnum-2};    else %new cut         C = [];    end        if DATA.verbose > 1        fprintf('Classify: %s Space %s\n',AllV.IDStr(DATA),sprintf('%d ',E.space));    end    if isfield(E,'usegmcluster') && E.usegmcluster == 1        xy = DATA.ndxy;        id = find(E.bestcl == 1);        nid = find(E.bestcl == 2);        e(1) = mean(DATA.energy(1,id));        e(2) = mean(DATA.energy(1,nid));        if diff(e) > 0            id = find(E.bestcl == 2);            nid = find(E.bestcl == 1);            cluster.sign = -1;        else            cluster.sign = 1;        end        r = xy(:,1);        cluster.space = E.space;        cluster.crit = (mean(xy(id,1))+mean(xy(nid,1)))/2;        x = AllV.FitGaussMeans(xy(DATA.uid,:),2,'clusterid',id);        cluster.shape = 3;    elseif E.shape == 2 || (E.space(1) == 6 && E.shape == 1) %line in n-D space        cluster.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),cluster.angle);        cluster.y = exy(:,2);        cluster.crit = mean(exy(:,1));        cluster.sign = E.sign;        cluster.space = E.space;        xy = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),cluster.angle);        r = xy(:,1);        if cluster.sign < 0            id = find(xy(DATA.uid,1) < cluster.crit);            nid = find(xy(DATA.uid,1) >= cluster.crit);        else            cluster.sign = 1;            id = find(xy(DATA.uid,1) > cluster.crit);            nid = find(xy(DATA.uid,1) <= cluster.crit);        end        aid = id;        id = DATA.uid(id);        nid = DATA.uid(nid);        e(1) = mean(DATA.energy(1,id));        e(2) = mean(DATA.energy(1,nid));        if e(2) > e(1) && forcesign == 0            cluster.sign = -cluster.sign;            cid = nid;            nid = id;            id = cid;        end        if quickmode.dips            cluster.bmc = AllV.BimodalCoeff(r,1.5);            cluster.dprime = AllV.CalcDprime(r(id),r(nid));            cluster.hdip = HartigansDipTest(sort(r));        end        x = AllV.FitGaussMeans(xy(DATA.uid,:),2,'clusterid',aid);%        [a,b] = FindDip(r,DATA.energy(1,:),'eval',cluster.crit);        if ~isfield(E,'gmdprime')            [a,b] = GMDip(xy(DATA.uid,1),DATA.energy(1,DATA.uid),'eval',cluster.crit,'label',DATA.idstr);            cluster.gmdprime = b.gmdprime;            cluster.mahal(4) = b.gmdprime;            cluster.autodipsize = b.dipsize;            cluster.dipsize = b.cdipsize;            cluster.gmfit1d = b.G{b.best};        end    elseif E.shape == 1 %line% for line cuts, rotate space so that x > crit defines Cluster% ClsuterDetails.xy will then have rotated values. Need rotating back in% Plotclsuters.        angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);        cluster.y = exy(:,2);        crit = mean(exy(:,1));        cluster.angle = angle;        cluster.crit = crit;        if isempty(p) %Var/E plot            y = DATA.spkvar(ispk,DATA.uid)./DATA.energy(1,DATA.uid);            xy = xyrotate(DATA.energy(1,:),y,angle);           cluster.space = [5];        elseif DATA.plottype == 2 %voltage pairs            AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');            xy = xyrotate(AllVoltages(p(1),p(2),:),AllVoltages(p(3),p(4),:),angle);           cluster.space = [2 p(1:4)];        elseif ismember(DATA.plottype,[3 4 7])            xy = xyrotate(DATA.TemplateScores(:,p(1)),DATA.TemplateScores(:,p(2)),angle);            cluster.space = [3 p(1:2)];        elseif E.space(1) == 10            xy = xyrotate(DATA.icas(:,p(1)),DATA.icas(:,p(2)),angle);            cluster.space = [E.space(1) p(1:2)];        elseif ismember(E.space(1),[13]) %James feature space            xy = xyrotate(DATA.spikefeatures(DATA.uid,p(1)),DATA.spikefeatures(DATA.uid,p(2)),angle);            cluster.space = [13 p(1:2)];        else            xy = xyrotate(DATA.pcs(:,p(1)),DATA.pcs(:,p(2)),angle);            cluster.space = [1 p(1:2)];        end        if forcesign == 1            id = find(xy(DATA.uid,1) > crit);            nid = find(xy(DATA.uid,1) <= crit);            cluster.sign = 1;        elseif forcesign == -1            id = find(xy(DATA.uid,1) < crit);            nid = find(xy(DATA.uid,1) >= crit);            cluster.sign = -1;        elseif crit < 0 && DATA.plottype == 1            id = find(xy(DATA.uid,1) < crit);            nid = find(xy(DATA.uid,1) >= crit);            cluster.sign = -1;        elseif crit <0 && mean(xy(:,1)) < 0            id = find(xy(DATA.uid,1) < crit);            nid = find(xy(DATA.uid,1) >= crit);            cluster.sign = -1;        else        id = find(xy(DATA.uid,1) > mean(exy(:,1)));        nid = find(xy(DATA.uid,1) <= mean(exy(:,1)));            cluster.sign = 1;        end        aid = id;        anid = nid;        id = DATA.uid(id);        nid = DATA.uid(nid);        e(1) = mean(DATA.energy(1,id));        e(2) = mean(DATA.energy(1,nid));        if e(2) > e(1) && forcesign == 0            cluster.sign = -cluster.sign;            cid = nid;            nid = id;            id = cid;            cid = anid;            anid = aid;            aid = cid;        end        if  length(id) < 2                    end        if quickmode.quickest == 1            cluster.dprime = NaN;            cluster.hdip = NaN;            cluster.dipsize = NaN;            cluster.autodipsize = NaN;            cluster.gmdprime = 0;        else%if quck ==2, do the 1d fit, so we can see fit quality                        if quickmode.fit1cut == 1 || quickmode.quick == 0            [a,b] = GMDip(xy(DATA.uid,:),DATA.energy(1,DATA.uid),'crit',crit,'label',DATA.idstr);            cluster.autodipsize = b.dipsize;            cluster.dipsize = b.cdipsize;            cluster.gmdprime = b.gmdprime;            cluster.gmfit1d = b.G{b.best};           end            if quickmode.quick == 0 || quickmode.fitallcut == 1                cluster.dprime = AllV.CalcDprime(xy(id,1),xy(nid,1));                cluster.hdip = HartigansDipTest(sort(xy(:,1)));                x = AllV.FitGaussMeans(xy(DATA.uid,:),2,'clusterid',aid);            end        end        r = xy(:,1);    else  %E.shape == 0   = ellpse        cluster.shape = 0;        xy = zeros(length(DATA.clst),2);        if E.space(1) == 6            xy(DATA.uid,:) = DATA.ndxy(DATA.uid,:);            cluster.space = E.space;        elseif isempty(p) %Var/E plot            xy(DATA.uid,2) = DATA.spkvar(ispk,DATA.uid)./DATA.energy(1,DATA.uid);            xy(DATA.uid,1) = DATA.energy(1,DATA.uid);            cluster.space = [5];        elseif ismember(E.space(1),[13]) %James feature space            xy(DATA.uid,1) = DATA.spikefeatures(DATA.uid,p(1));            xy(DATA.uid,2) = DATA.spikefeatures(DATA.uid,p(2));            cluster.space = [13 p(1:2)];        elseif ismember(E.space(1),[3 4])            xy(DATA.uid,1) = DATA.TemplateScores(DATA.uid,p(1));            xy(DATA.uid,2) = DATA.TemplateScores(DATA.uid,p(2));            cluster.space = [3 p(1:2)];        elseif E.space(1) == 2            AllVoltages = AllV.GetDataStruct(DATA, 'AllVoltages');            if DATA.verbose > 1                fprintf('%s AllVoltages %s(%d)\n',AllV.IDStr(DATA),sprintf('%d ',size(AllVoltages)),DATA.toplevel);            end            xy(DATA.uid,1) = AllVoltages(p(1),p(2),DATA.uid);            xy(DATA.uid,2) = AllVoltages(p(3),p(4),DATA.uid);            cluster.space = [2 p(1:4)];        elseif E.space(1) == 10            xy(DATA.uid,1) = DATA.icas(DATA.uid,p(1));            xy(DATA.uid,2) = DATA.icas(DATA.uid,p(2));            cluster.space = [E.space(1) p(1:2)];        else            xy(DATA.uid,1) = DATA.pcs(DATA.uid,p(1));            xy(DATA.uid,2) = DATA.pcs(DATA.uid,p(2));            cluster.space = [1 p(1:2)];        end        if isfield(E,'angle')            cluster.angle = E.angle;        else            cluster.angle = 0;        end        r = AllV.CalcClusterDistance(cluster, xy);        id = find(r(DATA.uid) < 1);        nid = find(r(DATA.uid) >1);        if quickmode.quick == 0 || quickmode.fitallcut == 1        x = AllV.FitGaussMeans(xy(DATA.uid,:),2,'clusterid',id);        end        if quickmode.quick == 0 || quickmode.fit1cut            if cluster.shape == 0                [a,b] = GMDip(AllV.Rprime(r(DATA.uid)),DATA.energy(1,DATA.uid),'crit',1,'label',DATA.idstr);            else                [a,b] = GMDip(r(DATA.uid),DATA.energy(1,DATA.uid),'crit',1,'label',DATA.idstr);            end            cluster.autodipsize = b.dipsize;            cluster.dipsize = b.cdipsize;            cluster.gmdprime = b.gmdprime;            cluster.gmfit1d = b.G{b.best};        end        id = DATA.uid(id);       nid = DATA.uid(nid);    end        %first reset any previous events to cluster 0    cl.clst(cl.clst == clnum) = 1;    %nid and clid refer to posttions in whole array, not just in uid    if length(DATA.uid) < length(cl.clst)        DATA.nid = nid;        DATA.clid = id;        X = setdiff(1:length(cl.clst),DATA.uid);        xid = find(cl.clst(X) > 0);        cl.clst(X(xid)) = 0;        cl.clst(id) = clnum;        if cluster.cluster == 1 %why did we do this - unsets clsuter 2 if its tehre            %            cl.clst(nid) = 1;        end    else        cl.clst(id) = clnum;        DATA.nid = nid;        DATA.clid = id;    end    DATA.clst = cl.clst;    cl.id = DATA.clid;    cl.nid = DATA.nid;        cluster.ncut = length(id);    cluster.clst = cl.clst;    if cluster.space(1) == 3        cluster.spacescale = DATA.TemplateScaling(cluster.space(2:end));    end        if quickmode.quick && quickmode.fitallcut ==0        cluster.bmc = 0;        cluster.mahal = [0 0 0 0];        if quickmode.fit2gauss            [dp, fits, dpdetails] = AllV.Fit2Gauss(cluster, r, DATA);            cluster.fitdprime = [dp dpdetails.fitpos dpdetails.dx dpdetails.minxpt];            cluster.fitdpparams = cat(1,fits{1}.params, fits{2}.params);        end    else        [dp, fits, dpdetails] = AllV.Fit2Gauss(cluster, r, DATA);        cluster.fitdprime = [dp dpdetails.fitpos dpdetails.dx dpdetails.minxpt];        cluster.fitdpparams = cat(1,fits{1}.params, fits{2}.params);        cluster.bmc = AllV.BimodalCoeff(r,1.5);        dip = MyDip(r);        cluster.mydipsize = dip.dipsize;        cluster.mydip = dip.x(dip.dip);;        a = AllV.FitGaussMeans(xy(DATA.uid,:),2);        cluster.mahal = [a.mahal a.dprime];        if isfield(cluster,'bestspace')            cluster.mahal(3) = cluster.bestspace(2);        else            cluster.mahal(3) = 0;        end        if ismember(E.shape,[0 1]) && x.obj.Converged >= 0            fprintf('Separation %.2f (manual) %.2f (auto)\n',x.mahal,a.mahal);            cluster.mahal(3) = x.mahal;            cluster.gmfit2dman = x.obj;        end        cluster.gmfit2d = a.obj;        if isfield(cluster,'gmdprime')            cluster.mahal(4) = cluster.gmdprime;        elseif isfield(E,'gmdprime')            cluster.mahal(4) = E.gmdprime;        elseif isfield(E,'gmfit1d') && isfield(E.gmfit1d,'mu')            cluster.mahal(4) = GMdprime(E.gmfit1d);        else            cluster.mahal(4) = 0;        end        cluster.times = DATA.t(id);        cluster.quick = 0;    end%don't rebuild MeanSpike if quick. But if there is no meanspike, it needs%buidling    if quickmode.quick == 0 || ~isfield(C,'MeanSpike') || quickmode.mean        cl.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc');        cluster.MeanSpike = cl.MeanSpike;        if DATA.plot.xyseq            C = cluster;            C.clst = cl.clst;            C.xy = xy;            AllV.PlotXYSequence(DATA,C);        end    end    if DATA.check.dropi(1) > 0        quickmode.dropi = 1;        quickmode.fit2gauss = 1;    end    cluster.quick = quickmode.quick; %if this is set, need to call again    cluster.r = r;    if quickmode.triggerhist        cluster =  AllV.PlotTriggerHist(DATA,cluster,quickmode);    end    cluster.minspke = prctile(DATA.energy(1,cl.id),1) .* 0.95;    cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),cl.nid),1) .* 0.95;    if recluster == 0  %If cluster  boundary has changed, will need new template%if reclsuter > 0 means this is applying an existing cluster, so if template has been %calcualted, no need to set this.        cluster.neednewtemplate = 1;    end    if isfield(DATA,'TemplateUsed') && ~isempty(DATA.TemplateUsed)         cluster.TemplateUsed = DATA.TemplateUsed;        if isfield(DATA,'DprimeUsed')            cluster.DprimeUsed = DATA.DprimeUsed;        end        if isfield(DATA,'mumeanUsed')            cluster.mumeanUsed = DATA.mumeanUsed;        end    end    if cluster.cluster > 1        clnum = cluster.cluster;        a = rmfields(cluster,'next'); %avoid recursion        cluster = DATA.cluster;        cluster.next{clnum-1} = a;    end     E.h = [];    if DATA.watchplots          if replot         DATA = AllV.ReplotPCs(DATA,cluster);        AllV.SetFigure(DATA.tag.spikes, DATA);        hold off;        DATA.spkst = DATA.uid;        AllV.QuickSpks(DATA, 1000);        if DATA.cluster.neednewtemplate == 0            AllV.PlotOneXY(DATA,DATA.xyplot.xy,'default');        end        if 1 % don't think we need this any more                    elseif (cluster.space(1) == 3) && ...                size(DATA.TemplateScores,2) > 12             AllV.SetFigure(DATA.tag.tmplscore, DATA);            subplot(1,1,1);            hold off;            plot(DATA.TemplateScores(cl.nid,1),DATA.TemplateScores(cl.nid,13),'.');            hold on;            plot(DATA.TemplateScores(cl.id,1),DATA.TemplateScores(cl.id,13),'r.');        elseif cluster.space(1) == 6            drawnow;            AllV.SetFigure(DATA.tag.tmplscore, DATA);            subplot(1,1,1);            AllV.PlotXY(DATA.ndxy,DATA.clst);            set(gca,'UserData',[NaN cluster.space]);        end        end        if DATA.plot.isi            AllV.PlotISI(DATA,0,0);        end    end    if DATA.auto.checkxcorr        cls = unique(DATA.clst);        if length(cls) > 2            AllV.CalcXcorr(DATA,[],'clusters');        end    elseif DATA.plot.xcorrprobes        AllV.CalcXcorr(DATA,[],'probes');    end    cluster.version = DATA.version;