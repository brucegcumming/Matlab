function res = SpoolFullV(V, varargin)
%AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%voltages, extracts segments triggered off one row, and plots PCS.
%try PCA on CSD/1st spatial derivative? 
%to build the V file, use MakeProbeIndex to build a list of files
% then use
%            FullV = PlotSpikeC(files,5,'probes',probelist,'sumv','submean','makev');
% to make FullV, the strucutre then used by AllVPcs
 
%BuildAllFullV calls this for all expts

% to fix try
% calculateing energy/spkvar for all probes is silly. Wastes memory and
% time.
%
%  Isues to explore. M170 Expt1 Bestspace mahal is 5 but and bmc is 00.32.
%  But bmc for PC 1,2 is 0.33.  Is this really the best space? check out
%  mahal distances for 1 and 2 D here. THE PC plot is a funny one witha
%  ring, and the cluster is stronges on other probes, so may not be soo
%  important.  
th(1) = 0;
setnspk = 0;
addmean = 1;

ispk = 1;
DATA.logfid = -1;
nprobepc = -1; % number of probes to include in in pc calculation
clplot = 0; %1 for density plot
plottype = 1; %plots made on first pass Could probably be 0 now...
plotv = 0;  %plot spikes superimposed on full voltage trace.
clusterprops = [];
spts = [-8:31];
smoothsd = 0; %smoothing just for trigger criterion
smoothv = 0; %smoothing of all traces for subsequent processing
tryall = 0;
DATA.dvdt = 0;
DATA.plotdvdt = 0;
DATA.autocutmode = 'mahal';
DATA.csd = 0;
DATA.name = [];
DATA.plotspk.probes = [];
DATA.plotspk.bytrial = 0;
DATA.plottype = 1;
DATA.tmpnorm = 1;
DATA.plotxcorr = 0;
DATA.usebmi = 0; %now do evertyhing on GM fit. Calculating old indices wastes a lot of time
DATA.lastcut = [];
DATA.SpaceTypes = {'Pcs' 'VarE'  'RawV' 'Template'};
DATA.elmousept.down = 0;
DATA.StdTemplate(1,:) = [   -0.0866   -0.0950   -0.1054   -0.1231   -0.1500   -0.2199   -0.4136   -0.7840 ...
      -1.2472   -1.6826   -1.8786   -1.6787   -1.2104   -0.6746   -0.2085    0.1502 ...
          0.4361    0.6637    0.8064    0.8326    0.7798    0.7036    0.6378    0.5895 ...
             0.5435    0.5011    0.4573    0.4182    0.3855    0.3550    0.3323...
             0.3 0.27 0.24 0.21 0.19 0.17 0.15 0.14 0.13]; 
DATA.StdTemplate(2,:) = [ -0.036 0.058 0.351 0.804 1.108 0.712 -0.372 -1.346 -1.804 -1.445 0.000 0.700 1.000 0.800 0.650 0.540 0.450 0.400 0.350 0.300 0.250 0.218 0.173 0.139 0.108 0.081 0.064 0.049 0.042 0.035 0.031 0.026 0.022 0.022 0.018 0.014 0.010 0.011 0.012 0.010 ];
DATA.usestdtemplates = 0;
thsign = 0;
calcpconly = 0;
calcclscores = 0;
muscale = 1;
newdata = 0;
addch = 0;
minenergy = 0;
minvar = 0;
oldscores = 0;
oldcluster = 0;
newdata = 0;
saveautocut = 0;
DATA.savespikes = 0;
DATA.watcharg = {};
DATA.watchplots = 0;
DATA.profligate = 0;
DATA.usegmcid = 0;
DATA.restricttimerange = [];
DATA.colors{1} = [0.5 0.5 0.5];
DATA.colors {2} = [1 0 0];
DATA.colors {3} = [0 1 0];
DATA.colors {4} = [0 0 1];
DATA.colors {5} = [0 1 1];
DATA.colors {6} = [1 1 0];
DATA.colors {7} = [1 0 1];
DATA.colors {8} = [0 1 0];
DATA.colors {9} = [0 1 0];
saveclusters = 0;
spkrate = 50;
autocutone = 0;
forcecluster = 0;
maxspksallowed = 600000;
recluster = [];
vt = [];
plotdprimemax = 0;  %old way to find boundaries. Really no good.
bmcrit = 0.21;
verbose = 1;
DATA.cstarttime = now;

spoolspikes = 0;
forcedrive = 'C:/bgc/data';
forcedrive = [];
forcename = [];
fullVname = [];
errs = {};
nerr = 0;

ttn = 1;
tt(ttn).time = now;
tt(ttn).str = 'Start';
ttn = ttn+1;
if length(varargin) && strcmp(varargin{end},'autocutall')
    autocutall = 1;
else
    autocutall = 0;
end

j = 1;
while j <= length(varargin)  %some varags must be parsed first
    if strcmp(varargin{j},'drive')
        j = j+1;
        if length(varargin{j}) <= 2
            if exist('Vall','var')
                if Vall.name(2) == ':';
                    Vall.name = [varargin{j} Vall.name(3:end)];
                else
                    Vall.name = [varargin{j} Vall.name];
                end
            end
        else
            forcedrive = varargin{j};
        end
    elseif strncmpi(varargin{j},'name',4)
        j = j+1;
        forcename = varargin{j};
    end
    j = j+1;
end

if ischar(V) 
    if exist(V,'file')
        fullVname = V;
       if ~autocutall
        tt(ttn).time = now;
        tt(ttn).str = sprintf('Loading %s',fullVname);
        ttn = ttn+1;

        if verbose
            tic;
            fprintf('Loading %s %s',fullVname,datestr(now,'HHMM:ss'));
        end
        load(V);
        if verbose
            fprintf(' took %.2f\n',toc);
        end
        maxl = size(FullV.V,2)-32;
        chspk = 1:size(FullV.V,1);
        V = FullV;
        clear FullV;
       else
           Vall = V;
       end
    else
    F = GetFigure(V);
    DATA = get(F,'UserData');
    vt = DATA.t;
    end
end
if isstruct(V)
    newdata = 1;
    Vall = V;
    DATA.name = Vall.name;
    if length(forcename)
        DATA.name = forcename;
    elseif length(forcedrive)
        DATA.name = regexprep(DATA.name,'[A-Z]:/Spike2/data',forcedrive);
    end
    DATA.Expt = [];
    if isfield(Vall,'exptno')
        DATA.Expt = LoadExpt(DATA,Vall.exptno);
        DATA.Expt.exptno = Vall.exptno;
        SetTrialList(DATA);
    end
    if isfield(Vall,'firstblk')
        DATA.Expt.blk = Vall.firstblk;
        DATA.Expt.exptno = DATA.Expt.exptno+0.1;
    else
        DATA.Expt.blk = 0;
    end
    DATA.starttime = Vall.start;
    if isfield(Vall,'samper')
    DATA.interval = Vall.samper;
    DATA.samplerate = 1./DATA.interval;
    else
    DATA.samplerate = 40000;
    DATA.interval = 1./DATA.samplerate;
    end
    if isfield(Vall,'blklen')
        DATA.blklen = Vall.blklen;
        DATA.blkstart = Vall.blkstart;
    end
    DATA.args = varargin;
    clear V;
elseif isfigure(V)
    DATA = get(V,'UserData');
    vt = DATA.t;
elseif ~isfield(DATA,'name')
    Vall.V = V;
    newdata = 1;
end



if exist('Vall','var') && ~ischar(Vall);
maxl = size(Vall.V,2)-32;
chspk = 1:size(Vall.V,1);
end
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'adcplot',3)
        DATA.plottype = 2;
    elseif strncmp(varargin{j},'autocutall',10)
        autocutall = 1;
    elseif strncmp(varargin{j},'autocut',7)
        autocutone = 1;
    elseif strncmp(varargin{j},'cutmode',7)
        j = j+1;
        DATA.autocutmode = varargin{j};
    elseif strncmp(varargin{j},'drive',5) %force drive letter
        j = j+1;
        if length(varargin{j}) == 1
        if Vall.name(2) == ':';
            Vall.name = [varargin{j} Vall.name(2:end)];
        else
            Vall.name = [varargin{j} Vall.name];
        end
        elseif length(varargin{j}) == 2
            if Vall.name(2) == ':';
                Vall.name = [varargin{j} Vall.name(3:end)];
            else
                Vall.name = [varargin{j} Vall.name];
            end
        else
            forcedrive = varargin{j};
        end
    elseif strncmp(varargin{j},'nomean',5)
        addmean = 0;
    elseif strncmp(varargin{j},'saveautocut',10)
        saveautocut = 1;
    elseif strncmp(varargin{j},'savespikes',6)
        DATA.savespikes = 1;
    elseif strncmp(varargin{j},'density',5)
        clplot = 1;
    elseif strncmpi(varargin{j},'clusters',4)
        if strncmpi(varargin{j},'clusterscores',10)
            calcclscores = 1;
        end
        if length(varargin) > j & iscell(varargin{j+1})
        j = j+1;
        Clusters = varargin{j};
        end
    elseif strncmpi(varargin{j},'setcluster',4)
        autocut = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j =j+1;
            clusterprops = varargin{j};
        end
    elseif strncmpi(varargin{j},'dvdt',4)
        DATA.dvdt = 1;
    elseif strncmpi(varargin{j},'accel',4)
        DATA.dvdt = 3;
    elseif strncmpi(varargin{j},'clear',4)
        muscale = 0.8;
    elseif strncmpi(varargin{j},'csda',4)
        DATA.csd = 2;
    elseif strncmpi(varargin{j},'csd',3)
        DATA.csd = 1;
    elseif strncmpi(varargin{j},'logfid',6)
        j = j+1;
        DATA.logfid = varargin{j};
    elseif strncmpi(varargin{j},'mine',4)
        j = j+1;
        minenergy = varargin{j};
    elseif strncmpi(varargin{j},'minvar',6)
        j = j+1;
        minvar = varargin{j};
    elseif strncmpi(varargin{j},'nowarn',6)
        warning('off','stats:gmdistribution:FailedToConverge');
        warning('off','stats:gmdistribution:MaxIterations');

    elseif strncmp(varargin{j},'tryall',5)
        tryall = 1;
    elseif strncmp(varargin{j},'tchan',5)
        j = j+1;
        ispk = varargin{j};
        newdata = 1;
        if length(ispk) > 1
            addch = 1;
            if strncmp(varargin{j-1},'tchanuse',8)
                addch = 2;
            end
        end
    elseif strncmp(varargin{j},'ichan',5)
        j = j+1;
        ispk = varargin{j};
        DATA.TemplateLabels = TemplateLabels(DATA,0);
    elseif strncmp(varargin{j},'muscale',6)
        j = j+1;
        muscale = varargin{j};
    elseif strncmp(varargin{j},'nprobepc',6)
        j = j+1;
        nprobepc = varargin{j};
    elseif strncmp(varargin{j},'name',4)
        j = j+1;
        DATA.name = varargin{j};
    elseif strncmp(varargin{j},'noplot',4)
        plottype = 0;
    elseif strncmp(varargin{j},'nspk',4)
        j = j+1;
        setnspk = varargin{j};
    elseif strncmp(varargin{j},'oldcl',4)
        oldcluster = 1;
        if  ~isempty(DATA.Clusters{DATA.probe(1)})
            C = DATA.Clusters{DATA.probe(1)};
            if C.mine > 0
                minenergy = C.mine;
                th = C.th;
            end
        end
    elseif strncmp(varargin{j},'oldt',4)
        oldscores = 1;
    elseif strncmp(varargin{j},'pcchan',5)
        j = j+1;
        chspk = varargin{j};
    elseif strncmp(varargin{j},'plotprobes',5)
        j = j+1;
        DATA.plotspk.probes = varargin{j};
    elseif strncmpi(varargin{j},'previous',4)
        j = j+1;
        DATA.lastcut = varargin{j};
    elseif strncmp(varargin{j},'plotv',5)
        plotv = 1;
        F = findobj('Tag','PCs','Type','Figure');
        if ~isempty(F)
            DATA = get(F,'UserData');
            PlotFullV(Vall.V, ispk, DATA);
            return;
        end
    elseif strncmp(varargin{j},'recut',5) %recut operates on whats in DATA
        DATA.cluster = DATA.Clusters{DATA.probe(1)};
        DATA.csd = DATA.cluster.csd;
        DATA.dvdt = DATA.cluster.dvdt;
        DATA = ReClassify(DATA,'newbound');
        set(DATA.toplevel,'UserData',DATA);
        return;
    elseif strncmp(varargin{j},'reclassify',6) %reclassify applies cluster space to new data (thr, crit, etc differen)
        recluster= 2;
        if length(varargin) > j && isfield(varargin{j+1},'dropi')
            j = j+1;
            DATA.cluster = varargin{j};
            forcecluster = 1;
        end
    elseif strncmp(varargin{j},'reapply',6) %applies cluster exactly, so should use same data
        recluster= 1;
        if length(varargin) > j && isfield(varargin{j+1},'dropi')
            j = j+1;
            DATA.cluster = varargin{j};
            forcecluster = 1;
        end
    elseif strncmp(varargin{j},'spool',4)
        spoolspikes = 1;
    elseif strncmp(varargin{j},'spts',4)
        j = j+1;
        spts = varargin{j};
    elseif strncmp(varargin{j},'vsmooth',4)
        j = j+1;
        smoothv = varargin{j};
    elseif strncmp(varargin{j},'smooth',4)
        j = j+1;
        smoothsd = varargin{j};
    elseif strncmp(varargin{j},'spkrate',4)
        j = j+1;
        spkrate = varargin{j};
    elseif strncmp(varargin{j},'th+',3)
        thsign = 1;
    elseif strncmp(varargin{j},'thboth',4)
        thsign = 2;
    elseif strncmp(varargin{j},'template',6) %reclassify applies cluster space to new data (thr, crit, etc differen)
       if strncmp(varargin{j},'templateshift',12) %move template to match probe
           forcecluster = 2;
       end
        j = j+1;
        recluster= 3;
        DATA.plottype = 3;
        if length(varargin{j}) == 1 && ~iscell(varargin{j})
            DATA.forceclusters{1} = varargin{j};
        else
            DATA.forceclusters = varargin{j};
        end
       if forcecluster == 2
        for k = 1:length(DATA.forceclusters)
           DATA.forceclusters{k}.MeanSpike.ms =  circshift(DATA.forceclusters{k}.MeanSpike.ms,ispk(1)-DATA.forceclusters{k}.probe);
        end
       end
        DATA.cluster = DATA.forceclusters{1};
        forcecluster = 1;
        plottype = 0;
    elseif strncmp(varargin{j},'bestspace',8)
        PCCluster(DATA, 0, 23);
        return;
    elseif strncmp(varargin{j},'autobestspace',8)
        PCCluster(DATA, 0, 26);
        return;
    elseif strncmp(varargin{j},'rpttemplate',6) %rebuild template scores, redo bestspace
        TemplatePlot(DATA);
        DATA = get(DATA.toplevel,'UserData');
        PCCluster(DATA, 0, 23);
        return;
    elseif strncmp(varargin{j},'saveclusters',10)
        saveclusters = 1;
    elseif strncmp(varargin{j},'templateline',10)
        PCCluster(DATA, 0, 6);
    elseif strncmp(varargin{j},'plottemplate',8)
        TemplatePlot(DATA);
    elseif strncmp(varargin{j},'threshold',3)
        j = j+1;
        th = varargin{j};
    elseif strncmp(varargin{j},'vpts',3)
        DATA.plottype =2;
    elseif strncmp(varargin{j},'watch',3)
        DATA.watcharg{1} = 'front';
        DATA.watchplots = 1;
    elseif strncmp(varargin{j},'xcorr',3)
        DATA.plotxcorr = 1;
    end
    j = j+1;
end

%want this after processing varargin so that can do autocutall
if autocutall
    if length(forcedrive) > 2
        Vall.name = regexprep(Vall.name,'[A-Z]:/Spike2/data',forcedrive);
    end
F = SetFigure('PCs');
res = AutoCutAll(ispk,  F, Vall, varargin(1:end-1));
return;
end

if ~isempty(vt) %Called with current figure, just to set a variable
    set(DATA.toplevel,'UserData',DATA);
    return;
end

if isfield(Vall,'t')
    vt = Vall.t;
elseif isfield(DATA,'blklen')
    first = 1;
    vt(sum(length(DATA.blklen))) = 0;
    for j = 1:length(DATA.blklen)
        last = first+DATA.blklen(j)-1;
        vt(first:last) = DATA.blkstart(j)+[1:DATA.blklen(j)].*DATA.interval;
        first = last+1;
    end
else
    vt= [1:size(Vall.V,2)] .* DATA.interval;
end

pres = {};
if length(ispk) > 4 && ~ autocutall
    SetFigure('PCs');
    args = {};
    clname = ClusterFile(DATA.name,DATA.Expt);
    if exist(clname,'file')
        load(clname);
        if calcclscores
            args = {args{:},'Clusterscores',Clusters};
        end
    end
    [nr,nc] = Nsubplots(length(ispk));
    for j = 1:length(ispk)
        res{j} = AllVPcs(Vall.V, varargin{:},'tchan',ispk(j),'noplot',args{:});
        if thsign == 2
            res{j}.pres = AllVPcs(Vall.V, varargin{:},'tchan',ispk(j),'th+','noplot',args{:});
        end
        GetFigure('PCs');
        subplot(nr,nc,j);
        id = [];
        PlotPCs(res{j}.pcs,1,2,clplot,res{j}.id,DATA.colors);
        title(sprintf('%.1f dp%.1f,%.1f',max(res{j}.dipvals)*100,res{j}.dp,res{j}.edp));
        drawnow;
    end
    if calcclscores
        hold off;
        for k =1:size(res{1}.Clusterscores,1)
            for j = 1:length(res)
                scores(k,j,1:size(res{j}.Clusterscores,2)) = smooth(res{j}.Clusterscores(k,:),50);
            end
        subplot(nr,nc,k);
        imagesc(squeeze(scores(k,:,:)));
        end
    end
    x = cat(1,res{:});
    GetFigure('Dprimes');
    hold off;
    plot([x.dp]);
    hold on;
    plot([x.edp],'r');
    if isfield(x,'pres')
        p = cat(1,x.pres);
        plot([x.dp],'o-');
        plot([x.edp],'ro-');
    end
    return;
end
if newdata && isfield(DATA,'name');
    DATA.msg = {};
    DATA.plotspk.submean = 0;
    DATA.plotspk.submax = 0;
    DATA.plotspk.submin = 0;
    DATA.plotspk.oneprobe = 0;
    DATA.plotspk.muscale = muscale;
    DATA.duration = size(Vall.V,2) .* DATA.interval;
    DATA.Clusters = {};
    if spkrate && setnspk == 0
        setnspk = round(DATA.duration .* spkrate);
    end
    
    afile = ClusterFile(DATA.name,DATA.Expt,'auto');
    DATA.probe = ispk;
    if exist(afile,'file')
        load(afile);
        DATA.AutoClusters = Clusters;
    else
        DATA.AutoClusters = {};
    end
    cfile = ClusterFile(DATA.name,DATA.Expt);
    if exist(cfile,'file')
        load(cfile);
        DATA.Clusters = Clusters;
        if isfield(DATA,'TemplateScores') %out of date
            DATA = rmfield(DATA,'TemplateScores');
        end
        if DATA.plottype == 3 & isfield(Clusters{DATA.probe(1)},'mean') %need to recalculate
            oldscores = 1;
        end
    end
    for j = 1:length(DATA.AutoClusters)
        if j > length(DATA.Clusters) ||  ~isfield(DATA.Clusters{j},'mahal')
            DATA.Clusters{j} = DATA.AutoClusters{j};
            DATA.Clusters{j}.auto = 1;
        end
    end
    if recluster
        if forcecluster == 0 %%used save cluster
            DATA.cluster = Clusters{DATA.probe(1)};
            if ~isfield(DATA.cluster,'minvar')
                DATA.cluster.minvar = 0;
            end
            th = DATA.cluster.Trigger;
        end
        tryall = 0;
        spts = DATA.cluster.spts;
        if isfield(DATA.cluster,'csd')
        DATA.csd = DATA.cluster.csd;
        DATA.dvdt = DATA.cluster.dvdt;
        end
        
        if recluster == 1 && DATA.cluster.minenergy > 0
            th = DATA.cluster.Trigger;
            minenergy = DATA.cluster.minenergy;
            minvar = DATA.cluster.minvar;
        end
    end
end




tt(ttn).times = now;
tt(ttn).strs = 'Start Trig';
ttn= ttn+1;
if verbose >1
    fprintf('Making Trig reference for %d %s\n',DATA.probe(1),datestr(now,'HHMM:ss'));
end
if length(ispk)  > 1 && addch
    rV = mean(Vall.V(ispk,:));
    if addch == 2
        Vall.V(ispk(1),:) = rV;
    end
else
    rV = Vall.V(ispk,:);
end

if length(smoothsd) > 1 %% explose effect of smoothing on the
if verbose >1
    fprintf('Smoothing  %s\n',datestr(now));
end
    colors = mycolors;
    GetFigure('Vhist');
    hold off;
    for j = 1:length(smoothsd)
        smv = smooth(Vall.V(ispk,:),smoothsd(j),'gauss');
        sgn = diff(sign(diff(smv,1,2)),1,2);
        id = find(sgn > 0 & V(ispk,2:end-1) < 0);
        [a,b] = smhist(Vall.V(ispk,id+1));
        plot(b,a,'color',colors{j});
        res.skew(j) = moment(Vall.V(ispk,id+1),3);
        hold on;
    end
    return;
elseif smoothsd
        smv = smooth(rV,smoothsd,'gauss');
else
    smv = rV;
end


if smoothv && newdata
    for j = 1:size(Vall.V,1)
        Vall.V(j,:) = smooth(Vall.V(j,:),smoothv,'gauss');
    end
end
for j = 1:size(rV,1)
sgn(j,:) = diff(sign(diff(smv(j,:),1,2)),1,2);
end
clear smv; 

tt(ttn).times = now;
tt(ttn).strs = 'Start Trig';
ttn= ttn+1;
if verbose > 1
    fprintf('%s %s\n',tt(ttn-1).strs,datestr(now,'HHMM:ss'));
end

for p = 1:size(Vall.V,1)
id = [];
autoth = 1;
rV = Vall.V(p,:);
sgn = diff(sign(diff(Vall.V(p,:),1,2)),1,2);
if th(1) < 0
    id = find(sgn(1,:) > 0)+1;
    xsd = std(rV(1,id));
    id = find(sgn(1,:) > 0 & rV(2:end-1) < th(1))+1;
elseif th(1) > 0
id = find(sgn(1,:) < 0 & rV(2:end-1) > th(1))+1;
else
    autoth = 1;
end
if isempty(id)
    if thsign == 1
    id = find(sgn(1,:) < 0)+1;
    prc = setnspk .* 100./length(id); % get 1000 spikes
    th(1) = prctile(rV(id),100-prc);
    id = id(rV(1,id+1) > th(1));
        xsd = std(rV(1,id));
    else
        id = find(sgn(1,:) > 0)+1;
        xsd = std(rV(1,id));
        prc = setnspk .* 100./length(id); % get 1000 spikes
        if prc > 100 %can happen if nspk > # minima
            th(1) = max(rV(1,id));
        else
            th(1) = prctile(rV(1,id),prc);
        end
        id = id(rV(1,id) < th(1));
        if size(rV,1) > 1
            xid = find(sgn(2,:) > 0)+1;
            prc = setnspk .* 100./length(xid); % get 1000 spikes
            th(1) = prctile(rV(2,xid),prc);
            xid = xid(rV(2,xid) < th(1));
            id = union(xid,id);
        end
    end
end
clear sgn;
id = id(id > -spts(1) & id < size(Vall.V,2)-spts(end));
ids{p} = id;
end
%remove any spikes at very beginning or end where there isn't enough
%data to include the whole spike

if th<0
res.xsd = xsd;
else
end
DATA.V = Vall.V;
DATA.vt = Vall.t;
DATA.ids = ids;
F = SetFigure('Spikes');
set(F,'UserData',DATA);
SpoolAllSpikes(DATA);
res.t = vt(id);
res.th = th;
DATA.Trigger = th;
DATA.autoth = autoth;
allid = [];

DATA.nprobes = size(Vall.V,1);


pcplots = [1 2; ...
           1 3; ...
           1 4; ...
           1 5; ...
           2 3; ...
           2 4; ...
           2 5; ...
           3 5];

    DATA.tmplots(1,:) = [1 3];
    DATA.tmplots(2,:) = [1 4];
    DATA.tmplots(3,:) = [1 2];
    DATA.tmplots(4,:) = [2 10];
    DATA.tmplots(5,:) = [1 8];
    DATA.tmplots(6,:) = [2 11];
    DATA.tmplots(7,:) = [2 12];
    DATA.tmplots(8,:) = [1 5];
    DATA.tmplots(9,:) = [2 13];
    DATA.tmplots(10,:) = [2 14];
    DATA.tmplots(11,:) = [2 15];
    DATA.tmplots(12,:) = [1 6];
    DATA.tmplots(13,:) = [6 11];
    DATA.tmplots(14,:) = [2 10];
    DATA.tmplots(15,:) = [2 10];
    DATA.tmplots(16,:) = [2 10];
    DATA.tmplots(17,:) = [1 2];
    DATA.tmplots(18,:) = [3 4];
    DATA.tmplots(19,:) = [1 3];
    DATA.tmplots(20,:) = [2 4];
    DATA.tmplots(21,:) = [1 4];
    DATA.tmplots(22,:) = [2 3];
    DATA.tmplots(23,:) = [1 2];
    DATA.tmplots(24,:) = [3 4];

    DATA.tmplspace(1,:)  = [1 8 10 5 12];

DATA.gmtypes = [1 0 2 3];
DATA.gmtypelabels = {'PCs', 'Var-E', 'ADC', 'Template'}; 


DATA.t = res.t;

DATA.spksperview = 100;
DATA.spts = spts;
DATA.spklst = 1:100;
DATA.pcprobes = chspk;

DATA.TemplateLabels = TemplateLabels(DATA,0);

DATA.toplevel = F;
set(F,'UserData',DATA);

function sz = memsize(X)
    x = whos('X');
    sz = x.bytes ./(1024 * 1024);
     

function S = SmallCluster(C)
%remove fields from C that use memory
S = C;
if isfield(S,'r')
    S = rmfield(S,'r');
end
if isfield(S,'xy')
    S = rmfield(S,'xy');
end

function PrintMsg(logfid, varargin)
                   
        s = sprintf(varargin{:});
        fprintf('%s\n',s);
        if logfid > 0
            fprintf(logfid,'%s\r\n',s);
        end

function res = AutoCutAll(ispk, toplevel, Vall, args)
    if length(ispk) < 3 && 0 %what was this for
        ispk = 1:size(Vall.V,1);
    end
    tstart = now;
    allcuts = {};
    logfile = ClusterFile(Vall.name,'log');
    logfid = fopen(logfile,'a');
    fprintf(logfid,'Start on %s at %s\r\n',Vall.name,datestr(now));
    for j = ispk
        istart = now;
        a =  AllVPcs(Vall ,args{:},'tchan',j,'logfid',logfid,'saveautocut');
        t = a.cluster.Trigger;
        spksd = std(a.cluster.MeanSpike.ms,0,2);
        [sv, svid] = max(spksd(a.chspk));
           a.cluster.good = GoodCluster(a.cluster);
           a.maxmean = a.chspk(svid);
           allcuts = {allcuts{:} a};
%
%  Explore lowering the trigger level if necessary. This can be for two
%  reasons:
%            1)  Dropping spikes. Produces a clipped distribution of
%            trigger point voltages, measured with dropi
%            2) cluster not dtopping too many spikes, but there are very
%            few non-cluster events, so the real cell is divided into tw
%            similar groups. NeedMore() Checks for this
%only explore lower triggers if the spike is biggest on this channel. Otherwise
%go to very low triggers to get all spikes that are really on another
%channel, and run out of memory
           evi = a.Evec.Eval(1)./sum(a.Evec.Eval);
           while NeedMore(a.cluster, a.Evec) && a.chspk(svid) == j && a.cluster.nspks < 500000
               a.cluster.needmore = NeedMore(a.cluster,a.Evec);
               PrintMsg(logfid, 'P%d: NeedMore%d %.2f(%.2f), bmi %.3f mahal%.2f, increasing events to %d\n',...
                   j,a.cluster.needmore,a.cluster.dropi(3),a.cluster.trigsd,a.cluster.bmc,a.cluster.mahal(1),a.cluster.nspks*2);
                 a = AllVPcs(Vall, args{:},'tchan',j, 'nspk',a.cluster.nspks * 2,'logfid', logfid, 'previous',PrevCluster(a.cluster),'saveautocut');
                 evi = a.Evec.Eval(1)./sum(a.Evec.Eval);
           end
           if GoodCluster(a.cluster)  && a.chspk(svid) == j
               dropi = a.cluster.dropi(3);
               nloop = 0;
           while ((a.cluster.dropi(3) < 2 && t < -a.xsd && GoodCluster(a.cluster) > 1 ...
                   && a.cluster.dropi(3) >= dropi)) && ...
                    a.maxspksused == 0 && nloop < 5
               nloop = nloop +1;
               dropi = a.cluster.dropi(3); %if dropi gets worse, not following a real cell.
               if a.cluster.dropi < 1 % just lower by 1SD. Estimate too bad to do more
                   t = a.cluster.Trigger + a.cluster.trigsd;
               else
                   t = a.cluster.Trigger + (2 - a.cluster.dropi(3)).*a.cluster.trigsd;
               end
               if t >= -a.xsd
                   t = -a.xsd; %effectively zero
               end
               PrintMsg(logfid,'P%d: DropSD %.2f(%.2f), bmi %.3f mahal%.2f, lowering threshold to %.3f, mine %.2f, took %.1f size %.1f\n',...
                   j,a.cluster.dropi(3),a.cluster.trigsd,a.cluster.bmc,a.cluster.mahal(1),t,a.cluster.minspke,mytoc(tstart),a.memsz);
               a.cluster.needmore = 4;
               a.cluster.first = PrevCluster(a.cluster);
               if isnan(a.cluster.minspke)
                   a.cluster.minspke = 0;
               end
               if a.cluster.space(1) == 6 && a.cluster.space(2) == 4
               a = AllVPcs(Vall, args{:},'tchan',j, 'thr',t,'mine',a.cluster.minspke,'minvar',a.cluster.minspkvar,'reclassify',a.cluster,'logfid',logfid,'saveclusters');
               else
               a = AllVPcs(Vall, args{:},'tchan',j, 'thr',t,'mine',a.cluster.minspke,'minvar',a.cluster.minspkvar,'logfid',logfid,'previous',a.cluster.first,'saveautocut');
               end
               spksd = std(a.cluster.MeanSpike.ms,0,2);
               [sv, svid] = max(spksd(a.chspk));
               a.cluster.good = GoodCluster(a.cluster);
               a.maxmean = a.chspk(svid);
               allcuts = {allcuts{:} a};
           end
           if nloop == 0
               fprintf('Good');
           end
           else
               fprintf('NoGood:');
           end
           a.cluster.good = GoodCluster(a.cluster);
           a.maxmean = a.chspk(svid);
           c = a.cluster;
           PrintMsg(logfid,'P%d(%d): %d/%d Spikes, dropi %.2f (%.2f) bmc %.3f, mahal %.1f,%.1f,%.1f G%d took %.2f size %.1f',...
               j,a.maxmean,c.ncut,c.nspks,c.dropi(3),c.dropi(4),c.bmc,c.mahal(1),c.mahal(2),c.mahal(4),GoodCluster(a.cluster),mytoc(istart),a.memsz);
           DATA.Clusters{j}.totaltime = mytoc(istart);
              
    end
        if logfid > 0
            fclose(logfid);
        end
        DATA = get(toplevel,'UserData');
        if DATA.plotxcorr
        a = PlotSpikeTimes(DATA.Clusters,'xcorr');
        for j = 1:length(DATA.Clusters)
%            DATA.Clusters{j}.synci = a.synci(j,:);
        end
        end
        SaveClusters(DATA,ClusterFile(DATA.name,DATA.Expt,'auto'));
        res.Clusters = DATA.Clusters;
        res.cuts = allcuts;
        
function first = PrevCluster(C)
    if isfield(C,'first')
        first.first = C.first;
    end
    if isfield(C,'needmore')
        first.needmore = C.needmore;
    end
    first.space = C.space;
    first.gmconverged =  C.gmfit1d.Converged;
    first.bestspace = C.bestspace;
    first.dropi = C.dropi;
    first.nspks = C.nspks;
    first.starttime = C.starttime;
    first.pcgms = C.pcgms;

        
function need = NeedMore(C, Evec)
%Cluster parameters can look bad if ALL of the triggered events are from a
%single cell, and tehre are no MU events. 
    if nargin > 1
        evi = Evec.Eval(1)./sum(Evec.Eval);
    else
        evi = 0;
    end

    szratio = std(C.MeanSpike.mu(C.probe,:))./std(C.MeanSpike.ms(C.probe,:));
    if C.dropi(3) > 1.7 && C.dropi(4) > 1
        need = 1;
    elseif C.dropi(3) > 1.9 && C.mahal(1) < 1.5 && szratio > 0.9
        need = 3;
    elseif C.MeanSpike.muxc > 0.8 && evi > 0.2
        need = 2;
    else
        need = 0;
    end
    
        
function C = CheckForMean(DATA,C)
        
    if ismember(C.space(1),[3 4]) | (C.space(1) ==6 && C.space(2) == 4) %need mean
        if ~isfield(C,'MeanSpike')
            if isempty(DATA.clid) && isfield(C,'r') && isfield(C,'sign')
                DATA.clid = find(C.r .* C.sign > C.crit(1) .* C.sign);
                DATA.nid = setdiff(1:size(DATA.AllV,3),DATA.clid)';
            end
            C.MeanSpike = PlotMeanSpike(DATA,'recalc');
        end
    end
    if length(C.mahal) > 3 & ~isfield(C,'gmdprime')
        C.gmdprime = C.mahal(4);
    end

function DATA = ReClassify(DATA, varargin)
    newbound = 0;
    fittemplate = 0;
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'newboundary',4)
            newbound = 1;
        elseif strncmpi(varargin{j},'template',4)
            fittemplate = 1;
        end
        j = j+1;
    end
    
    DATA.cluster = CheckForMean(DATA,DATA.cluster);
    oldms = DATA.cluster.MeanSpike.ms;
    if fittemplate
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
        [a,b, DATA.xy, details] = TemplateSpace(DATA,'template');
        fprintf('Template GM %.2f 6D %.2f 2d, %.2f1D\n',a(1),a(2),details.gmdprime)
        DATA.cluster.gmfit = b;
        DATA.cluster.gmdprime = details.gmdprime;
        DATA.cluster.sign = details.cluster.sign;
        DATA.cluster.crit = details.cluster.crit;
        DATA.cluster.shape = 2;
        DATA.cluster.angle = 0;
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.xy = DATA.xy; %best ND;
        DATA.cluster.xy(:,3) = DATA.TemplateScores(:,2); %sum of r
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.bestcl = details.bestcl;
        DATA.cluster.usegmcluster = details.usegmcluster;

    elseif ismember(DATA.cluster.space(1),[3 4]) & ~isfield(DATA,'TemplateScores')
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
   
    elseif DATA.cluster.space(1) == 6 && DATA.cluster.space(2) == 4
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
        DATA.plottype = 3;
        [a,b, DATA.xy, details] = BestSpace(DATA,'template');
        DATA = ReplotPCs(DATA,[]);
        DATA.cluster.gmfit = b;
    elseif DATA.cluster.space(1) == 6
        DATA.xy = ProjectND(DATA,DATA.cluster.space(2),DATA.cluster.gmfit);
    end
        if newbound
%            [dip, details] = FindDip(DATA.xy(:,1),DATA.energy(DATA.probe(1),:),'gmix');
            [dip, details] = GMDip(DATA.xy,DATA.energy(DATA.probe(1),:));
            DATA.cluster.crit = dip(1);
            DATA.cluster.sign = details.sign;
            DATA.cluster.gmdipres = details.dipres;
            DATA.cluster.gmdprime = details.gmdprime;
            DATA.cluster.autodipsize = details.dipsize;
            DATA.cluster.dipsize = details.cdipsize;
        end
    [cl, cluster, DATA.xy] = ClassifySpikes(DATA,DATA.cluster);
    if newbound || fittemplate%if boundary changed, need to record new criterion, threshold/dropi etc.
        DATA.cluster = cluster;
    end
    DATA.Cluster{DATA.probe(1)}.mahal = cluster.mahal;
    DATA.clusterboundary = BoundaryFromCluster([],cluster);
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    DATA.MeanSpike = cl.MeanSpike;
    DATA.cluster.MeanSpike = cl.MeanSpike;
    DATA.cluster.minspke = prctile(DATA.energy(DATA.probe(1),DATA.clid),1) .* 0.95;
    DATA.cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),DATA.clid),1) .* 0.95;
    if fittemplate
        DATA.cluster.xy = DATA.xy;
        DATA.cluster.xy(:,3) = DATA.TemplateScores(:,2);
        xc = corrcoef(oldms(:),DATA.cluster.MeanSpike.ms(:));
        DATA.cluster.templatexc = xc(1,2);
    end
    
function DATA = CalcTemplatesFromMean(DATA, MeanSpike);
        
            
    DATA.TemplateLabels = TemplateLabels(DATA,0);
    Scores = CalcScores(DATA,MeanSpike);
    ispk = find(DATA.chspk == DATA.probe(1));

    if isfield(DATA,'TemplateScores')
    oldscores = DATA.TemplateScores;
    end
    if size(Scores,1) > 1
        DATA.TemplateScores(:,1)= Scores(ispk,1,:);
        DATA.TemplateScores(:,8)= Scores(ispk,2,:);
    end
    DATA.TemplateScores(:,2)= sum(Scores(:,1,:));
    DATA.TemplateScores(:,3)= Scores(1,1,:);

    if size(Scores,2) > 3
        DATA.TemplateScores(:,5)= Scores(ispk,4,:);
    end

    if ispk > 1 && ispk == size(Scores,1);
        DATA.TemplateScores(:,4)= Scores(end-1,1,:);
    else
        DATA.TemplateScores(:,4)= Scores(end,1,:);
    end
    DATA.TemplateScores(:,9) = Scores(ispk,3,:);
    DATA.TemplateScores(:,10)= sum(Scores(:,2,:));
    DATA.TemplateScores(:,12)= sum(Scores(:,3,:));
    DATA.TemplateScores(:,11)= sum(Scores(:,5,:));
    DATA.TemplateScores(:,6)= Scores(ispk,5,:);
    DATA.TemplateScores(:,7)= Scores(1,5,:);
    if DATA.tmpnorm
    for j = 1:size(DATA.TemplateScores,2)
        DATA.TemplateScores(:,j) = DATA.TemplateScores(:,j)./std(DATA.TemplateScores(:,j));
    end
    end
%    DATA.tmpdips = CalculateTemplateDips(DATA);
    DATA.tmpdips = zeros(1,8);

    
    
function tt = TimeMark(tt, str)

    ttn = length(tt)+1;
    tt(ttn).time = now;
    tt(ttn).str = str;


function good = GoodCluster(C)
    good = 0;
    if C.bmc > 0.3 || C.mahal(1) > 6
        good = 3;
    elseif C.bmc > 0.23 || C.mahal(1) > 3
        good = 2;
    elseif C.bmc > 0.21 || C.mahal(1) > 2
        good = 1;
    end

function C = ClusterFromBoundary(E, C)
    C.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
    exy = xyrotate(E.pos([1 3]),E.pos([2 4]),C.angle);
    C.crit = mean(exy(:,1));
    C.len = abs(diff(E.pos([1 3])) + i * diff(E.pos([2 4])));
if isfield(E,'shape')
    C.shape = E.shape;
end
        
function E = BoundaryFromCluster(E, C)
if isfield(C,'space')
    E.pcplot = C.space(2:end);
end
if isfield(C,'y')
    y = C.y;
elseif isfield(E,'pos')
    l = abs(diff(E.pos([1 3])) + i * diff(E.pos([2 4])));
    y = [-l/2 l/2];
elseif isfield(C,'len')
    l = C.len;
    y = [-l/2 l/2];
else
y = [-1 1];
end
if C.shape == 2 %line
E.xyr(1) = C.crit(1) .* cos(-C.angle);
E.xyr(2) = C.crit(1) .* sin(-C.angle);
exy = xyrotate([C.crit(1) C.crit(1)],y,-C.angle);
E.pos(1) = exy(1,1);
E.pos(2) = exy(1,2);
E.pos(3) = exy(2,1);
E.pos(4) = exy(2,2);
elseif isfield(C,'xyr')
    E.pos(1) = C.xyr(1) - C.xyr(3);
    E.pos(2) = C.xyr(2) - C.xyr(4);
    E.pos(3) = C.xyr(1)+ C.xyr(3);
    E.pos(4) = C.xyr(2) +C.xyr(4);
end


angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));

if isfield(C,'sign')
E.sign = C.sign;
else
    E.Sign = 0;
end
if isfield(C,'shape')
    E.shape = C.shape;
    if C.shape == 2
        E.space = C.space;
    elseif C.shape == 1 %% not sure if this needs difffernt treatment
        E.space = C.space;
    end
end
if ~isfield(E,'h')
    E.h = [];
end
if ~isfield(E,'shape')
    E.shape = 1;
end
f = {'angle' 'firstspace' 'auto', 'took','bestspace','ctime','xyr', 'bestcl' 'usegmcluster' 'gmdprime' 'gmdipres' 'gmfit' 'gmfit1d' 'bestll'};
for j = 1:length(f)
    if isfield(C,f{j})
        E.(f{j}) = C.(f{j});
    end
end

function [C, Evec, pcs, dip, chspk, errs] = CalcPCs(DATA, nprobepc)

errs = {};
uid = DATA.uid;
if length(nprobepc) > 1 || nprobepc(1) >= 0
    if length(nprobepc) > 1
        chspk = DATA.probe(1) + nprobepc;
    else
        chspk = DATA.probe(1)-nprobepc:1:DATA.probe(1)+nprobepc;
    end
    if DATA.csd == 2
        chspk = chspk-2;
        chspk = chspk(chspk > 0 & chspk < size(DATA.AllV,1)-1);
        if isempty(chspk)
            chspk = DATA.probe(1);
        end
    elseif DATA.csd
        chspk = chspk-1;
        chspk = chspk(chspk > 0 & chspk < size(DATA.AllV,1));
    else
        chspk = chspk(chspk > 0 & chspk <= size(DATA.AllV,1));
    end
end


if DATA.csd 
    if DATA.csd == 2
    csd = diff(DATA.AllV,2,1);
    else
    csd = diff(DATA.AllV,1,1);
    end
    TV = csd(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,csd(chspk(j),:,uid));
    end
elseif DATA.dvdt == 2
    TV = DATA.AllV(chspk(1),:,uid);
    TV = cat(2,TV,diff(DATA.AllV(chspk(1),:,uid),1,2));
    for j = 2:length(chspk)
        TV = cat(2,TV,DATA.AllV(chspk(j),:,uid));
        TV = cat(2,TV,diff(DATA.AllV(chspk(j),:,uid),1,2));
    end
elseif DATA.dvdt == 3
    TV = diff(DATA.AllV(chspk(1),:,uid),2,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(DATA.AllV(chspk(j),:,uid),2,2));
    end
elseif DATA.dvdt
    TV = diff(DATA.AllV(chspk(1),:,uid),1,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(DATA.AllV(chspk(j),:,uid),1,2));
    end
else
    TV = DATA.AllV(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,DATA.AllV(chspk(j),:,uid));
    end
end
TV = squeeze(TV)';

C = cov(TV);
[pc, E] = eig(C);
[a,b] = max(diag(E));
Evec.Eval = diag(E);
if b == 1
    fprintf('Max Eigernvector First\n');
    errs = {errs{:} 'Max Eigernvector First\n'};
else
    pc = fliplr(pc); %put largest first;
    Evec.Eval = Evec.Eval(end:-1:1);
end
pcs = TV*pc;
Evec.Evec = pc;
pcs = pcs;

if nargout > 3
    if DATA.usebmi
    for j = 1:10
%        dip(j) = HartigansDipTest(sort(pcs(:,j)));
        dip(j) = BimodalCoeff(pcs(:,j),1.5);
    end
    p = DATA.pcplots;
    for j =1:length(p)
        [as(j),bs(j)] = BestAngle(pcs(:,p(j,1)),pcs(:,p(j,2)),1);
    end
    dip = bs;
    [P, dip(length(p)+1)] = GMfit(pcs(:,1:4),2,1);
    else
    [P, dip(1)] = GMfit(pcs(:,1:4),2,1);
    end
end


function PlotFullV(V, ispk, DATA)

SetFigure('FullV');
hold off;
y = 0;
for p = 1:length(ispk)
plot(V(ispk(p),:)+y);
hold on;
id = DATA.clid;
spts = DATA.spts;
for j = 1:length(id);
    plot([DATA.t(id(j)):DATA.t(id(j))+length(spts)-1]+spts(1),DATA.AllV(ispk(p),:,id(j))+y,'r');
end
nid = DATA.nid;
for j = 1:length(nid);
    plot([DATA.t(nid(j)):DATA.t(nid(j))+length(spts)-1]+spts(1),DATA.AllV(ispk(p),:,nid(j))+y,'g');
end
y = y + max(V(ispk(p),:));
end
hold off;


function distance  = gmdistance(G)

    D = mahal(G,G.mu);
    distance = sqrt(2./((1./D(1,2))+(1./D(2,1))));

function [theta, c, details] = BestGMAngle(x,y, test, varargin)
  a = 0:pi/36:pi;
  if diff(size(x)) > 0
      G = gmdistribution.fit([x' y'],2,'Options',statset('MaxIter',1000));
  else
      G = gmdistribution.fit([x y],2,'Options',statset('MaxIter',1000));
  end
  for j = 1:length(a)
      xy = xyrotate(x,y,a(j));
      G = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',1000)); %2 Gaussians, 1 dimension
      details.mahal(j) = abs(diff(G.mu))./sqrt(mean(G.Sigma));
  end
 [c,id] = max(details.mahal);
 details.mahal2d = gmdprime(G);
 theta = a(id);
 b = theta-pi/36:pi/360:theta+pi/36;
 for j = 1:length(b)
      xy = xyrotate(x,y,b(j));
      G = gmdistribution.fit(xy(:,1),2,'Options',statset('MaxIter',1000)); %2 Gaussians, 1 dimension
      finemahal(j) = abs(diff(G.mu))./sqrt(mean(G.Sigma));
 end
mahals = cat(2,details.mahal, finemahal);
[details.angles,id] = sort(cat(2,a,b));
details.mahal = mahals(id);
[c,id] = max(details.mahal);
theta = details.angles(id);
  
  
function [theta, c, details] = BestAngleGM(xy, G, varargin)
%
% Find best angle to cut in a 2D space, using 1-D Gaussian Mixtures

  a = 0:pi/36:pi;

  if isempty(G)
      G = GMFit(xy,2,1);
  end
  id = cluster(G, xy);
  
  c = 0;
  for j = 1:length(a)
      XY = xyrotate(xy(:,1),xy(:,2),a(j));
      [aa,bb] = GMDip(XY,[],'idlist',id);
      details.d(j) = bb.mahal(bb.best);
      if details.d(j) > c
          c = details.d(j);
          bestfit = bb.G{bb.best};
      end
  end
  [c, j] = max(details.d);
  details.besti = j;
  theta = a(j);
  details.angles = a;
  details.gmfit = bestfit;

      details.xy = xyrotate(xy(:,1),xy(:,2),theta);
      details.dip = HartigansDipTest(sort(details.xy(:,1)));


    function [theta, c, details] = BestAngle(x,y, test, varargin)
%
% Find best angle to cut in a 2D space. With very skewed distributison
% (Energy, ADC value where it is clipped), the bimodality coeffeicient is
% misleading. But otherwise itts smoother and more reliable. 

if bitand(test,8)
    domix = 1;
else
    domix = 0;
end
  a = 0:pi/36:pi;
  
  details.mahal = NaN;
  if domix
  if std(x-y) > std(x+y)/10000 && std(y) > 0 && std(x) > 0
  if diff(size(x)) > 0
  G = gmdistribution.fit([x' y'],2,'Options',statset('MaxIter',1000));
  else
  G = gmdistribution.fit([x y],2,'Options',statset('MaxIter',1000));
  end
  details.mahal = gmdistance(G);
  end
  end
  
  for j = 1:length(a)
      xy = xyrotate(x,y,a(j));
      if bitand(test,2)
      dip(j) = HartigansDipTest(sort(xy(:,1)));
      end
      if bitand(test,4)
%      [aa,bb] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
      [aa,bb] = GMDip(xy,[]);
      mydip(j) = bb.dipsize(1);
      end
      skews(j) = skewness(xy(:,1));
      kurts(j) = kurtosis(xy(:,1));
      coeff(j) = BimodalCoeff(xy(:,1),1.5);
  end
  
  details.bmc = coeff;
  if bitand(test,1)
      details.coeff = coeff;
  elseif bitand(test,2)
      details.hdip = dip;
      details.coeff = dip;
  end
  [c, j] = max(details.coeff);
  details.besti = j;
  theta = a(j);
  details.angles = a;
  
  %if using bimodality coeff to find best angle, calc Hartigan for best
  %angle so can compare with other measures.
  if test == 1
      xy = xyrotate(x,y,theta);
      details.dip = HartigansDipTest(sort(xy(:,1)));
  end

 function CheckClusters(Clusters, str)

     for j = 1:length(Clusters)
         if ~isfield(Clusters{j},'mahal') %empty really
             missing(j) = 1;
         else
             missing(j) = 0;
         end
         if ~isfield(Clusters{j},'MeanSpike')
             missing(j) = 2;
         end
     end
 nerr = sum(missing);
 if nerr
     errordlg(sprintf('Missing %d Clusters at %s',nerr,str));
 end
 nm = sum(missing == 2);
 if nm
     errordlg(sprintf('Missing %d MeanSpikes at %s',nerr,str));
 end
            
            
  
function [Clusters, ClusterDetails, id] = SaveClusters(DATA, outname)
    oldname = get(DATA.toplevel,'name');
    dname = strrep(outname,'.mat','Details.mat');
    set(DATA.toplevel,'name',sprintf('Saving %s',outname));
    drawnow;
    CheckClusters(DATA.Clusters,'Save');
    if exist(outname,'file')
        load(outname);
        for j = 1:length(DATA.Clusters);
            if j > length(Clusters) || isempty(Clusters{j})
                Clusters{j} = DATA.Clusters{j};
            elseif ~isfield(Clusters{j},'mahal') %empty really
                Clusters{j} = DATA.Clusters{j};
            elseif isfield(Clusters{j},'mean') %old
                Clusters{j} = rmfield(Clusters{j},'mean');
            elseif isfield(Clusters{j},'r') 
                Clusters{j} = rmfield(Clusters{j},'r');
            elseif isfield(Clusters{j},'clst') 
                Clusters{j} = rmfield(Clusters{j},'clst');
            end
            if ~isfield(Clusters{j},'probe')
                Clusters{j}.probe = j;
            end
            if ~isfield(Clusters{j},'mahal') %still empty
                if isfield(DATA,'AutoClusters') && length(DATA.AutoClusters) >= j
                    Clusters{j} = DATA.AutoClusters{j};
                    fprintf('Cluster %d was empty - reverting to Auto\n',j);
                    errordlg(sprintf('Cluster %d was empty Reloaded AutoCluster\n',j));                    
                else
                    errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',j));                    
                end
                    
            end
        end
    end
    if exist(dname,'file')
        load(dname);
    end
    p = DATA.probe(1);
    
    Clusters{DATA.probe(1)}.times = DATA.t(DATA.uid);
    Clusters{DATA.probe(1)}.dpsum = sum(abs(DATA.MeanSpike.vdprime(DATA.chspk,:)),2);
    if ~isempty(DATA.restricttimerange)
        Clusters{DATA.probe(1)}.restricttimerange = DATA.restricttimerange;
    end
    if DATA.savespikes
        Clusters{DATA.probe(1)}.spkfile = SpkFilename(DATA);
    end
    ClusterDetails{DATA.probe(1)}.Evec = DATA.Evec;
    f = fields(DATA.cluster);
    for j = 1:length(f)
        if strcmp(f{j},'r') || strcmp(f{j},'clst')  %xy, clst go into ClusterDetails
        else
        Clusters{DATA.probe(1)}.(f{j}) = DATA.cluster.(f{j});
        end
    end
    Clusters{p}.clst  = DATA.clst(DATA.uid);
    Clusters{DATA.probe(1)}.auto = 0;
    Clusters{p}.pcmean = mean(DATA.pcs(:,1:4)); %to check for sign reversal
    savexy = 2;
    
    Clusters{DATA.probe(1)}.savetime(1) = now;
    if savexy && isfield(DATA,'xy')
        if savexy ==2
            ClusterDetails{DATA.probe(1)}.xy = DATA.xy(DATA.uid,:);
            ClusterDetails{DATA.probe(1)}.t = DATA.t(DATA.uid);
            ClusterDetails{DATA.probe(1)}.clst = DATA.clst(DATA.uid);
            id = DATA.uid;
        else
        minpts = min([ceil(DATA.duration * 50) size(DATA.AllV,3)]);
        minpts = min([minpts length(DATA.uid)]);
        if length(DATA.clid) > length(DATA.uid)/2
            npts = length(DATA.xy);
        elseif length(DATA.clid) < length(DATA.uid)/10;
            npts = length(DATA.uid)/5;
        else
            npts = length(DATA.clid).*2;
        end
        npts = max([ceil(npts) minpts]);
        if mean(DATA.xy(DATA.clid,1)) > mean(DATA.xy(DATA.nid,1)) 
            [a,id] = sort(DATA.xy(DATA.uid,1),'descend');
        else
            [a,id] = sort(DATA.xy(DATA.uid,1));
        end        
        id = sort(id(1:npts));
        if length(id) > size(DATA.AllV,3);
            id = id(id < size(DATA.AllV,3));
        end
        id = DATA.uid(id);
        ClusterDetails{DATA.probe(1)}.xy = DATA.xy(id,:);
        ClusterDetails{DATA.probe(1)}.t = DATA.t(id);
        ClusterDetails{DATA.probe(1)}.clst = DATA.clst(id);
        end
        ClusterDetails{DATA.probe(1)}.ctime = Clusters{DATA.probe(1)}.ctime;
    else
        id = [];
    end
    Clusters{DATA.probe(1)}.savetime(2) = now;
    nerr = 0;
    for j = 1:length(Clusters)
        if ~isfield(Clusters{j},'mahal')
            nerr = nerr+1;
        end
    end
    if nerr > 0
        errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',j));
    end

    save(outname,'Clusters');
    save(dname,'ClusterDetails');
    Clusters{DATA.probe(1)}.savetime(3) = now;
    set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));

    
function C = StripClusters(Clusters)
        
for j = 1:length(Clusters)
    C{j} = Clusters{j};
    if isfield(C{j},'xy')
        C{j} = rmfield(C{j},'xy');
    end
end



function name = SpkFilename(DATA)
    [a,b] = fileparts(DATA.name);
    if isdir(DATA.name)
        a = DATA.name;
    end
    spkdir = [a '/Spikes'];
    xs='';
    if rem(DATA.Expt.exptno,1) > 0.001
        xs = 'a';
    end
    name = sprintf('%s/%s.p%dt%d%s.mat',spkdir,b,DATA.probe(1),floor(DATA.Expt.exptno),xs);
        
function SaveSpikes(DATA, id)
    name = SpkFilename(DATA);
    a = fileparts(name);
    if ~exist(a,'dir')
        mkdir(a);
    end
    Spikes.values = squeeze(DATA.AllV(DATA.probe(1),:, id));
    if size(Spikes.values,2) > 100
        Spikes.values = Spikes.values';
    end
    Spikes.times = reshape(DATA.t(id),length(id),1);
    xy = DATA.xy(id,1);
    Spikes.codes = zeros(size(xy));
    if isfield(DATA,'clst')
        Spikes.codes = DATA.clst-1;
    elseif DATA.cluster.sign < 0
        Spikes.codes(xy < DATA.cluster.crit) = 1;
    else
        Spikes.codes(xy > DATA.cluster.crit) = 1;
    end
    Spikes.codes = reshape(Spikes.codes,length(Spikes.codes),1);
    Spikes.Header.ctime = now;
    save(name,'Spikes');
    
function tcut = IsTemplateCut(E)    
   tcut = 0;
   if isfield(E,'space') && E.space(1) ==6 && E.space(2) == 4
       tcut = 1;
   elseif E.plottype == 3
       tcut = 1;
   end
    
function [E, cluster] = CutAndSave(DATA, varargin)

    nosave = 0;
    args = {};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'nosave',4)
            nosave = 1;
        else
            args = {args{:} varargin{j}};
        end
        j = j+1;
    end
    [E, Scores, tmpdips, xy] = AutoCut(DATA, args{:});

    if ~isempty(Scores) && (~isfield(DATA,'TemplateScores') || E.plottype == 3)
        DATA.TemplateScores = Scores;
        DATA.tmpdips = tmpdips;
    elseif E.newscores
        DATA = get(DATA.toplevel,'UserData');
    end
    DATA.xy = xy;
    DATA.cboundary = E;
    if IsTemplateCut(E)
        DATA.plottype = 3;
    else
        DATA.plottype = E.plottype;
    end
    if DATA.watchplots
    PlotHistogram(DATA,E);
    end
    [cl, DATA.cluster, DATA.xy] = ClassifySpikes(DATA,E);
    if  ~isfield(DATA.cluster,'dropi')
        fprintf('no dropi calculated');
    end
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    DATA.MeanSpike = cl.MeanSpike';
    DATA.cluster.MeanSpike = DATA.MeanSpike;
    DATA.cluster.errs = DATA.errs;
    DATA.cluster.starttime = DATA.cstarttime;
    DATA.Clusters{DATA.probe(1)} = DATA.cluster;
    if nosave == 0
        outname = ClusterFile(DATA.name,DATA.Expt,'auto');
        [DATA.Clusters, DATA.ClusterDetails, DATA.savespkid] =  SaveClusters(DATA, outname);
        if length(DATA.savespkid) > length(DATA.xy)
            fprintf('Save Spike Mismatch %d vs %d\n',length(DATA.savespkid),length(DATA.xy));
        end
    end
    set(DATA.toplevel,'UserData',DATA);
    cluster = DATA.cluster;
    
    
function [distance, obj, xy, details] = TemplateSpace(DATA, varargin)
    recalc = 0;
    nc = 2;
    ntr = 1;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'recalc',5)
            recalc = 1;
        end
        j = j+1;
    end
    if isfield(DATA,'clid')  && ~isempty(DATA.clid)
        idlist = ones(1,length(DATA.TemplateScores));
        idlist(DATA.clid) = 2;
    else
        idlist = [];
    end
    if isfield(DATA.cluster,'space') && DATA.cluster.space(1) == 3
        space = DATA.cluster.space;
    else
        space = [3 1 8];
    end
    [objs{2}, distance(2), a] = GMfit(DATA.TemplateScores(:,space(2:3)),nc,ntr,'idlist',idlist);
        [objs{1}, distance(1), a] = GMfit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),nc,ntr,'idlist',idlist);
%        objs{5} = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),3,'Options',statset('MaxIter',1000));
        details.Converged(1) = objs{1}.Converged;
        if ~isfield(DATA,'clid') || isempty(DATA.clid) || recalc
            xy = ProjectND(DATA, 4, objs{1});
%        [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
        [C.crit, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
        details.usegmcluster = 0;
        C.sign = details.sign;
        if details.sign >= 0
            DATA.clid = find(xy(:,1) > C.crit(1));
            DATA.nid = find(xy(:,1) <= C.crit(1));
        else
            DATA.clid = find(xy(:,1) < C.crit(1));
            DATA.nid = find(xy(:,1) >= C.crit(1));
        end
        else
            C.sign = 0;
        end
        obj = objs{1};
        details.cluster = C;
        details.Gs = objs;
        if distance(2) > distance(1)
        details.bestcl = cluster(objs{2},DATA.TemplateScores(:,space(2:3)));
        else
        details.bestcl = cluster(objs{1},DATA.TemplateScores(:,DATA.tmplspace(1,:)));
        end
        if max(distance) > details.gmdprime * 1.5 %1D cut is very poor
            details.usegmcluster = 1;
        end
        if DATA.logfid > 0
            fprintf(DATA.logfid,'P%d, T%d at %s\r\n',DATA.probe(1),DATA.cluster.probe,datestr(now));
        end
 
function [distance, obj, xy, details] = BestSpace(DATA, varargin)

    newtemplate = 0;
    nloops = 0; %to test multiple fits for consistency
    nr = 1;
    ntr = 1;
    nc = 2;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'newtemplate',7)
            newtemplate = 1;
            nr=5;
        elseif strncmpi(varargin{j},'template',5)
            ntr=5; %redoing a template fit. Make sure its right.
        elseif strncmpi(varargin{j},'ncells',5)
            j = j+1;
            nc = varargin{j};
        end
        j = j+1;
    end
    p= DATA.probe(1);
    uid = DATA.uid;
    %when doing this from scratch (remaking template from other cut), use
    %several starts to make sure the PC cut is best. 
    %Now uses GMfit that set a sensible start point.
    [P, distance(1), a] = GMfit(DATA.pcs(uid,DATA.pcspace),nc,1);
    objs{1} = P;
    details.Converged(1) = P.Converged;
    ll(1) = P.NlogL;
    xy(:,1) = DATA.energy(p,uid);
    xy(:,2) = DATA.spkvar(p,uid)./DATA.energy(p,uid);
    
    [E, distance(2)] = GMfit(xy,nc,1);
    details.Converged(2) = E.Converged;
    objs{2} = E;
    D = mahal(E,E.mu);
    distance(2) = sqrt(1./(2./D(1,2)+1./D(2,1)));
    ll(2) = E.NlogL;
    
    [V, distance(3), a] = GMfit(squeeze(DATA.AllV(p,DATA.vspace,uid))',nc,1);
    objs{3} = V;
    details.Converged(3) = V.Converged;

    ll(3) = V.NlogL;
    quick = 1;
    if newtemplate || ~isfield(DATA,'TemplateScores')
        [d, best] = max(distance);
        [xy, details.cid] = ProjectND(DATA, best, objs{best});
%        [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
%could use gmdistribution/cluster to assign to groups here.
        if quick & length(unique(details.cid)) > 1
            DATA.clid = find(details.cid  == 2);
            DATA.nid = find(details.cid  == 1);
        else
        [cluster.crit, details] = GMDip(xy(uid,:),DATA.energy(DATA.probe(1),uid));
        cluster.sign = details.sign;
        if details.sign >= 0
            DATA.clid = find(xy(:,1) > cluster.crit(1));
            DATA.nid = find(xy(:,1) <= cluster.crit(1));
        else
            DATA.clid = find(xy(:,1) < cluster.crit(1));
            DATA.nid = find(xy(:,1) >= cluster.crit(1));
        end
        end
        DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
        TemplatePlot(DATA,'nodip','usemean');
        DATA = get(DATA.toplevel,'UserData');
        imean = DATA.MeanSpike;
        ntr = 1;
    end
    if isfield(DATA,'TemplateScores')

        [objs{4}, distance(4), a] = GMfit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,ntr);
%        objs{5} = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),3,'Options',statset('MaxIter',1000));
        D = mahal(objs{4},objs{4}.mu);
        distance(4) = sqrt(2./(1./D(1,2)+1./D(2,1)));
        ll(4) = objs{4}.NlogL;
        details.Converged(4) = objs{4}.Converged;

        tic;
        for j = 1:nloops
            T{j}= gmdistribution.fit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,'Options',statset('MaxIter',1000));
            ll(j) = T{j}.NlogL;
        end
        if nloops
            toc
            tic;
            Tn= gmdistribution.fit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,'Options',statset('MaxIter',1000),'Replicates',nloops);
            toc;
        end
%        [a,b,c] = cluster(objs{4},DATA.TemplateScores(:,DATA.tmplspace(1,:)));
    end
    
    bestll = min(ll);
    [d, best] = max(distance);
    if best == 3 && ll(3) > 2 * bestll
        safeid = [1 2 4];
        [d, best] = max(distance(safeid));
        best = safeid(best);
        err = sprintf('Using Space %d (%.1f, %.1f) not Space 3 (%.1f. %.1f)\r\n',best,distance(best),ll(best),distance(3),ll(3));
        fprintf('%s',err);
        if DATA.logfid > 0
            fprintf(DATA.logfid,'%s',err);
        end
    end
    if best == 4
        DATA.cluster.space = [6 4];
        if ~isfield(DATA,'clid') || isempty(DATA.clid)
        xy = ProjectND(DATA, best, objs{best});
%        [cluster.crit, details] =
%        FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
        [cluster.crit, details] = GMDip(xy(uid,:),DATA.energy(DATA.probe(1),uid));
        cluster.sign = details.sign;
        if details.sign >= 0
            DATA.clid = find(xy(:,1) > cluster.crit(1));
            DATA.nid = find(xy(:,1) <= cluster.crit(1));
        else
            DATA.clid = find(xy(:,1) < cluster.crit(1));
            DATA.nid = find(xy(:,1) >= cluster.crit(1));
        end
        end
        olddistance = distance;

        objs{4} = IterateTemplateFit(DATA, objs{best});
        D = mahal(objs{4},objs{4}.mu);
         distance(4) = sqrt(2./(1./D(1,2)+1./D(2,1)));
        DATA = get(DATA.toplevel,'UserData');
    end
    
    obj = objs{best};
   [xy, details.cid] = ProjectND(DATA, best, objs{best});
   e(1) = mean(DATA.energy(p,details.cid == 1));
   e(2) = mean(DATA.energy(p,details.cid == 2));
   if e(1) > e(2)
       details.cid = 3 - details.cid;
   end


%    [a,b,c] = BestAngle(xy(:,1),xy(:,2), 3); %Should not be necessary.
    details.ll = ll;
%    details.bestangle = a;
%    xy = xyrotate(xy(:,1),xy(:,2),a);


function [G, D, all] = GMfit(X, nd, nr, varargin)
    idlist = [];
    if nr < 1
        nr = 1;
    end
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'idlist',5)
            j =j+1;
            idlist = varargin{j};
        end
        j = j+1;
    end
%include my own starting point
      C = cov(X);
      pvar = diag(C);
        [E,V] = eig(C);
        pc = E(:,end);
        pcb = E(:,end-1);
        pc = pc./max(abs(pc));
        pcb = pcb./max(abs(pcb));
        for j = 1:size(X,2)
            S.mu(1,j) = mean(X(:,j)) + pc(j);
            S.mu(2,j) = mean(X(:,j)) - pc(j);
            if nd == 3
                S.mu(3,j) = mean(X(:,j)) + pcb(j);
            end
        end
        for j = 1:nd
        S.Sigma(:,:,j) = C./sqrt(2);
        end
        all.guess = S;

    
    for j = 1:nr
        try
       all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000));
       all.ds(j) = gmdprime(all.Gs{j});
        catch
            fprintf('GM Fit fail\n');
            all.Gs{j} = S;
            all.Gs{j}.Converged = -1;
            all.Gs{j}.NlogL = NaN;
        end
    end
    j = j+1;
    try
    all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000),'Start',S);
    all.ds(j) = gmdprime(all.Gs{j});
    catch
        fprintf('GM Fit (My Start) fail\n');
        all.Gs{j} = S;
        all.Gs{j}.Converged = -1;
        all.Gs{j}.NlogL = NaN;
        all.ds(j) = 0;
    end
    if length(idlist) == length(X)
        j = j+1;
        try
            id = find(idlist == 2);
            nid = find(idlist==1);
            S.mu(1,:) = mean(X(nid,:));
            S.mu(2,:) = mean(X(id,:));
            S.Sigma(:,:,1) = cov(X(nid,:));
            S.Sigma(:,:,2) = cov(X(id,:));
            S.PComponents(1) = length(nid)./size(X,1);
            S.PComponents(2) = length(id)./size(X,1);
            all.Gs{j} = gmdistribution.fit(X,nd,'Options',statset('MaxIter',1000),'Start',S);
            all.ds(j) = gmdprime(all.Gs{j});
        catch
            fprintf('GM Fit (Start with Classification) fail\n');
            all.Gs{j} = S;
            all.Gs{j}.Converged = -1;
            all.Gs{j}.NlogL = NaN;
            all.ds(j) = 0;
        end
    end
    if all.ds(j) > max(all.ds(1:j-1)) * 1.1
        [a,b] = max(all.ds(1:j-1));
           fprintf('Best manual start %.2f vs %.2f (%d of %d)\n',all.ds(j),a,b,j-1);
    end
    [D,b] = max(all.ds);
    G = all.Gs{b};


function G = IterateTemplateFit(DATA, G)
j = 1;
xc = [0 0];
DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
T{1} = DATA.MeanSpike.ms;
nc = size(G.mu,1);
    while j < 10 && xc(1,2) < 0.5
    xy = ProjectND(DATA, 4, G);
%    [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
    [cluster.crit, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
    cluster.sign = details.sign;
    if details.sign >= 0
        DATA.clid = find(xy(:,1) > cluster.crit(1));
        DATA.nid = find(xy(:,1) <= cluster.crit(1));
    else
        DATA.clid = find(xy(:,1) < cluster.crit(1));
        DATA.nid = find(xy(:,1) >= cluster.crit(1));
    end
    DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
    DATA.TemplateScores = TemplatePlot(DATA,'nodip');
    T{j+1} = DATA.MeanSpike.ms;
    [G, distance(j),a] = GMfit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),nc,1);
%    G = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),2,'Options',statset('MaxIter',1000));
%    D = mahal(G,G.mu);
%    distance(j) = sqrt(2./(1./D(1,2)+1./D(2,1)));
    xc = corrcoef(T{j+1}(:),T{j}(:));
    xcs(j) = xc(1,2);
    j = j+1;
    end
    set(DATA.toplevel,'UserData',DATA);
        
function [xy, cid] = ProjectND(DATA, best, obj)
    pnorm = 0;
    if size(obj.mu,1) > 2
        v = range(obj.mu);
    else
        v = diff(obj.mu);
    end
        if best == 1
        sd = std(DATA.pcs(:,[1 2 3 4]));
        v = v./sd;
        if pnorm
        S(:,1) = DATA.pcs(:,1) ./ sd(1);
        S(:,2) = DATA.pcs(:,2) ./ sd(2);
        S(:,3) = DATA.pcs(:,3) ./ sd(3);
        S(:,4) = DATA.pcs(:,4) ./ sd(4);
        else
            S = DATA.pcs(:,1:4);
        end
        xy(:,1) = v * S';
    elseif best == 2
        xy(:,1) = DATA.energy(DATA.probe(1),:);
        xy(:,2) = DATA.spkvar(DATA.probe(1),:)./DATA.energy(DATA.probe(1),:);
        S = xy;
        xy(:,1) = v * S';
    elseif best == 3
        S = squeeze(DATA.AllV(DATA.probe(1),DATA.vspace,:))';
        xy(:,1) = v * S';
    elseif best == 4
        if pnorm
        sd = std(DATA.TemplateScores(:,DATA.tmplspace(1,:)));
        v = v./sd;
        for j = 1:size(DATA.tmplspace,2)
            S(:,j) = DATA.TemplateScores(:,DATA.tmplspace(1,j)) ./ sd(j);
        end
        else
            S = DATA.TemplateScores(:,DATA.tmplspace(1,:));
        end
        xy(:,1) = v([1 5]) * S(:,[1 5])';
        xy(:,1) = v * S';
    end
    
    ov = v;
    ov(1) = v(2);
    ov(2) = -v(1);
    if length(v) > 3
    ov(3) = v(4);
    ov(4) = -v(3);
    end
    if length(v) > 10 %make this 5 to activate
    ov(5) = v(6);
    ov(4) = -v(5);
    end
    xy(:,2) = ov * S';
    cid = cluster(obj, S);% mahal distance is unsigned, not so useful
    if length(cid) > length(DATA.uid)
        id = setdiff(1:length(cid),DATA.uid);
        cid(id) = 0;
    end
    if length(unique(cid)) == 1
        fprintf('GM clustering only 1 group\n');
    end
   
    
function [E, Scores, tbs, xy, details]  = AutoCut(DATA, varargin)
usev = 0;
refine = 0;
usegm = 0;
tbs = [];
Scores = [];
newtemplate = 0;
newDATA = 0;

if strcmp(DATA.autocutmode,'mahal')
    usegm = 1;
end
E.cutmode = DATA.autocutmode; 
E.newscores = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'refine',4)
        refine  = 1;
    elseif strncmpi(varargin{j},'mahal',4)
        usegm = 1;
        E.cutmode = 'mahal';
    end
j = j+1;
end
ctime = now;
if usegm
    [bd,obj,xy,c] = BestSpace(DATA,'newtemplate');
    E.newscores = 1;
    [a,b] = max(bd);
    E.bestspace = [a b];
    E.bestll = c.ll;
    E.bestcl = c.cid;
    E.bestfitres = c.Converged;
    if b == 4
        [gm,d]= max(bd(1:3));
        s = sprintf('GM %.2f for space %d (%.2f for %d)dt%d csd%d',a,b,gm,d,DATA.dvdt,DATA.csd);
    else
        s = sprintf('GM %.2f for space %s dt%d csd%d',a,DATA.SpaceTypes{b},DATA.dvdt,DATA.csd);
    end
    fprintf('%s\n',s);
    if DATA.logfid > 0
        fprintf(DATA.logfid,'P%d (S%d) %s %d events at %s\r\n',DATA.probe(1),DATA.savespikes,s,size(DATA.AllV,3),datestr(now));
    end
%    [dip, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:),'gmix');
    [dip, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
   [E.gmfit2d, d] = GMfit(xy,2,1,'idlist',c.cid);
   if d > 2  && d > details.gmdprime * 1.5 %% 2D significantly better - check rotation
       [a, gm,  dipres] = BestAngleGM(xy, E.gmfit2d);
       if gm > 2 && gm > details.gmdprime * 1.5
           xy = dipres.xy; %rotated values
           details.dipres = dipres.gmfit;
           details.gmdprime = gm;
           [dip, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
       end
   end
    
    E.space = [6 b];
    E.bestd = bd;
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.gmdip = dip;
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    E.pos = [crit min(xy(:,2)) crit max(xy(:,2))];
    E.plottype = DATA.plottype;
    E.gmfit = obj;
    E.gmfit1d = details.G{details.best};
    E.gmdipres = details.dipres;
    E.gmdprime = details.gmdprime;
    E.autodipsize = details.dipsize;
    E.dipsize = details.cdipsize;
    details.newDATA = 1;
    return;
end
    p = DATA.pcplots;
    for j =1:length(p)
        [as(j),bs(j),c] = BestAngle(DATA.pcs(:,p(j,1)),DATA.pcs(:,p(j,2)),1);
        gd(j) = c.mahal;
    end
    n = length(as);
    [x,y] = GetClusterXYData(DATA,[]);
    n = n+1;
    vare = n;
    [as(vare), bs(vare),c] = BestAngle(x,y,1);
    bs(vare) = bs(vare).* 0.7;  %% only use varE if substantially better
    gd(vare) = c.mahal;
    n = n+1;

    [bd,obj,xy] = BestSpace(DATA);
    bs(n) = BimodalCoeff(xy(:,1));
    [gd(n), besttype] = max(bd);
    as(n) = 0;
    E.bestspace(1) = bs(n);
    E.bestspace(2) = gd(n);
    E.bestd = bd;
    n = n+1;

    
    if usev
    p = DATA.vpts;
    for j =1:length(p)
        [as(j+n),bs(j+n),c] = BestAngle(DATA.AllV(p(j,1),p(j,2),:),DATA.AllV(p(j,3),p(j,4),:),2);
        bs(j+n) = c.bmc(c.besti);
    end
    end
    [a,j] = max(bs);
    pcbii = a;
    cluster.angle = as(j);
    if j == vare %x,y already made
        p = DATA.vpts;
        cluster.space = [];
        E.plottype = DATA.plottype;
    elseif j == vare+1  %BestSpace
        p = DATA.vpts;
        cluster.space = [6 besttype];
        cluster.angle = 0;
        E.plottype = DATA.plottype;
        E.shape = 2;
        E.space = cluster.space;
        x = xy(:,1);
        y = xy(:,2);
    elseif j > size(DATA.pcplots,1)
        j = j-8;
        p = DATA.vpts;
        cluster.space = [2 DATA.vpts(j,:)];
        x = DATA.AllV(p(j,1),p(j,2),:);
        y = DATA.AllV(p(j,3),p(j,4),:);
        E.plottype = 2;
    else
        p= DATA.pcplots;
        cluster.space = [1 DATA.pcplots(j,:)];
        x = DATA.pcs(:,p(j,1));
        y = DATA.pcs(:,p(j,2));
        E.plottype = 1;
    end
    xy = xyrotate(x,y,cluster.angle);
%    [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
    [cluster.crit, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
    cluster.sign = details.sign;
    if refine
        if j == 9  %used var/e
            bettercrit = 1;
        else
            bettercrit = 1.2; %prefer PC cut if equal
        end
        if details.sign >= 0
        DATA.clid = find(xy(:,1) > cluster.crit(1));
        DATA.nid = find(xy(:,1) <= cluster.crit(1));
        else
        DATA.clid = find(xy(:,1) < cluster.crit(1));
        DATA.nid = find(xy(:,1) >= cluster.crit(1));
        end
        Scores = TemplatePlot(DATA);
        E.newscores = 1;
        cluster.firstspace = cluster.space;
        cluster.firstbmi = bs(j);
        if length(bd) < 4  %no template for the BestSpace calc; Do again
            DATA.TemplateScores = Scores;
            [bd, obj, bxy] = BestSpace(DATA);
            [gd(vare+1), besttype] = max(bd);
            bmi =  BimodalCoeff(bxy(:,1));
            E.bestspace(1) = bmi;
            E.bestspace(2) = gd(vare+1);
            E.bestd = bd;
            bs(vare+1) = bmi;
            if bmi > pcbii
                xy = bxy;
                cluster.space = [6 besttype];
                cluster.angle = 0;
                E.plottype = DATA.plottype;
                E.shape = 2;
                E.space = cluster.space;
                x = xy(:,1);
                y = xy(:,2);
                pcbii = bmi;
%                [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
                [cluster.crit, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
                cluster.sign = details.sign;
            end
        end
        p = DATA.tmplots;
        for j =1:8
            [tas(j),tbs(j),c] = BestAngle(Scores(:,p(j,1)),Scores(:,p(j,2)),1);
            tgd(j) = c.mahal;
        end
        [a,j] = max(tbs);
        if a > pcbii * bettercrit
            cluster.angle = tas(j);
            cluster.firstspace = cluster.space;
            cluster.space = [3 p(j,:)];
            x = Scores(:,p(j,1));
            y = Scores(:,p(j,2));
            E.plottype = 3;
            E.shape = 1;
            xy = xyrotate(x,y,cluster.angle);
%            [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
              [cluster.crit, details] = GMDip(xy,DATA.energy(DATA.probe(1),:));
            cluster.sign = details.sign;
        end
    end
    cluster.auto = 1;
    E.autotook = mytoc(ctime);
    if length(cluster.space) > 1 && cluster.space(1) == 6 && cluster.space(2) == 4
        E.plottype = 3;
    end
    details.newDATA = newDATA;
    E = BoundaryFromCluster(E,cluster);
    

function C = CutAndPlot(x,y, energy)
    [as,bs] = BestAngle(x,y,1);
    C.angle = as;
    xy = xyrotate(x,y,C.angle);
%    [crit,b] = FindDip(xy(:,1),energy,'plot');
    [crit,b] = GMDip(xy,energy,'plot');
    C.crit = crit(1);
    C.hdip = as;
    C.bmc = BimodalCoeff(xy(:,1),1.5);
   hold off; 
   id = find(xy(:,1) > crit(1));
   nid = find(xy(:,1) <= crit(1));
   plot(x(id),y(id),'r.','markersize',1);
   hold on;
   plot(x(nid),y(nid),'.','markersize',1);
   smw = round(size(xy,1)./500);
   [pdf, xv] = smhist(xy(:,1),smw);
   if mean(xv < 0)
       xv = -xv;
   end
   xv = (xv-min(xv)) .* range(x)./range(xv);
   plot(xv, pdf .* max(y)./max(pdf));
   title(sprintf('Dips %.4f,%.3f',C.hdip,C.bmc));
    
function C = OptimizeVarE(DATA)

    
   GetFigure('VarE');
   subplot(2,2,1);
    x = DATA.energy(DATA.probe(1),:);
    y = DATA.spkvar(DATA.probe(1),:)./DATA.energy(DATA.probe(1),:);
    Cs(1) = CutAndPlot(x,y,DATA.energy(DATA.probe(1),:));
    drawnow;

    subplot(2,2,2);
    x = DATA.energy(DATA.probe(1),:).^2;
    y = DATA.spkvar(DATA.probe(1),:).^2 ./DATA.energy(DATA.probe(1),:).^2;
    Cs(2) = CutAndPlot(x,y,DATA.energy(DATA.probe(1),:));
    drawnow;
    C = Cs(1);
    C.sign = 0;
    


    
function C = OptimizeBoundary(DATA)
    
    space = DATA.cluster.space;
    if space(1) == 5
        C = OptimizeVarE(DATA);
        return;
    elseif space(1) == 6  %% cut in > 2 dimensions.
        x = DATA.xy(:,1);
        y = DATA.xy(:,2);
    else
    xi = space(2);
    yi = space(3);
    end
    if space(1) == 1
        x = DATA.pcs(:,xi);
        y = DATA.pcs(:,yi);
    elseif space(1) == 2
        xi = space(3);
        yi = space(5);
        x = squeeze(DATA.AllV(space(2),xi,:));
        y = squeeze(DATA.AllV(space(4),yi,:));
    elseif space(1) == 3
        x = DATA.TemplateScores(:,xi);
        y = DATA.TemplateScores(:,yi);
    end
    a = -pi/2:pi/36:pi/2; % use this range because this is what atan returns;
    
    for j = 1:length(a)
        xy = xyrotate(x,y,a(j));
        dip(j) = HartigansDipTest(sort(xy(:,1)));
        [aa,bb] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:));
        mydip(j) = bb.dipsize(1);
        coeff(j) = BimodalCoeff(xy(:,1),1.5);
    end
    if space(1) == 2
        [dipval,b] = max(dip);
        dipval = coeff(b);
    else
        [dipval,b] = max(coeff);
    end
    C.angle = a(b);
    xy = xyrotate(x,y,C.angle);
    

    [crit,b] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:),'gmix','plot');
    C.crit = crit(1);
    C.sign = b.sign;
    GetFigure('Dips')
    hold off;
    plot(a,dip./max(dip));
    hold on;
    plot(a,mydip./max(mydip),'r');
    plot(a,coeff./max(coeff),'g');

function c = BimodalCoeff(x, e)
    e = 1.3;
    c = (1+skewness(x).^2)./((kurtosis(x).^e)+3);
    
function ProbeMenu(a,b, fcn)

    onoff = {'off' 'on'};
[DATA, F] = GetDataFromFig(a);
if ismember(fcn,[1 2 3])
    [C, DATA.Evec, DATA.pcs, DATA.dipvals, DATA.chspk] = CalcPCs(DATA,fcn-1);
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
elseif fcn == 4
    DATA.dvdt = ~DATA.dvdt;
    set(a,'Checked',onoff{DATA.plotdvdt});
end
    
    
function res = PlotClusters(a,b,fcn)
[DATA, F] = GetDataFromFig(a);
SetFigure('Clusters');
res = [];
if fcn == 1
    PlotSpikeTimes(DATA.Clusters);
elseif fcn == 2
    res = PlotSpikeTimes(DATA.Clusters,'xcorr');
    for j = 1:length(DATA.Clusters)
        DATA.Clusters{j}.synci = res.synci(j,:);
    end
    DATA.xcorrs = res;
    set(DATA.toplevel,'UserData',DATA);
elseif fcn == 3
    PlotSpikeTimes(DATA.Clusters,'probequality');
end

    
function PCCluster(a,b, fcn)

[DATA, F] = GetDataFromFig(a);
onoff = {'off' 'on'};

if ismember(fcn,[1 2 3 5 6])  %if intereactive, look at plots when change cuts
    DATA.watchplots = 1;
end
if fcn == 1
    if isstruct(a)
        E = b;
    else
    SetEllipseDrawing(DATA, 0); %called from gui
    return;
    end
    DATA.usegmcid = 0;
elseif fcn == 32
    if isstruct(a)
        E = b;
    else
        SetEllipseDrawing(DATA, 0,'cluster',2); %called from gui
    return;
    end
    DATA.usegmcid = 0;
elseif fcn == 2
    it = GetFigure('VarE',DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r');
    E.pcplot = [];

elseif fcn == 3 %replot
    PlotMeanSpike(DATA);
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 4 %toggle density plot
    DATA.clplot = ~DATA.clplot;
    set(a,'Checked',onoff{DATA.clplot+1});
    ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif ismember(fcn, [5 32])
    if fcn == 5
        TemplatePlot(DATA);
    else
        TemplatePlot(DATA,'stdtemplate','nodip');
    end
    DATA = get(DATA.toplevel,'UserData');
    DATA.plottype = 3;
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 6
    it = GetFigure('TemplateScores',DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    E.pcplot = get(gca,'UserData');
    PlotHistogram(DATA,E);
    return;
elseif fcn == 7
    it = GetFigure('VarE',DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    E.pcplot = [];
    DATA.cboundary = PlotHistogram(DATA,E);
elseif ismember(fcn, [8 11])
    DATA.watchplots = 1;
    it = GetFigure('PCs',DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    if isempty(E)
        return;
    end
    E.pcplot = get(gca,'UserData');
    DATA.cboundary = PlotHistogram(DATA,E);
elseif fcn == 9
    DATA.plottype = 1;
    set(DATA.toplevel,'UserData',DATA);
    DATA = ReplotPCs(DATA,[]);
    return;
elseif fcn == 10
    DATA.plottype = 2;
    set(DATA.toplevel,'UserData',DATA);
    DATA = ReplotPCs(DATA,[]);
    return;
elseif ismember(fcn,[12 25])
    oldname = get(DATA.toplevel,'name');
    outname = ClusterFile(DATA.name,DATA.Expt);
    [DATA.Clusters, DATA.ClusterDetails, id] = SaveClusters(DATA,outname);
    if fcn == 25 %save Spikes
        SaveSpikes(DATA,id);
    end
    set(DATA.toplevel,'UserData', DATA);
    return;

elseif fcn == 13
    if ~isfield(DATA,'Vall') %can't do this from Gui
        return;
    end
    if DATA.plottype > 2
        DATA.plottype = 1;
        set(DATA.toplevel,'UserData',DATA);
    end
    name = get(DATA.toplevel,'Name');
    set(DATA.toplevel,'Name',sprintf('Triggering on Probe %d',DATA.probe(1)));
    drawnow;
    AllVPcs(DATA.toplevel,'tchan',DATA.probe(1)+1,DATA.args{3:end});
    set(DATA.toplevel,'Name',name);
    return;
elseif fcn == 14
    DATA.plottype = 3;
    set(DATA.toplevel,'UserData',DATA);
    if isfield(DATA,'TemplateScores') && sum(DATA.TemplateScores(:,12)) > 0
    DATA = ReplotPCs(DATA,BoundaryFromCluster([],DATA.cluster));
    else
    TemplatePlot(DATA);
    end
    return;
elseif fcn == 15
    DATA.plottype = 4;
    set(DATA.toplevel,'UserData',DATA);
    DATA = ReplotPCs(DATA,[]);
    return;
elseif fcn == 17
    DATA.Clusters{DATA.probe(1)} = [];
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 18  %optimnize boundary line in current plot
    C = OptimizeBoundary(DATA);
%    DATA.Clusters{DATA.probe(1)}.angle = C.angle;
%    DATA.Clusters{DATA.probe(1)}.crit = C.crit;
    DATA.cluster.angle = C.angle;
    DATA.cluster.crit = C.crit;
    DATA.cluster.sign = C.sign;
    E = BoundaryFromCluster([],DATA.cluster);
    [cl, cluster, DATA.xy] = ClassifySpikes(DATA,E);
    E.sign = cluster.sign;
    DATA.cluster = rmfield(cluster,'r');
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    PlotHistogram(DATA,E);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 19  %make Automatic cut
    set(DATA.toplevel,'Name','Thinknig');
    drawnow;
    DATA.watchplots = 1; 
    [E, scores, dips, xy, details] = AutoCut(DATA,'refine');
    if details.newDATA
        DATA = get(DATA.toplevel,'UserData');
    end
    DATA.xy = xy;
%If template scores have been calculated, keep them
    if ~isempty(scores)
        DATA.TemplateScores = scores;
        DATA.tmpdips = dips;
    end
    set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));
    DATA.plottype = E.plottype;
    DATA.cboundary = E;
    PlotHistogram(DATA,E);
elseif fcn == 20
    DATA.plottype = 6;
    set(DATA.toplevel,'UserData',DATA);
    DATA = ReplotPCs(DATA,[]);
    return;
elseif fcn == 21
    if isfield(DATA,'Clusters') && length(DATA.Clusters) >= DATA.probe(1) && ~isempty(DATA.Clusters{DATA.probe(1)})
        DATA.cluster = DATA.Clusters{DATA.probe(1)};
        DATA = CheckTemplates(DATA,DATA.cluster);
        if DATA.cluster.space(1) == 6
                   [DATA. xy, DATA.gmcid] = ProjectND(DATA, DATA.cluster.space(2), DATA.cluster.gmfit);
                   DATA.cluster.bestcl = DATA.gmcid;
        end
        DATA.watchplots = 1;  %must be in gui
        [cl, cluster] = ClassifySpikes(DATA,DATA.cluster);
        DATA.cluster.MeanSpike = cl.MeanSpike;
        DATA.clid = cl.id;
        DATA.clst = cl.clst;
        PlotHistogram(DATA,[]);
    end
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 22 %evaluate Gaussian means fit in the cluster space
    a = FitGaussMeans(DATA.xy,2,'verbose');
    DATA.Clusters{DATA.probe(1)}.mahal = [a.mahal a.dprime];
    GetFigure('Hist');
    subplot(2,1,2);
    hold on;
    ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif ismember(fcn,[23 26 27]) %Check spacee out with guassian mixture model.
    usegm = 0;
    ntest = 1; % set high to test that minima are reliably found
 
    if usegm
        [E, scores, dips, DATA.xy] = AutoCut(DATA, 'usegm');
    else
        if fcn == 26
            for j = 1:ntest
            [d, obj, xy, details] = BestSpace(DATA,'newtemplate');
            ds(j) = max(d);
            if j > 1
                [a,b] = max(d);
                fprintf('%d:%.2f(%d)\n',j,max(d),b)
            end
            end
            DATA = get(DATA.toplevel,'UserData');
        elseif fcn == 27
            [d, obj, xy, details] = BestSpace(DATA, 'ncells',3);
        else
            [d, obj, xy, details] = BestSpace(DATA);
        end
    if ~isfield(DATA,'TemplateScores')
        DATA = get(DATA.toplevel,'UserData');
    end
    DATA.xy = xy;
    E.bestcl = details.cid;
    E.gmfit = obj;
    [a,b] = max(d);
    if b == 4
        [c,d]= max(d(1:3));
        fprintf('GM %.2f for space %d (%.2f for %d)dt%d csd%d\n',a,b,c,d,DATA.dvdt,DATA.csd);
    else
        fprintf('GM %.2f for space %d dt%d csd%d\n',a,b,DATA.dvdt,DATA.csd);
    end

%    [dip, details] = FindDip(DATA.xy(:,1),DATA.energy(DATA.probe(1),:),'gmix');
    [dip, details] = GMDip(DATA.xy,DATA.energy(DATA.probe(1),:),'gmix');
    E.space = [6 b];
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    E.pos = [crit min(DATA.xy(:,2)) crit max(DATA.xy(:,2))];
        
    end
    DATA.gmcid = E.bestcl;
    DATA.cboundary = E;
    DATA.usegmcid = 1;
    [cl, cluster] = ClassifySpikes(DATA,E);
    PlotHistogram(DATA,E,'plotgm');
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    cluster.MeanSpike = cl.MeanSpike;
    DATA.cluster = cluster;
    DATA.cluster.gmfit = obj;
    DATA.MeanSpike = cl.MeanSpike;
    if DATA.gmtypes(b)
        DATA.plottype = DATA.gmtypes(b);
        DATA = ReplotPCs(DATA,E);
    else
        GetFigure('VarE',DATA.watcharg{:});
        subplot(1,1,1);
        PlotVarE(DATA);
        hold on;
        ezcontour(@(x,y)pdf(obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
    end
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 24
    GetFigure('Cluster',DATA.watcharg{:});
    xy = xyrotate(DATA.xy(:,1),DATA.xy(:,2),-DATA.cluster.angle);
    hold off;
    x = mean(xy);
    if DATA.plotspk.muscale < 1
    a = prctile(abs(DATA.xy(:,1) - DATA.cluster.crit),2);
    id = find(abs(DATA.xy(:,1) - DATA.cluster.crit) < a);
    r = rand(size(id)).*a;
    xid = find(abs(DATA.xy(id,1) - DATA.cluster.crit) > r);
    xid = id(xid);
    nid = setdiff(DATA.nid,xid);
    clid = setdiff(DATA.clid,xid);
    else
        nid = DATA.nid;
        clid = DATA.clid;
    end
    ms = 5;
    plot(xy(nid,1),xy(nid,2),'.','markersize',ms,'color',[0.0 0.0 0.0]);
    hold on;
    plot(xy(clid,1),xy(clid,2),'r.','markersize',ms);
%   plot(xy(xid,1),xy(xid,2),'g.','markersize',ms);
    return;
elseif fcn == 26
    PlotHistogram(DATA, DATA.cboundary,'plotgm');
    return;
elseif fcn == 28 %K means
    DATA.gmcid = kmeans(DATA.pcs(:,1:4),2,'Distance','correlation');
elseif fcn == 29 %Use ND classification
    DATA.usegmcid = ~DATA.usegmcid;
    set(a,'Checked',onoff{DATA.usegmcid+1});
elseif fcn == 30 %Classify in template space
        [a,b, DATA.xy, details] = TemplateSpace(DATA,'template','recalc');
        DATA.cluster.gmfit = b;
        DATA.cluster.gmdprime = details.gmdprime;
        fprintf('Template GM %.2f 6D %.2f 2d, %.2f1D\n',a(1),a(2),details.gmdprime)
        DATA.cluster.sign = details.cluster.sign;
        DATA.cluster.crit = details.cluster.crit;
        DATA.cluster.shape = 2;
        DATA.cluster.space = [2 4];
        DATA.cluster.angle = 0;
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.bestcl = details.bestcl;
        DATA.cluster.usegmcluster = details.usegmcluster;
else
    for j = 1:size(DATA.pcplots,1);
        yl = prctile(DATA.pcs(:,DATA.pcplots(j,2)),[100-fcn fcn]);
        xl = prctile(DATA.pcs(:,DATA.pcplots(j,1)),[100-fcn fcn]);
        GetFigure('PCs');
        subplot(2,4,j);
        set(gca,'xlim',xl,'ylim',yl);
    end
    return;
end

y = max(max(abs(DATA.AllV(DATA.probe(1),:,:))));
DATA.voffset = [1:size(DATA.AllV,1)].*y;
if ~exist('E','var')
    E = BoundaryFromCluster([],DATA.cluster);
end
if isfield(E,'bestcl')
    DATA.gmcid = E.bestcl;
end
[cl, cluster, DATA.xy] = ClassifySpikes(DATA,E);
if isfield(cluster,'r')
    DATA.distance = cluster.r;
cluster = rmfield(cluster,'r');
end
DATA.cluster = cluster;
DATA.cluster.MeanSpike = cl.MeanSpike;
if fcn == 11 
    if isfield(DATA,'distance')
    GetFigure('TemplateScores');
    plot(DATA.distance,r,'.','m=arkersize',1);
    end
end
DATA.clid = cl.id;
DATA.nid = cl.nid;
    DATA.clst = cl.clst;

DATA.MeanSpike = cl.MeanSpike;
if ismember(fcn, [1 2 6 7 8 11 30])%cluster cut called from GUI, so replot
    DATA = ReplotPCs(DATA,E);
    if isfield(cluster,'gmfit2dman') && cluster.gmfit2dman.Converged > -1
        E.gmfit2dman = cluster.gmfit2dman;
        E.mahal = cluster.mahal;
    end
    PlotHistogram(DATA, E);
end
set(DATA.toplevel,'UserData',DATA);

function DATA = CheckTemplates(DATA, C)
if ~ismember(C.space(1),[3 4 6]) || isfield(DATA,'TemplateScores');
    return;
end
DATA.TemplateLabels = TemplateLabels(DATA,0);
Scores = CalcScores(DATA,C.MeanSpike);
if size(Scores,1) > 1
    DATA.TemplateScores(:,1)= Scores(2,1,:);
    DATA.TemplateScores(:,8)= Scores(2,2,:);
end
DATA.TemplateScores(:,2)= sum(Scores(:,1,:));
DATA.TemplateScores(:,3)= Scores(1,1,:);
if size(Scores,1) > 2
    DATA.TemplateScores(:,4)= Scores(3,1,:);
end

DATA.TemplateScores(:,10)= sum(Scores(:,2,:));
DATA.TemplateScores(:,12)= 0;
DATA.tmpdips = CalculateTemplateDips(DATA);

    function Ex = LoadExpt(DATA, ei)
        Ex = [];
        ei = floor(ei); %plain .mat files are not split just because FullV files are
        smrname = regexprep(DATA.name,'lem/M([0-9]*)','$0/lemM$1');
        exfile = [smrname '.' num2str(ei) 'idx.mat'];
        if exist(exfile,'file')
            fprintf('Loading %s\n',exfile);
            load(exfile);
            id = find(Expt.Trials.Result == 1);
            for t = length(id):-1:1
                Ex.Trials(t).Start = Expt.Trials.Start(id(t));
                Ex.Trials(t).End = Expt.Trials.End(id(t));
                Ex.Trials(t).ed = Expt.Trials.ed(id(t));
                Ex.Trials(t).id = Expt.Trials.id(id(t));
                Ex.Trials(t).Trial = id(t);
            end
            Ex.Header.expname = ExptList(end).expname;
            Ex.Stimvals.ed = mean(Expt.Trials.ed(id));
        end
        


    
function res = FitGaussMeans(X,N, varargin)
    verbose = 0;
    S = [];
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'verbose',4)
            verbose = 1;
            if length(varargin) > j && isnumeric(varargin{j+1})
                j = j+1;
                verbose = j;
            end
        elseif strncmpi(varargin{j},'clusterid',9)
            j = j+1;
            id = varargin{j};
            nid = setdiff(1:size(X,1),id);
                S.mu(1,:) = mean(X(id,:),1);
                S.mu(2,:) = mean(X(nid,:),1);
                S.Sigma(:,:,1) = cov(X(id,:));
                S.Sigma(:,:,2) = cov(X(nid,:));
                S.PComponents(1) = length(id)./size(X,1);
                S.PComponents(2) = length(nid)./size(X,1);
        end
        j = j+1;
    end

    try
        if isempty(S)
            G = gmdistribution.fit(X,N,'Options',statset('MaxIter',1000));
        else
            G = gmdistribution.fit(X,N,'Options',statset('MaxIter',1000),'Start',S);
        end
    res.obj = G;
    distance = mahal(G,G.mu);
    distance = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    for j = 1:size(G.Sigma,3)
    sigmas(j) = sqrt(sum(diag(G.Sigma(:,:,j)))); %var for this dimemsion
    end
    nsd = diff(G.mu)./sigmas;
    dprime = sqrt(sum(nsd.^2));
    res.mahal = distance;
    res.dprime = dprime;
    if verbose >1
        fprintf('Distance %.2f (%.2f)\n',distance,dprime);
    end
    catch
        res.mahal = 0;
        res.dprime = 0;
        res.obj.mu = zeros(size(X,2),N);
        res.obj.Converged = -1;
    end
    
    
    
    
function d = gmdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
    distance = mahal(G,G.mu);
    d = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    
function [x,y] = GetClusterXYData(DATA, p)
    if isempty(p) %Var/E plot
        y = DATA.spkvar(DATA.probe(1),:)./DATA.energy(DATA.probe(1),:);
        x = DATA.energy(DATA.probe(1),:);
    elseif DATA.plottype == 2 %voltage pairs
        x =DATA.AllV(p(1),p(2),:);
        y = DATA.AllV(p(3),p(4),:);
    elseif DATA.plottype == 3
        x = DATA.TemplateScores(:,p(1));
        y = DATA.TemplateScores(:,p(2));
    else
        x = DATA.pcs(:,p(1));
        y = DATA.pcs(:,p(2));
    end
x = squeeze(x);
y = squeeze(y);

function Cut = PlotHistogram(DATA, E, varargin)
    
    plotdip = 0;
    checkdprimes = 0;
    plotgm = 0;
   replotpts = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plotdip',5)
            plotdip = 1;
        elseif strncmpi(varargin{j},'plotgmdetails',10)
            plotgm = 2;
        elseif strncmpi(varargin{j},'plotgm',5)
            plotgm = 1;
        end
        j = j+1;
    end

% E.shape == 2 means that DATA.xy has already been rotated, so don't rotate here.
% if classify spike hasn't been called yet, this doesn't work.  Need to
% check for this....
if isempty(E)
    E = BoundaryFromCluster([],DATA.cluster);
end
if E.shape(1) == 2
    xy = DATA.xy;
    crit = E.xyr(1);
    allx = DATA.xy(:,1);
    ally = DATA.xy(:,2);
    exy = [E.pos(1) E.pos(2); E.pos(3) E.pos(4)];
    angle = 0;
else
    pos = E.pos;
angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));

[allx, ally] = GetClusterXYData(DATA,E.pcplot);
xy = xyrotate(allx,ally,angle);
exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);
crit = mean(exy(:,1));
end
if ~isfield(E,'sign')
    E.sign = 0;
end
if E.sign > 0
clid = find(xy(:,1) > crit);
nid = find(xy(:,1) <= crit);
else
clid = find(xy(:,1) < crit);
nid = find(xy(:,1) >= crit);
end
if replotpts
hold off;
plot(allx,ally,'.','markersize',1);
hold on;
plot(E.pos([1 3]),E.pos([2 4]),'r-');
plot(mean(allx(clid)),mean(ally(clid)),'r+');
plot(mean(allx(nid)),mean(ally(nid)),'r+');
end
cfig = SetFigure('Hist',DATA.watcharg{:});
subplot(2,1,2);
xid = [];
showgmfit = 0;
if isfield(E,'bestcl') & length(E.bestcl) == length(xy)
    fprintf('XY using GM clustering\n');
    nid = find(E.bestcl == 1);
    clid = find(E.bestcl == 2);
    xid=find(E.bestcl ==3);
    diffid = find(E.bestcl ==1 & xy(:,1) .* E.sign > crit .* E.sign);
    diffcid = find(E.bestcl ==2 & xy(:,1) .* E.sign < crit .* E.sign);
elseif isfield(E,'gmfit2dman')
    showgmfit = 1;
    fprintf('XY using GM 2D clustering\n');
if E.shape == 0
    E.bestcl = cluster(E.gmfit2dman,DATA.xy);
elseif E.shape == 1
    E.bestcl = cluster(E.gmfit2dman,xy);
else
    E.bestcl = cluster(E.gmfit2dman,cat(2,allx,ally));
end
    nid = find(E.bestcl == 1);
    clid = find(E.bestcl == 2);
    xid=find(E.bestcl ==3);
    diffid = find(E.bestcl ==1 & xy(:,1) .* E.sign > crit .* E.sign);
    diffcid = find(E.bestcl ==2 & xy(:,1) .* E.sign < crit .* E.sign);
else
    diffcid = [];
end
hold off;
plot(xy(nid,1),xy(nid,2),'.','markersize',1);
hold on;
plot(xy(clid,1),xy(clid,2),'r.','markersize',1);
plot(xy(xid,1),xy(xid,2),'g.','markersize',1);
if showgmfit
    plot(E.gmfit2dman.mu(1,1),E.gmfit2dman.mu(1,2),'c+','linewidth',2);
    plot(E.gmfit2dman.mu(2,1),E.gmfit2dman.mu(2,2),'c+','linewidth',2);
end
if length(diffcid) && DATA.watchplots
PlotSpikes(DATA,cat(1,diffid, diffcid),'fixy');
figure(cfig);
end
if plotgm
   a = FitGaussMeans(xy,2);
   ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
   if isfield(E,'bestspace')
       str = sprintf(' Mahal %.2f (%.2f)',a.mahal,E.bestspace(1));
   else
       str = sprintf(' Mahal %.2f (%.2f)',a.mahal);
   end
elseif isfield(E,'bestspace') 
    str = sprintf('Mahal %.2f',E.bestspace(1));
elseif isfield(E,'mahal')
    str = sprintf('Mahal %.2f.%.2f',E.mahal(1),E.mahal(2));
else
    str = [];
end

if E.shape == 2
    title(sprintf('ND: %s%s',DATA.gmtypelabels{E.space(2)},str));
elseif isempty(E.pcplot)
    title(sprintf('Var-E %s',str));
else
    title(sprintf('%d: %dvs%d%s',DATA.plottype,E.pcplot(1),E.pcplot(2),str));
end



hdat = get(gcf,'UserData');
hdat.elmousept.pos(1) = crit;
hdat.elmousept.pos(3) = crit;
hdat.elmousept.pos(2) = exy(1,2);
hdat.elmousept.pos(4) = exy(2,2);
hdat.elmousept.shape = 1;
hdat.elmousept.down = 0;
hdat.elmousept.h = plot(exy(:,1),exy(:,2),'r-');

dp = (mean(xy(clid,1))-mean(xy(nid,1)))./sqrt(mean([var(xy(clid,1)) var(xy(nid,1))]));
hold off;
subplot(2,1,1);
hold off; 
[a,b] = hist(xy(:,1),500);
bar(b,a,1);
hold on;
plot(exy(:,1),get(gca,'ylim'),'r-');
%[dip, mydip] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:),'eval',crit,'plot','gmix');
if plotgm == 2
    [dip, mydip] = GMDip(xy, DATA.energy(DATA.probe(1),:),'plot','crit',crit);
    if mydip.converged(1) == 0
        FindDip(xy(:,1),DATA.energy(DATA.probe(1),:),'eval',crit,'plot','gmix');
    end
elseif isfield(E,'gmfit1d')
    [dip, mydip] = GMDip(xy, E.gmfit1d);
    mydip.sign = E.sign;
else
    [dip, mydip] = GMDip(xy, DATA.energy(DATA.probe(1),:));
end
if E.sign == 0
    E.sign  = mydip.sign;
end
mycrit = dip;
if isfield(mydip,'gxy') %pdf of GM fit in 1D
    scale = trapz(b,a);
    plot(mydip.gxy(:,1),sum(mydip.gxy(:,[2 3]),2).*scale,'r');
    plot(mydip.gxy(:,1),mydip.gxy(:,2).*scale,'g');
    plot(mydip.gxy(:,1),mydip.gxy(:,3).*scale,'g');
end
if mydip.type == 1
    plot([dip(1) dip(1)],get(gca,'ylim'),'g-');
    plot([dip(2) dip(2)],get(gca,'ylim'),'m-');
    plot([dip(3) dip(3)],get(gca,'ylim'),'m--');
    if length(dip) > 3
    plot([dip(4) dip(4)],get(gca,'ylim'),'c-');
    end
else
plot([dip(1) dip(1)],get(gca,'ylim'),'g-');
plot([dip(2) dip(2)],get(gca,'ylim'),'g-');
end
if length(dip) >4 && ~isnan(dip(5))
plot([dip(5) dip(5)],get(gca,'ylim'),'m-');
plot([dip(6) dip(6)],get(gca,'ylim'),'m--');
end    
if checkdprimes
[a,b] = MaxDprime(xy(:,1));
yl = get(gca,'ylim');
scale = yl(2)./10;
plot(b.crit,abs(b.dps).*scale,'r');
end
dip = HartigansDipTest(sort(xy(:,1))).*100;
bii = BimodalCoeff(xy(:,1),1.5);
title(sprintf('P%d Dip %.1f(%.1f,%.2f gm%.2f)',DATA.probe(1),dip,mydip.dipsize(1),bii,mydip.gmdprime));
DATA.clid = clid;
DATA.nid = nid;
Cut = E;
Cut.dip = mycrit;
Cut.angle = angle;
Cut.crit = [crit mycrit];
Cut.hdip = dip;
Cut.mydip = [mydip.cdipsize mydip.dipsize];
Cut.space = [DATA.plottype E.pcplot];
Cut.shape = E.shape;
Cut.y = exy(:,2);
DATA.cluster.shape = 1;
hdat.cluster = Cut;
C = ClusterFromBoundary(E, Cut);
newE = [];
newE = BoundaryFromCluster(newE, C);

set(gcf,'UserData',hdat);
%DATA = ReplotPCs(DATA,E);
%DrawLine(newE);

function F = SetFigure(lb, varargin)

[F, isnew] = GetFigure(lb,varargin{:});
if isnew 
    if strcmp(lb,'PCs')
        hm = uimenu(F,'Label','Cluster','Tag','ClusterMenu');
        uimenu(hm,'Label','Draw &Ellipse','Callback',{@PCCluster, 1});
        uimenu(hm,'Label','Draw &Ellipse 2','Callback',{@PCCluster, 32});
        uimenu(hm,'Label','VarE Ellipse','Callback',{@PCCluster, 2});
        uimenu(hm,'Label','Density','Callback',{@PCCluster, 4});
        uimenu(hm,'Label','ReDraw','Callback',{@PCCluster, 3});
        uimenu(hm,'Label','Template','Callback',{@PCCluster, 5});
        uimenu(hm,'Label','StdTemplate','Callback',{@PCCluster, 32});
        uimenu(hm,'Label','TemplateLine','Callback',{@PCCluster, 6});
        uimenu(hm,'Label','Var E Line','Callback',{@PCCluster, 7});
        uimenu(hm,'Label','PC &Line','Callback',{@PCCluster, 8});
        uimenu(hm,'Label','PC Line2','Callback',{@PCCluster, 11});
        uimenu(hm,'Label','Optimize Line','Callback',{@PCCluster, 18});
        uimenu(hm,'Label','AutoMatic Cut','Callback',{@PCCluster, 19});
        uimenu(hm,'Label','GaussMeans','Callback',{@PCCluster, 22});
        uimenu(hm,'Label','BestSpace','Callback',{@PCCluster, 23});
        uimenu(hm,'Label','BestSpace (scratch)','Callback',{@PCCluster, 26});
        uimenu(hm,'Label','BestSpace (2 cells)','Callback',{@PCCluster, 27});
        uimenu(hm,'Label','K means','Callback',{@PCCluster, 28});
        uimenu(hm,'Label','Use ND clid','Callback',{@PCCluster, 29});
        uimenu(hm,'Label','TemplateSpace','Callback',{@PCCluster, 30});
        uimenu(hm,'Label','PCSpace','Callback',{@PCCluster, 31});
        
        sm = uimenu(hm,'Label','Plot','Tag','ClusterMenu');
        uimenu(sm,'Label','&PCS','Callback',{@PCCluster, 9});
        uimenu(sm,'Label','&ADC','Callback',{@PCCluster, 10});
        uimenu(sm,'Label','&Templates','Callback',{@PCCluster, 14});
        uimenu(sm,'Label','Templates2','Callback',{@PCCluster, 15});
        uimenu(sm,'Label','&dvdt','Callback',{@PCCluster, 20});
        uimenu(sm,'Label','&Cluster','Callback',{@PCCluster, 24});
        uimenu(sm,'Label','&Histogram','Callback',{@PCCluster, 26});
        sm = uimenu(hm,'Label','Axes:tight','Tag','ClusterMenu');
        uimenu(sm,'Label','100%','Callback',{@PCCluster, 100});
        uimenu(sm,'Label','99%','Callback',{@PCCluster, 99});
        sm = uimenu(hm,'Label','PlotClusters','Tag','PlotClusters');
        uimenu(sm,'Label','SpkTimes','Callback',{@PlotClusters, 1});
        uimenu(sm,'Label','xCorr','Callback',{@PlotClusters, 2});
        uimenu(sm,'Label','Quality','Callback',{@PlotClusters, 3});
        uimenu(hm,'Label','Delete','Callback',{@PCCluster, 17});
        uimenu(hm,'Label','Revert','Callback',{@PCCluster, 21});
        uimenu(hm,'Label','&Save','Callback',{@PCCluster, 12});
        uimenu(hm,'Label','&Save Spikes','Callback',{@PCCluster, 25});
        uimenu(hm,'Label','&Xcorr','Callback',{@CalcXcorr, 25});
        hm = uimenu(F,'Label','Probes','Tag','ProbeMenu');
        uimenu(hm,'Label','1','Callback',{@ProbeMenu, 1});
        uimenu(hm,'Label','3','Callback',{@ProbeMenu, 2});
        uimenu(hm,'Label','5','Callback',{@ProbeMenu, 3});
        uimenu(hm,'Label','dvdt','Callback',{@ProbeMenu, 4});
        uimenu(hm,'Label','csd','Callback',{@ProbeMenu, 5});
        uimenu(hm,'Label','Next','Tag','NextButton','Callback',{@PCCluster, 13});
        set(F, 'KeyPressFcn',@PCKeyPressed);
    elseif strcmp(lb,'FullV')
%            set(F, 'WindowScrollWheelFcn',@ScrollV);
    elseif strcmp(lb,'Hist')
        DATA.parentfigtag = 'PCs';
        set(F,'UserData',DATA);
        set(F, 'KeyPressFcn',@HistKeyPressed);
        set(F, 'WindowButtonDownFcn',@HistButtonPressed);
        set(F, 'WindowButtonMotionFcn',@HistButtonDragged);
        set(F, 'WindowButtonUpFcn',@HistButtonReleased);
        hm = uimenu(F,'Label','Plot','Tag','ClusterMenu');
        uimenu(hm,'Label','Dip Criterion','Callback',{@HistMenu, 1});
        uimenu(hm,'Label','mahal/angle','Callback',{@HistMenu, 2});
        uimenu(hm,'Label','Hist Replot','Callback',{@HistMenu, 3});
    elseif strcmp(lb,'Clusters')
        DATA.parentfigtag = 'PCs';
        hm = uimenu(F,'Label','Plot','Tag','ClusterMenu');
        set(F,'UserData',DATA);
        uimenu(hm,'Label','Times','Callback',{@PlotCluster, 1});
        uimenu(hm,'Label','xcorr','Callback',{@PlotCluster, 2});
        uimenu(hm,'Label','xcorr adjacent','Callback',{@PlotCluster, 3});
        uimenu(hm,'Label','xcorr 2sep','Callback',{@PlotCluster, 4});
        uimenu(hm,'Label','Quality Scatter','Callback',{@PlotCluster, 5});
        uimenu(hm,'Label','QualityProbes','Callback',{@PlotCluster, 6});
    elseif strcmp(lb,'Spikes')
        DATA.parentfigtag = 'PCs';
        set(F,'UserData',DATA);
        set(F, 'WindowScrollWheelFcn',@ScrollSpikes);
        set(F, 'KeyPressFcn',@KeyPressed);

        hm = uimenu(F,'Label','Scroll','Tag','ClusterMenu');
        uimenu(hm,'Label','More Spikes','Callback',{@SpikeDraw, 1});
        uimenu(hm,'Label','Fewer Spikes','Callback',{@SpikeDraw, 2});
        uimenu(hm,'Label','Spool Spikes','Callback',{@SpikeDraw, 3});
        uimenu(hm,'Label','Spool by Trial','Callback',{@SpikeDraw, 10});
        uimenu(hm,'Label','dVdt','Callback',{@SpikeDraw, 4});
        uimenu(hm,'Label','csd','Callback',{@SpikeDraw, 5});
        uimenu(hm,'Label','Subtact mean','Callback',{@SpikeDraw, 6});
        uimenu(hm,'Label','Subtract peak','Callback',{@SpikeDraw, 7});
        uimenu(hm,'Label','Subtract Min','Callback',{@SpikeDraw, 8});
        uimenu(hm,'Label','Main Probe Only','Callback',{@SpikeDraw, 9});        
        uimenu(hm,'Label','Exclude Later Spikes','Callback',{@SpikeDraw, 11});        
        uimenu(hm,'Label','Exclude Earlier Spikes','Callback',{@SpikeDraw, 12});        
        bp = [0.95 0.95 0.1 0.05];
        uicontrol(F,'style','pop','string','1|2|3','Units','Normalized','Position',bp,'Tag','ChooseTrial',...
            'Callback', @SelectTrial);
    elseif strcmp(lb,'TemplateScores')
        DATA.parentfigtag = 'PCs';
        set(F,'UserData',DATA);
        hm = uimenu(F,'Label','Scroll','Tag','Classify');
        uimenu(hm,'Label','Line','Callback',{@PCCluster, 6});
    end
end


function CalcXcorr(a,b,fcn)
DATA = GetDataFromFig(src);
ta = DATA.t(DATA.clst ==1);        
tb = DATA.t(DATA.clst ==2);
xc = xcorrtimes(ta,tb);

function PCKeyPressed(src, ks)

DATA = GetDataFromFig(src);
spk = [0 0];
changed = 1;


if strmatch(ks.Key,'p')
    PCCluster(DATA,[],9); %plot PCS
elseif strmatch(ks.Key,'t') 
    PCCluster(DATA,[],14); %plot Templates
elseif strmatch(ks.Key,'a') 
    PCCluster(DATA,[],10); %plot ADCs
elseif strmatch(ks.Key,'d') 
    PCCluster(DATA,[],20); %plot ADCs
elseif strmatch(ks.Key,'s') 
    PCCluster(DATA,[],12); %save cluster
elseif strmatch(ks.Key,'l') 
    PCCluster(DATA,[],8); %cut with line
elseif strmatch(ks.Key,'e') 
    PCCluster(DATA,[],1); %ellipse
elseif strmatch(ks.Key,'n')
    if strmatch('shift',ks.Modifier)
    PCCluster(DATA,[],13); %Next Probe
    end
elseif strmatch(ks.Key,'v') 
    SetFigure('Spikes',DATA.watcharg{:});
    SpoolSpikes(DATA); %ellipse
    figure(DATA.toplevel);
else
    changed = 0;
end

function SelectTrial(src, b)
DATA = GetDataFromFig(src);
id = get(src,'value');
PlotTrialSpikes(DATA,id); %Sets UserDATA

function SetTrialList(DATA)
    it = findobj('Tag','ChooseTrial');
    if ~isempty(it) & isfield(DATA.Expt,'Trials')
    set(it,'string',sprintf('%d|',[DATA.Expt.Trials.Trial]),'value',1);
    end
        
        

function HistKeyPressed(src, ks)

DATA = GetDataFromFig(src);
spk = [0 0];
changed = 1;


if strmatch(ks.Key,'rightarrow')
    if DATA.cboundary.shape == 2
        DATA.xy = xyrotate(DATA.xy(:,1),DATA.xy(:,2),5 * pi/180);
    else
        DATA.cboundary.pos = RotateLine(DATA.cboundary.pos,5 * pi/180);
        if isfield(DATA.cboundary,'axis')
            axes(DATA.cboundary.axis);
        end
    end
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plotgmdetails');
    GetFigure('Hist',DATA.watcharg{:});
elseif strmatch(ks.Key,'leftarrow')
    if DATA.cboundary.shape == 2
        DATA.xy = xyrotate(DATA.xy(:,1),DATA.xy(:,2),-5 * pi/180);
    else
        DATA.cboundary.pos = RotateLine(DATA.cboundary.pos,-5 * pi/180);
        if isfield(DATA.cboundary,'axis')
            axis(DATA.cboundary.axis);
        end
    end
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plotgmdetails');
    GetFigure('Hist',DATA.watcharg{:});
elseif strmatch(ks.Key,'d')
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plot');
else
    changed = 0;
end
if changed
set(DATA.toplevel,'UserData',DATA);
end

    function pos = RotateLine(pos, da)
    oldangle =  atan(diff(pos([1 3]))/diff(pos([2 4])));
    angle = oldangle+da;
    len = abs(diff(pos([1 3])) + i * diff(pos([2 4])))/2;
    x = mean(pos([1 3]));
    y = mean(pos([2 4]));
    pos(1) = x + len* sin(angle);
    pos(3) = x - len* sin(angle);
    pos(2) = y + len* cos(angle);
    pos(4) = y - len* cos(angle);

 function KeyPressed(src, ks)

DATA = GetDataFromFig(src);
spk = [0 0];


if isempty(DATA.spklst)
    DATA.spklst = 1:100;
end
w = DATA.spksperview;


if strmatch(ks.Key,'delete') 
    DeleteCluster(mousept.cluster, a);
    mousept.mode = 0;
    mousept.angle = 0;
    if mousept.lasth & ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
elseif strmatch(ks.Key,'rightarrow')
    DATA.stopspool = 1;
    w = DATA.spksperview;
    spk = minmax(DATA.spklst);
    nt = DATA.currenttrial+1;
    spk(1) = spk(2);
    spk(2) = spk(2)+w;
elseif strmatch(ks.Key,'leftarrow')
    DATA.stopspool = 1;
    w = DATA.spksperview;
    spk = minmax(DATA.spklst);
    nt = DATA.currenttrial-1;
    spk(2) = spk(1);
    spk(1) = spk(1)-w;
elseif strmatch(ks.Key,'downarrow')
    DATA.stopspool = 1;
    DATA.spklst = DATA.spklst+1;
    PlotSpikes(DATA,DATA.spklst(1),'fixy');
    spk = [0 0];
elseif strmatch(ks.Key,'uparrow')
    DATA.stopspool = 1;
    DATA.spklst = DATA.spklst-1;
    PlotSpikes(DATA,DATA.spklst(1),'fixy');
    spk = [0 0];
elseif strmatch(ks.Key,'add')
    it = findobj(a,'Tag','Clusterid');
    c = get(it,'value');
    set(it,'value',c+1);
elseif strmatch(ks.Key,'subtract')
    it = findobj(a,'Tag','Clusterid');
    c = get(it,'value');
    if c > 1 
    set(it,'value',c-1);
    end
elseif strmatch(ks.Key,'space')
    DATA.stopspool = 1;
elseif ks.Key == 'r'
    mousept.angle = mousept.angle+0.02;
    mousept.mode = 10;
    if mousept.lasth & ishandle(mousept.lasth)
    delete(mousept.lasth);
    end
   mousept = myellipse(mousept,[0 0; 0 0 ]);
end
spklst = [];
if DATA.plotspk.bytrial
    PlotTrialSpikes(DATA,nt);
    return;
end
if spk(2) > 0 && spk(1) < size(DATA.AllV,3)
    spklst = spk(1):spk(2);
    DATA.spklst = spklst;
    PlotSpikes(DATA,DATA.spklst,'fixy');
    set(DATA.toplevel,'UserData',DATA);
end
function PlotTrialSpikes(DATA, nt)
    
    if nt > 0 && nt <= length(DATA.Expt.Trials)
        spklst = find(DATA.t .*10000 > DATA.Expt.Trials(nt).Start & DATA.t.*10000 < DATA.Expt.Trials(nt).End);
       DATA.currenttrial = nt;
    DATA.spklst = spklst;
    PlotSpikes(DATA,DATA.spklst,'fixy');
    set(DATA.toplevel,'UserData',DATA);
    end


function PlotCluster(a,b, mode)
DATA = GetDataFromFig(a);

if mode == 1
    PlotSpikeTimes(DATA.Clusters);
elseif mode == 2
    PlotSpikeTimes(DATA.Clusters,'xcorr');
elseif mode == 3
    PlotSpikeTimes(DATA.Clusters,'xcorrstep',1);
elseif mode == 4
    PlotSpikeTimes(DATA.Clusters,'xcorrstep',2);
elseif mode == 5
    PlotSpikeTimes(DATA.Clusters,'dips');
elseif mode == 6
    PlotSpikeTimes(DATA.Clusters,'probequality');
end

function HistMenu(a,b, mode)
DATA = GetDataFromFig(a);
if mode == 1
%    dip = FindDip(DATA.xy(:,1),DATA.energy(DATA.probe(1),:),'plot');
    dip = GMDip(DATA.xy,DATA.energy(DATA.probe(1),:),'plot');
elseif mode == 2
    [a,b,c] = BestGMAngle(DATA.xy(:,1),DATA.xy(:,2));
    GetFigure('dips',DATA.watcharg{:});
    hold off; 
    plot(c.angles,c.mahal);
    title(sprintf('Mahal 1D %.3f at %.2fdeg , 2D %.3f',b, a .*180/pi,c.mahal2d));
elseif mode == 3
    PlotHistogram(DATA, DATA.cboundary,'plotgmdetails');
end
    
function SpikeDraw(a,b, mode)
DATA = GetDataFromFig(a);
onoff = {'off' 'on'};
if mode == 2
    DATA.spksperview = round(DATA.spksperview/2);
elseif mode == 1
    DATA.spksperview = DATA.spksperview *2;
elseif mode == 3
   SpoolSpikes(DATA);
elseif mode == 4
    DATA.plotdvdt = ~DATA.plotdvdt;
    set(a,'Checked',onoff{DATA.plotdvdt+1});
    PlotMeanSpike(DATA);
elseif mode == 5
    DATA.plotcsd = ~DATA.plotcsd;
    set(a,'Checked',onoff{DATA.plotcsd+1});
    PlotMeanSpike(DATA);
elseif mode == 6
    DATA.plotspk.submean = ~DATA.plotspk.submean;
    set(a,'Checked',onoff{DATA.plotspk.submean+1});
    PlotMeanSpike(DATA);
elseif mode == 7
    DATA.plotspk.submax = ~DATA.plotspk.submax;
    set(a,'Checked',onoff{DATA.plotspk.submax+1});
    PlotMeanSpike(DATA);
elseif mode == 8
    DATA.plotspk.submin = ~DATA.plotspk.submin;
    set(a,'Checked',onoff{DATA.plotspk.submim+1});
    PlotMeanSpike(DATA);
elseif mode == 9
    DATA.plotspk.oneprobe = ~DATA.plotspk.oneprobe;
    set(a,'Checked',onoff{DATA.plotspk.oneprobe+1});
    PlotMeanSpike(DATA);
elseif mode == 10
    DATA.plotspk.bytrial = ~DATA.plotspk.bytrial;
    set(a,'Checked',onoff{DATA.plotspk.bytrial+1});
    SpoolSpikes(DATA);
elseif mode == 11  %restrict time range
    if ~isempty(DATA.restricttimerange)
        t(1) = DATA.restricttimerange(1);
    else
        t(1) = DATA.t(1)-0.01;
    end
    t(2) = DATA.t(DATA.spklst(end))+0.01;
    DATA = RestrictTimeRange(DATA,t);
elseif mode == 12 %restrict time range
    if ~isempty(DATA.restricttimerange)
        t(2)= DATA.restricttimerange(2);
    else
        t(2) = DATA.t(end)+0.01;
    end
    t(1)= DATA.t(DATA.spklst(1))-0.01;
    DATA = RestrictTimeRange(DATA,t);
end
if ismember(mode,[1 2 4 5 6 7 8 9])
    DATA.spklst = DATA.spklst(1):DATA.spklst(1)+DATA.spksperview-1;
    PlotSpikes(DATA,DATA.spklst);
    set(DATA.toplevel,'UserData',DATA);
end

function SpoolAllSpikes(DATA, varargin)
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'front',3)
            SetFigure('Spikes','front');
        end
            
        j = j+1;
    end
SetTrialList(DATA);
set(gca,'Position',[0.05 0.05 0.9 0.9]);
        hold off;
npts = 32;
for p = 1:size(DATA.V,1)
    DATA.handles(p,1) = plot([-8 32],[0 0], 'color',[0.5 0.5 0.5]);
    hold on;
    DATA.handles(p,2) = plot([-8 32],[0 0], 'color','r');
    xoff = floor((p-1)/6) .* npts;
    yoff = (rem(p-1,6)) * 3;
    text(xoff,yoff+1,sprintf('%d:%.1f',p,DATA.Clusters{p}.mahal(1)));
end
set(gca,'xlim',[-8 120],'ylim',[-3 16]); 
th =  title('Trial');
ts = now;
for nt = 1:length(DATA.Expt.Trials)
    T = DATA.Expt.Trials(nt);
    for p = 1:size(DATA.V,1)
        spklist = find(DATA.vt(DATA.ids{p}) .*10000 > DATA.Expt.Trials(nt).Start & DATA.vt(DATA.ids{p}).*10000 < DATA.Expt.Trials(nt).End);
        xoff = floor((p-1)/6) .* npts;
        yoff = (rem(p-1,6)) * 3;
        PlotProbeSpikes(DATA, p, DATA.ids{p}(spklist),[npts 8],[xoff yoff]);
    end
        drawnow;
end
fprintf('Took %.2f\n',mytoc(ts));
    DATA.currenttrial = nt;

function PlotProbeSpikes(DATA, p, spklist,npts,offset)
    j = 1;
    xoff = offset(1);
    yoff = offset(2);
  hold on;
  if isempty(spklist)
      set(DATA.handles(p,1),'Ydata',[0 0]+yoff,'Xdata',[npts(2) npts(1)+npts(2)]+xoff);
      return;
  end
  
  x = [1:npts(1)]-npts(2);
  X = repmat(x,length(spklist),1)';
  id = repmat(spklist,npts(1),1) + X;
  nV = DATA.V(p,id);
  nV(npts(1)+1:npts(1):end) = NaN;
  X(:,end) = NaN;
  set(DATA.handles(p,1),'Ydata',nV(:)+yoff,'Xdata',X(:)+xoff);
  tid = find(ismember(DATA.vt(spklist),DATA.Clusters{p}.times));
  if length(tid)
      X = repmat(x,length(tid),1)';
      id = repmat(spklist(tid),npts(1),1) + X;
      nV = DATA.V(p,id);
      nV(npts(1)+1:npts(1):end) = NaN;
      X(:,end) = NaN;
      set(DATA.handles(p,2),'Ydata',nV(:)+yoff,'Xdata',X(:)+xoff);
  end

function SpoolSpikes(DATA, varargin)
    
    SpoolAllSpikes(DATA,varargin);
    return;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'front',3)
            SetFigure('Spikes','front');
        end
            
        j = j+1;
    end
    
    if ~isfield(DATA,'plotcsd');
        DATA.plotcsd = 0;
    end
DATA.stopspool = 0;
SetTrialList(DATA);

spklist = 1:DATA.spksperview;
for j = DATA.plotspk.probes
if DATA.plotcsd
maxV(j,:,:) = max(max(DATA.AllCSD(j,:,:),[],3),[],2);
minV(j,:,:) = min(min(DATA.AllCSD(j,:,:),[],3),[],2);
else
maxV(j,:,:) = max(max(DATA.AllV(j,:,:),[],3),[],2);
minV(j,:,:) = min(min(DATA.AllV(j,:,:),[],3),[],2);
end
end
DATA.voffset = cumsum(max(abs([maxV minV]),[],2));
DATA.spkyrange(1) = min(DATA.voffset(DATA.plotspk.probes)-DATA.voffset(DATA.probe(1))) + minV(DATA.probe(1));
DATA.spkyrange(2) = max(DATA.voffset(DATA.plotspk.probes)-DATA.voffset(DATA.probe(1))) + maxV(DATA.probe(1));
set(DATA.toplevel,'UserData',DATA);

go = 1;
nt = 0;
if DATA.plotspk.bytrial
    id = find([DATA.Expt.Trials.Start] > DATA.t(DATA.uid(1)).*10000);
    if length(id)
        nt = id(1)-1;
    else
        nt = 0;
    end
end
if isempty(DATA.restricttimerange)
    maxtime = NaN;
    mintime = NaN;
else
    maxtime = DATA.restricttimerange(2) * 10000;
    mintime = DATA.restricttimerange(1) * 10000;
end
while (isempty(spklist) || spklist(end) < size(DATA.AllV,3)) && go
    if DATA.plotspk.bytrial
        nt = nt+1;
        spklist = find(DATA.t .*10000 > DATA.Expt.Trials(nt).Start & DATA.t.*10000 < DATA.Expt.Trials(nt).End);
        if nt >= length(DATA.Expt.Trials) || DATA.Expt.Trials(nt).Start > maxtime
            go = 0;
        elseif DATA.Expt.Trials(nt).Start < mintime
            go = 2;
        else
            go = 1;
        end
    else
        spklist = spklist + DATA.spksperview;
    end
    if go == 1
    DATA.currenttrial = nt;
    PlotSpikes(DATA,spklist,'fix',DATA.spkyrange);
    drawnow;
    DATA = get(DATA.toplevel,'UserDATA');
    end
    if DATA.stopspool
        go = 0;
        DATA.spklst = spklist;
    end
end
DATA.currenttrial = nt;
DATA.stopspool = 0;
set(DATA.toplevel,'UserData',DATA);

function ScrollV(src, evnt)
DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

DATA = GetDataFromFig(src);
xl = get(gca,'xlim');
w = diff(xl);
if sign(evnt.VerticalScrollCount) > 0
    xl(1) = xl(2);
    xl(2) = xl(2)+w;
    set(gca,'xlim',xl);
elseif sign(evnt.VerticalScrollCount) < 0
    xl(2) = xl(1);
    xl(1) = xl(1)-w;
    set(gca,'xlim',xl);
end

function ScrollSpikes(src, evnt)
DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

DATA = GetDataFromFig(src);
if isempty(DATA.spklst)
    DATA.spklst = 1:100;
end
if DATA.plotspk.bytrial
    nt = DATA.currenttrial + evnt.VerticalScrollCount;
    PlotTrialSpikes(DATA,nt);
    return;
else
spk = minmax(DATA.spklst);
w = DATA.spksperview;
if sign(evnt.VerticalScrollCount) > 0
    spk(1) = spk(2);
    spk(2) = spk(2)+w;
elseif sign(evnt.VerticalScrollCount) < 0 
    spk(2) = spk(1);
    spk(1) = spk(1)-w;
end
end

if spk(2) > 0 && spk(1) < size(DATA.AllV,3)
DATA.spklst = spk(1):spk(2);
PlotSpikes(DATA,DATA.spklst);
end
set(DATA.toplevel,'UserData',DATA);
    
        

function PlotPCs(pcs, a,b, type, id, colors, varargin)
ptsz = 1;
fixrange  = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fixrange',5)
        fixrange = 1;
    end
    j = j+1;
end
if type == 0
    hold off;
    if fixrange
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
    end
    
    if length(id) == size(pcs,1) 
        nc = unique(id);
        nc = nc(nc > 0);
        if length(nc) < 8
        for j = 1:length(nc)
            nid = find(id == nc(j));
            plot(pcs(nid,a),pcs(nid,b),'.','markersize',ptsz,'color',colors{nc(j)});
            hold on;
        end
        else
            plot(pcs(id,a),pcs(id,b),'r.','markersize',ptsz);
        end
    elseif length(id)
        nid = setdiff(1:size(pcs,1),id);
        if length(nid)
        plot(pcs(nid,a),pcs(nid,b),'.','markersize',ptsz);
        hold on;
        end
        plot(pcs(id,a),pcs(id,b),'r.','markersize',ptsz);
        hold on;
    else
        plot(pcs(:,a),pcs(:,b),'.','markersize',ptsz);
    end
    if fixrange
        set(gca,'xlim',xl);
        set(gca,'ylim',yl);
    end
else
    DensityPlot(pcs(:,a),pcs(:,b),'sd',[2 2],'ynormal');
end


function PlotVals(DATA, a,b, type, id, colors, varargin)
ptsz = 1;
fixrange  = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fixrange',5)
        fixrange = 1;
    end
    j = j+1;
end
if type == 0
    hold off;
    if fixrange
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
    end
    if length(id) == size(DATA.AllV,3)
        nc = unique(id);
        nc = nc(nc > 0);
        for j = 1:length(nc)
            ids{j} = find(id == nc(j));
        end
    elseif length(id)
        ids{2} = id;
        ids{1} =setdiff(1:size(DATA.AllV,3),id);
    else 
        ids{1} = 1:size(DATA.AllV,3);
    end
    if DATA.plotcsd
        for j = 1:length(ids)
            plot(squeeze(DATA.AllCSD(a(1),a(2),ids{j})),squeeze(DATA.AllCSD(b(1),b(2),ids{j})),'.','markersize',ptsz,'color',colors{nc(j)});
            hold on;
        end
    elseif DATA.plottype == 6
        for j = 1:length(ids)
            plot(squeeze(DATA.dV(a(1),a(2),ids{j})),squeeze(DATA.dV(b(1),b(2),ids{j})),'.','markersize',ptsz,'color',colors{j});
            hold on;
        end        
    else
        for j = 1:length(ids)
        plot(squeeze(DATA.AllV(a(1),a(2),ids{j})),squeeze(DATA.AllV(b(1),b(2),ids{j})),'.','markersize',ptsz,'color',colors{j});
        hold on;
        end
    end
    title(sprintf('%d:%d vs %d:%d',a(1),a(2),b(1),b(2)));
    if fixrange
        set(gca,'xlim',xl);
        set(gca,'ylim',yl);
    end
else
    DensityPlot(DATA.AllV(a(1),a(2),:),DATA.AllV(a(1),a(2),:),'sd',[2 2]);
end

function [cl, cluster, xy] = ClassifySpikes(DATA, E)

    if ~isfield(E,'pcplot') && isfield(E,'space')
        cluster = E;
        DATA.plottype = cluster.space(1);
        E = BoundaryFromCluster([],cluster);
    end
    p = E.pcplot;
    cx = E.xyr(1);
    cy = E.xyr(2);
    if length(E.xyr) > 3
    rx = E.xyr(3);
    ry = E.xyr(4);
    end
    cluster.xyr = E.xyr;
    ispk = DATA.probe(1);
    cluster.nspks = size(DATA.AllV,3);
    cluster.shape = E.shape;
    cluster.minenergy = DATA.minenergy;
    cluster.minvar = DATA.minvar;
    cluster.Trigger = DATA.Trigger;
    cluster.spts = DATA.spts;
    cluster.dvdt = DATA.dvdt;
    cluster.csd = DATA.csd;
    cluster.ctime = now;
    cluster.eveci = DATA.Evec.Eval(1)./sum(DATA.Evec.Eval);
    cluster.pcgms = DATA.dipvals;
    cluster.probe = DATA.probe(1);
    if isfield(DATA.Expt,'exptno')
        cluster.exptno = DATA.Expt.exptno;
    end
    f = {'bestspace' 'bestd' 'auto' 'autotook' 'firstspace' 'firstbmi' 'bestll' 'gmdip' 'gmdipres' 'gmfit1d' 'gmfit'};
    for j = 1:length(f)
        if isfield(E,f{j})
            cluster.(f{j}) = E.(f{j});
        end
    end
    if ~isempty(DATA.lastcut)
        cluster.first = DATA.lastcut;
    end
    if isfield(DATA,'clst')
        cl.clst = DATA.clst;
    end
    if isfield(E,'cluster')
        cluster.cluster = E.cluster;
    else
        cluster.cluster = 1;
    end

    if isfield(E,'usegmcluster') && E.usegmcluster == 1
        xy = DATA.xy;
        id = find(E.bestcl == 1);
        nid = find(E.bestcl == 2);
        e(1) = mean(DATA.energy(DATA.probe(1),id));
        e(2) = mean(DATA.energy(DATA.probe(1),nid));
        if diff(e) > 1
            id = find(E.bestcl == 2);
            nid = find(E.bestcl == 1);
        end
        r = xy(:,1);
    elseif E.shape == 2 %n-D space
        cluster.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),cluster.angle);
        cluster.y = exy(:,2);
        cluster.crit = mean(exy(:,1));
        cluster.sign = E.sign;
        cluster.space = E.space;
        xy = xyrotate(DATA.xy(:,1),DATA.xy(:,2),cluster.angle);
        r = xy(:,1);
        if cluster.sign < 0
            id = find(xy(:,1) < cluster.crit);
            nid = find(xy(:,1) >= cluster.crit);
        else
            cluster.sign = 1;
            id = find(xy(:,1) > cluster.crit);
            nid = find(xy(:,1) <= cluster.crit);
        end
        e(1) = mean(DATA.energy(DATA.probe(1),id));
        e(2) = mean(DATA.energy(DATA.probe(1),nid));
        if e(2) > e(1);
            cluster.sign = -cluster.sign;
            cid = nid;
            nid = id;
            id = cid;
        end
        cluster.bmc = BimodalCoeff(r,1.5);
        cluster.dprime = CalcDprime(r(id),r(nid));
        cluster.hdip = HartigansDipTest(sort(r));
%        [a,b] = FindDip(r,DATA.energy(DATA.probe(1),:),'eval',cluster.crit);
        if ~isfield(E,'gmdprime')
            [a,b] = GMDip(r,DATA.energy(DATA.probe(1),:),'eval',cluster.crit);
            cluster.gmdprime = b.gmdprime;
            cluster.autodipsize = b.dipsize;
            cluster.dipsize = b.cdipsize;
        end
    elseif E.shape == 1
        angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);
        cluster.y = exy(:,2);
        crit = mean(exy(:,1));
        cluster.angle = angle;
        cluster.crit = crit;
        if isempty(p) %Var/E plot
            y = DATA.spkvar(ispk,:)./DATA.energy(ispk,:);
            xy = xyrotate(DATA.energy(ispk,:),y,angle);
           cluster.space = [5];
        elseif DATA.plottype == 2 %voltage pairs
        xy = xyrotate(DATA.AllV(p(1),p(2),:),DATA.AllV(p(3),p(4),:),angle);
           cluster.space = [2 p(1:4)];
        elseif ismember(DATA.plottype,[3 4])
            xy = xyrotate(DATA.TemplateScores(:,p(1)),DATA.TemplateScores(:,p(2)),angle);
            cluster.space = [3 p(1:2)];
        else
            xy = xyrotate(DATA.pcs(:,p(1)),DATA.pcs(:,p(2)),angle);
            cluster.space = [1 p(1:2)];
        end
        if crit < 0 && DATA.plottype == 1
            id = find(xy(:,1) < crit);
            nid = find(xy(:,1) >= crit);
            cluster.sign = -1;
        elseif crit <0 && mean(xy(:,1)) < 0
            id = find(xy(:,1) < crit);
            nid = find(xy(:,1) >= crit);
            cluster.sign = -1;
        else
        id = find(xy(:,1) > mean(exy(:,1)));
        nid = find(xy(:,1) <= mean(exy(:,1)));
            cluster.sign = 1;
        end
        e(1) = mean(DATA.energy(DATA.probe(1),id));
        e(2) = mean(DATA.energy(DATA.probe(1),nid));
        if e(2) > e(1);
            cluster.sign = -cluster.sign;
            cid = nid;
            nid = id;
            id = cid;
        end
        if  length(id) < 2
            
        end
        cluster.dprime = CalcDprime(xy(id,1),xy(nid,1));
        cluster.hdip = HartigansDipTest(sort(xy(:,1)));
%        [a,b] = FindDip(xy(:,1),DATA.energy(DATA.probe(1),:),'eval',crit);
        [a,b] = GMDip(xy,DATA.energy(DATA.probe(1),:));
        cluster.autodipsize = b.dipsize;
        cluster.dipsize = b.cdipsize;
        x = FitGaussMeans(xy,2,'clusterid',id);
        r = xy(:,1);
    else  %E.shape == 0
    if isempty(p) %Var/E plot
        xy(:,2) = DATA.spkvar(ispk,:)./DATA.energy(ispk,:);
        xy(:,1) = DATA.energy(ispk,:);
        cluster.space = [5];
    elseif DATA.plottype == 3
        xy(:,1) = DATA.TemplateScores(:,p(1));
        xy(:,2) = DATA.TemplateScores(:,p(2));
        cluster.space = [3 p(1:2)];
    elseif DATA.plottype == 2
        xy(:,1) = DATA.AllV(p(1),p(2),:);
        xy(:,2) = DATA.AllV(p(3),p(4),:);
        cluster.space = [2 p(1:4)];
    else
        xy(:,1) = DATA.pcs(:,p(1));
        xy(:,2) = DATA.pcs(:,p(2));
        cluster.space = [1 p(1:2)];
    end
    if isfield(E,'angle')
        cluster.angle = E.angle;
    else
        cluster.angle = 0;
    end
    xy = xyrotate(xy(:,1)-cx,xy(:,2)-cy,cluster.angle);
    r = ((xy(:,1))./rx).^2 + ((xy(:,2))./ry).^2;
    id = find(r < 1);
    nid = find(r>1);
    x = FitGaussMeans(xy,2,'clusterid',id);
    end
    cluster.bmc = BimodalCoeff(r,1.5);
    a = FitGaussMeans(xy,2);
    cluster.mahal = [a.mahal a.dprime];
    if isfield(cluster,'bestspace')
        cluster.mahal(3) = cluster.bestspace(2);
    else
        cluster.mahal(3) = 0;
    end
    if ismember(E.shape,[0 1]) && x.obj.Converged >= 0
        fprintf('Separation %.2f (manual) %.2f (auto)\n',x.mahal,a.mahal);
        cluster.mahal(3) = x.mahal;
        cluster.gmfit2dman = x.obj;
    end
    cluster.gmfit2d = a.obj;
    if isfield(E,'gmdprime')
        cluster.mahal(4) = E.gmdprime;
    else
        cluster.mahal(4) = 0; 
    end
    if cluster.cluster == 2
        cl.clst(id) = 3;
    else
    cl.clst(id) = 2;
    cl.clst(nid) = 1;
    end
    if length(DATA.uid) < length(cl.clst);
    DATA.nid = DATA.uid(ismember(DATA.uid, nid));
    DATA.clid = DATA.uid(ismember(DATA.uid, id));
    x = setdiff(1:length(cl.clst),DATA.uid);
    cl.clst(x) = 0;
    else
    DATA.nid = nid;
    DATA.clid = id;
    end
    DATA.clst = cl.clst;
    cl.id = DATA.clid;
    cl.nid = DATA.nid;
    
    cluster.ncut = length(id);
    E.h = [];
    if DATA.watchplots
        DATA = ReplotPCs(DATA,E);
        SetFigure('Spikes');
        hold off;
        DATA.spkst = DATA.uid;
        step = max([1 round(length(DATA.spklst)/1000)]);
        PlotSpikes(DATA,DATA.spklst(1:step:end));
    end
    cl.MeanSpike = PlotMeanSpike(DATA,'recalc');
    cluster.r = r;
 GetFigure('Vhist');
 subplot(1,1,1);
 t = find(DATA.spts == 0);
 if length(DATA.clid) < 10
     nbins(1) = 1;
     cluster.dropi = [0 0 0 0];
     cluster.trigsd = 0;
     cluster.minspke = min(DATA.energy(DATA.probe(1),cl.id))
     cluster.minspkvar = min(DATA.spkvar(DATA.probe(1),cl.id));
 elseif length(DATA.clid) < 100
     nbins(1) = 10;
 else
     nbins(1) = round(length(DATA.clid)./20);
 end
 if length(nid) < 10
     nbins(2) = 1;
 elseif length(nid) < 100
     nbins(2) = 10;
 else
     nbins(2) = round(length(nid)./20);
 end
     
 if nbins(1) > 2
     hold off;
%     V = DATA.AllV(DATA.probe(1),t,id);
     V = DATA.rV(id);
     [a,b] = hist(V,nbins(1));
     [c,d] = hist(DATA.rV(nid),nbins(2));
     bar(b,a,1);
     hold on;
     plot(d,c,'color',[0.5 0.5 0.5]);
     if DATA.Trigger(1) < 0
         tid = find(b <= DATA.Trigger(1));
         ntid = find(d <= DATA.Trigger(1));
     else
         tid = find(b >=DATA.Trigger(1));
         ntid = find(d >= DATA.Trigger(1));
     end
     gfit = FitGauss(b(tid),a(tid));
     ngfit = FitGauss(d(ntid),c(ntid));
     hold on;
     plot(b(tid),gfit.fitted,'r');
     nclose = sum(abs(V-DATA.Trigger(1)) < std(V)/10);
     cluster.dropi(1) = nclose./length(DATA.clid);
     cluster.dropi(2) = abs(mean(V)-DATA.Trigger(1))./std(V);
     if DATA.Trigger(1) < 0
         cluster.dropi(3) = (DATA.Trigger(1)-gfit.mean)./gfit.sd;
         cluster.dropi(4) = (DATA.Trigger(1)-ngfit.mean)./ngfit.sd;
     else
         cluster.dropi(3) = (gfit.mean-DATA.Trigger(1))./gfit.sd;
         cluster.dropi(4) = (DATA.Trigger(1)-ngfit.mean)./ngfit.sd;
     end
     cluster.trigsd = abs(gfit.sd);
     cluster.minspke = prctile(DATA.energy(DATA.probe(1),cl.id),1) .* 0.95;
     cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),cl.nid),1) .* 0.95;

     title(sprintf('P%d %d/%d spikes drop %.3f,%.2f',DATA.probe(1),length(id),length(id)+length(nid),cluster.dropi(1),cluster.dropi(3)));
     hold on;
     plot([DATA.Trigger(1) DATA.Trigger(1)],get(gca,'ylim'),'r');
 end
 
 if cluster.cluster == 2
     a = cluster;
     cluster = DATA.cluster;
     cluster.next = a;
 end
    


function MeanSpike = PlotMeanSpike(DATA, varargin)
    recalc = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'recalc',4)
            recalc = 1;
        end
        j = j+1;
    end
    GetFigure('MeanSpike');
    subplot(1,2,1);
    hold off;
    id = DATA.clid;
    nid = DATA.nid;
    if ~isfield(DATA,'MeanSpike') || recalc
        for j = 1:size(DATA.AllV,1)
            ms(j,:) = mean(DATA.AllV(j,:,id),3);
            mu(j,:) = mean(DATA.AllV(j,:,nid),3);
        end
        MeanSpike.ms = ms;
        MeanSpike.mu = mu;
        xc = corrcoef(ms(:),mu(:));
        MeanSpike.muxc = xc(1,2);
    else
        ms = DATA.MeanSpike.ms;
        mu = DATA.MeanSpike.mu;
    end

    chspk = DATA.probe(1)+ [-1:1]; 
    chspk = DATA.chspk;
    dj = 0;
%do csd first, then can do dvdt to CSD 
    if DATA.plotcsd
       ms = (diff(ms,2,1));
       mu = (diff(mu,2,1));
       csd = diff(DATA.AllV,2,1);
        chspk = chspk-1;
        dj = 1;
    end
    if DATA.plotdvdt
       ms = (diff(ms,1,2));
       mu = (diff(mu,1,2));
    end
    imagesc(ms);
    subplot(1,2,2);
    
    chspk = chspk(chspk >0 & chspk <= size(ms,1));
    

    hold off;
    voff = DATA.voffset - DATA.voffset(DATA.probe(1));
    for j = chspk
        plot(mu(j,:)+voff(j)/5,'color',[0.5 0.5 0.5]);
        hold on;
        plot(ms(j,:)+voff(j)/5,'r');
        if DATA.plotdvdt && DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(csd(j,:,id),1,2),[],3) var(diff(csd(j,:,nid),1,2),[],3)]));
        elseif DATA.plotdvdt
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(DATA.AllV(j,:,nid),1,2),[],3) var(diff(DATA.AllV(j,:,id),1,2),[],3)]));
        elseif DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(csd(j,:,id),[],3) var(csd(j,:,nid),[],3)]));
        else
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(DATA.AllV(j,:,nid),[],3) var(DATA.AllV(j,:,id),[],3)]));
        end
        if j == DATA.probe(1)
        plot(abs(dp),'m');
        else
        plot(abs(dp),'g');
        end
        text(size(ms,2),voff(j)/5,sprintf('%d',j+dj));
        peaks = find(diff(sign(diff(abs(dp)))) > 0);
        [a,b] = sort(abs(dp(peaks)),'descend');
        MeanSpike.dpmax(j,1:length(b)) = peaks(b);
        MeanSpike.dp(j,:) = dp;
        if recalc
        MeanSpike.vdprime(j,:) = dp;
        end
    end
    set(gca,'xlim',[1 size(ms,2)+1]);

function Scores = CalcScores(DATA, MeanSpike)
% this need to be modified to work like TemplatePlot now, with
%T being for all probes, but only calculating scores for chspk
 
if isfield(MeanSpike,'ms')
    T = MeanSpike.ms;
    mT = MeanSpike.mu;
else
    T = MeanSpike;
end
    if DATA.chspk(end) == DATA.nprobes
    else
        dspk = [DATA.chspk DATA.chspk(end) + 1];
    end
    if DATA.chspk(end) == DATA.nprobes
        dspk = [DATA.chspk(1) - 1 DATA.chspk];
        csdspk = [DATA.chspk(1) - 2 DATA.chspk(1) - 1 DATA.chspk];
    elseif DATA.chspk(1) == 1
        dspk = [DATA.chspk DATA.chspk(end) + 1];
        csdspk = [DATA.chspk DATA.chspk(end) + 1 DATA.chspk(end) +2];
    else
        dspk = [DATA.chspk DATA.chspk(end) + 1];
        csdspk = [DATA.chspk(1) - 1 DATA.chspk DATA.chspk(end) + 1];
    end

    dv = diff(DATA.AllV(dspk,:,:),1,1);
    mdv = diff(T(dspk,:,:),1,1);
    csd = diff(DATA.AllV(csdspk,:,:),2,1);
    mcsd = diff(T(csdspk,:,:),2,1);
    meanV = repmat(mean(DATA.AllV(DATA.chspk,:,:),2),[1  size(DATA.AllV,2) 1]);
    vid = 1:min([size(T,2) size(DATA.AllV,2)]); %% in case old cut used different range
    for j = 1:length(DATA.chspk)
        c = DATA.chspk(j);
        T(c,vid) = T(c,vid) - mean(T(c,vid));
        Scores(j,1,:) = squeeze(T(c,vid)) * squeeze(DATA.AllV(c,vid,:) - meanV(j,vid,:));
        Scores(j,3,:) = squeeze(mdv(j,:)) * squeeze(dv(j,:,:));
        Scores(j,4,:) = squeeze(mcsd(j,:)) * squeeze(csd(j,:,:));
        Scores(j,2,:) = squeeze(diff(T(c,vid),1,2)) * squeeze(diff(DATA.AllV(c,vid,:),1,2));
        if size(MeanSpike.vdprime,1) >= c
        dp = MeanSpike.vdprime(c,:);
        Scores(j,5,:) = dp * squeeze(DATA.AllV(c,:,:) - repmat(MeanSpike.mu(c,:),[1 1 size(DATA.AllV,3)]));
        else
        Scores(j,5,:) = 0;
        end
    end

function Labels = TemplateLabels(DATA, usestd)
    ispk = DATA.probe(1);
    chspk = DATA.probe(1)+ [-1:1];
    chspk = chspk(chspk >0 & chspk <= size(DATA.AllV,1));
    if usestd
    Labels{1} = sprintf('1r');
    Labels{2} = sprintf('2r');
    Labels{3} = sprintf('1dt');
    Labels{4} = sprintf('2dt');
    Labels{5} = sprintf('?');
    Labels{6} = sprintf('?');
    Labels{7} = sprintf('?');
    Labels{8} = sprintf('2dt');
    Labels{9} = sprintf('?');
    Labels{10} = sprintf('2dt');
    Labels{11} = sprintf('?');
    Labels{12} = sprintf('?');
    return;
    end
    Labels{1} = sprintf('%d:r',ispk);
    Labels{2} = sprintf('sum');
    Labels{8} = sprintf('%d:dt',ispk);
    Labels{9} = sprintf('%d:dy',ispk);
    Labels{5} = sprintf('%d:csd',chspk(1));
    Labels{6} = sprintf('%d:dp',ispk);
    Labels{7} = sprintf('%d:dp',chspk(1));
    Labels{3} = sprintf('%d:r',chspk(1));
    if length(chspk) > 2
        Labels{4} = sprintf('%d:r',chspk(3));
    else
        if max(chspk)  == size(DATA.AllV,1)
            xspk = min(chspk)-1;
        else
            xspk = max(chspk)+1;
        end
        Labels{4} = sprintf('%d:r',xspk);
    end
    Labels{10} = sprintf('sumdt');
    Labels{11} = sprintf('sumdp',ispk);
    Labels{12} = sprintf('sumdy',ispk);
    Labels{13} = sprintf('PC1',ispk);
    Labels{14} = sprintf('PC2',ispk);
    Labels{15} = sprintf('PC3',ispk);

function out = TemplatePlot(DATA, varargin)

    projectout = 0;
    calcdips = 1;
    usemean = 0;
    usestd = 0;
    j = 1;
    while  j <= length(varargin)
        if strncmpi(varargin{j},'nodip',4)
            calcdips = 0;
        elseif strncmpi(varargin{j},'projectout',7)
            projectout = 1;
        elseif strncmpi(varargin{j},'usemean',7)
            usemean = 1;
        elseif strncmpi(varargin{j},'stdtemplate',7)
            usestd = 1;
        end
        j = j+1;
    end
    
    DATA.usestdtemplates = usestd;
    set(DATA.toplevel,'Name','Calculating Templates...');
    drawnow;
    id = DATA.clid;
    nid = DATA.nid;
    if usemean 
        ms = DATA.MeanSpike.ms;
        mu = DATA.MeanSpike.mu;
    else
    for j = size(DATA.AllV,1):-1:1
        ms(j,:) = mean(DATA.AllV(j,:,id),3);
        mu(j,:) = mean(DATA.AllV(j,:,nid),3);
        mu(j,:) = mu(j,:)-mean(mu(j,:));
        ms(j,:) = ms(j,:)-mean(ms(j,:));
%        xc(j,1,:) = squeeze(TemplateScores(j,1,:))./(DATA.spkvar(j,:)' .* std(ms(j,:)));
%        xc(j,2,:) = squeeze(TemplateScores(j,2,:))./(DATA.spkvar(j,:)' .* std(mu(j,:)));
    end
    end

    ispk = find(DATA.chspk == DATA.probe(1));
    if usestd
        j = DATA.probe(1);
        meanV = repmat(mean(DATA.AllV(j,:,:),2),[1  size(DATA.AllV,2) 1]);
        TemplateScores(ispk,1,:) = DATA.StdTemplate(1,:) * squeeze(DATA.AllV(j,:,:) - meanV);
        TemplateScores(ispk,2,:) = DATA.StdTemplate(2,:) * squeeze(DATA.AllV(j,:,:) - meanV);
        TemplateScores(ispk,3,:) = squeeze(diff(DATA.StdTemplate(1,:),1,2)) * squeeze(diff(DATA.AllV(j,:,:),1,2));
        TemplateScores(ispk,7,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(DATA.AllV(j,:,:),1,2));
        TemplateScores(ispk,8,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(DATA.AllV(j,:,:),1,2));
    else
    for k = length(DATA.chspk):-1:1
        j = DATA.chspk(k);
        meanV = repmat(mean(DATA.AllV(j,:,:),2),[1  size(DATA.AllV,2) 1]);
        TemplateScores(k,1,:) = squeeze(ms(j,:)) * squeeze(DATA.AllV(j,:,:) - meanV);
        TemplateScores(k,2,:) = squeeze(mu(j,:)) * squeeze(DATA.AllV(j,:,:) - meanV);
        dp(j,:) = (mean(DATA.AllV(j,:,id),3)-mean(DATA.AllV(j,:,nid),3))./sqrt(mean([var(DATA.AllV(j,:,nid),[],3) var(DATA.AllV(j,:,id),[],3)]));
    end
    clear meanV;
    dpcrit = 1;
    for k = length(DATA.chspk):-1:1
        j = DATA.chspk(k);
%        dp * squeeze(DATA.AllV(j,:,:) - repmat(mu(j,:),[1 1 size(DATA.AllV,3)]));
        TemplateScores(k,3,:) = squeeze(dp(j,:)) * squeeze(DATA.AllV(j,:,:) - repmat(mu(j,:),[1 1 size(DATA.AllV,3)]));
        id = find(abs(dp(j,:)) > dpcrit);
        if length(id)
            TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(DATA.AllV(j,id,:));
            TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(DATA.AllV(j,id,:));
        else
            [a,b] = sort(dp(j,:),'descend');
            id = b(1:2);
            TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(DATA.AllV(j,id,:));
            TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(DATA.AllV(j,id,:));
        end
        TemplateScores(k,6,:) = squeeze(diff(ms(j,:),1,2)) * squeeze(diff(DATA.AllV(j,:,:),1,2));
    end
    if max(DATA.chspk) < DATA.nprobes
        dyspk = [DATA.chspk DATA.chspk(end)+1];
    else
        dyspk = [DATA.chspk(1)-1 DATA.chspk];
    end
    dv = diff(DATA.AllV(dyspk,:,:),1,1);  %diff along length
    mdv = diff(ms(dyspk,:),1,1);
    for k = 1:size(dv,1)
        TemplateScores(k,7,:) = squeeze(mdv(k,:)) * squeeze(dv(k,:,:));
    end
    clear dv;
    clear mdv;

    if min(dyspk) > 1
        csdspk = [dyspk(1)-1 dyspk];
    else
        csdspk = [dyspk dyspk(end)+1];
    end

    csd = diff(DATA.AllV(csdspk,:,:),2,1);
    mcsd = diff(ms(csdspk,:),2,1);
    for j = 1:size(csd,1)
        TemplateScores(j,8,:) = squeeze(mcsd(j,:)) * squeeze(csd(j,:,:));            
    end
    clear csd;
    clear mcsd;
    end

    if length(DATA.chspk) > 2
    chspk = DATA.chspk;
    else
    chspk = DATA.probe(1)+ [-1:1];
    end
    chspk = chspk(chspk >0 & chspk <= size(DATA.AllV,1));
    if min(chspk) > 1
    csdspk = chspk-1;
    else
        csdspk = chspk;
    end
    
    if max(chspk)  == size(DATA.AllV,1)
        xspk = min(chspk)-1;
    else
        xspk = max(chspk)+1;
    end
        
    if projectout  %doesn't seem much use...., and uses too much memory
%project out the templates, then redo the pca
    mg = sum(ms.^2,2);
    G = sum(DATA.AllV.*repmat(ms,[1 1 size(DATA.AllV,3)]),2)./repmat(mg,[1 1 size(DATA.AllV,3)]);
    nv = DATA.AllV - repmat(ms,[1 1 size(DATA.AllV,3)]) .* repmat(G,[1 size(DATA.AllV,2) 1]);
    TV = nv(chspk(1),:,:);
    for j = 2:length(chspk)
        TV = cat(2,TV,nv(chspk(j),:,:));
    end
    TV = squeeze(TV)';
    [pc, E] = eig(cov(TV));
    pc = fliplr(pc); %put largest first;
    pcs = TV*pc;
    end
    TMPL.pcs(:,1) = sum(TemplateScores(:,1,:));
    ispk = find(DATA.chspk == DATA.probe(1));

    if projectout == 2
        TMPL.pcs(:,2:9) = pcs(:,1:8);
    elseif usestd
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);  %1r
        TMPL.pcs(:,2) = TemplateScores(ispk,2,:); %
        TMPL.pcs(:,3) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,4) = TemplateScores(ispk,7,:); 
        TMPL.pcs(:,5) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,6) = TemplateScores(ispk,2,:);
        TMPL.pcs(:,7) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,8) = TemplateScores(ispk,7,:); %2dt 
        TMPL.pcs(:,9) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,10) = TemplateScores(ispk,7,:); %2dt
        TMPL.pcs(:,11) = TemplateScores(ispk,3,:);
        TMPL.pcs(:,12) = TemplateScores(ispk,7,:); 
    else
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,2) = sum(TemplateScores(:,1,:));
        TMPL.pcs(:,8) = TemplateScores(ispk,6,:); %dvdt
        TMPL.pcs(:,9) = TemplateScores(ispk,7,:); %dvdy
        TMPL.pcs(:,5) = TemplateScores(ispk,8,:); %csd
        TMPL.pcs(:,6) = TemplateScores(ispk,3,:); %dprime weighted
        TMPL.pcs(:,7) = TemplateScores(1,3,:); %dprime weighted
%This needs work if length(chspk) > 3
        if length(chspk) > 1
            TMPL.pcs(:,3) = TemplateScores(1,1,:);
        end
        if usestd
            TMPL.pcs(:,4) = TemplateScores(ispk,1,:);
        elseif length(chspk) > 2
            if ispk == 3
            TMPL.pcs(:,4) = TemplateScores(2,1,:);
            else
            TMPL.pcs(:,4) = TemplateScores(3,1,:);
            end
        else
            TMPL.pcs(:,4) = TemplateScores(end,1,:);
        end
        TMPL.pcs(:,10) = sum(TemplateScores(:,6,:));
        TMPL.pcs(:,11) = sum(TemplateScores(:,3,:)); %dprime weighted
        TMPL.pcs(:,12) = sum(TemplateScores(:,7,:)); %sum dy
        if projectout
            TMPL.pcs(:,13:15) = pcs(:,1:3);
        end
    end
    DATA.TemplateLabels = TemplateLabels(DATA, usestd);
    TMPL.dvdt = 0;
    TMPL.csd = 0;
    TMPL.clid = id;
    TMPL.nid = nid;
    TMPL.toplevel = DATA.toplevel;
    TMPL.pcplots = DATA.pcplots;
    TMPL.clplot = DATA.clplot;
    TMPL.plottype = 1;
    DATA.TemplateScores = TMPL.pcs;
    for j = 1:size(DATA.TemplateScores,2)
        DATA.TemplateScores(:,j) = DATA.TemplateScores(:,j)./std(DATA.TemplateScores(:,j));
    end
    
    if calcdips
        for j = 1:5
            TMPL.dipvals(j) = HartigansDipTest(sort(TMPL.pcs(:,j)));
        end
        DATA.tmpdips = CalculateTemplateDips(DATA);
    theta = 0:pi/36:pi * 35/36;
    for j = 1:length(theta)
        for k = 2:5;
        xy = xyrotate(TMPL.pcs(:,1),TMPL.pcs(:,k),theta(j));
        rdip(j,k) = HartigansDipTest(xy(:,1));
        end
    end
    else
        DATA.tmpdips = zeros(1,8);
    end

    if ~ismember(DATA.plottype,[3 4])
    DATA.plottype = 3;
    end

    if DATA.watchplots
    SetFigure('TemplateScores');
    subplot(1,1,1);
    PlotTemplateScores(TMPL,TemplateScores, [1 2]);


    DATA = ReplotPCs(DATA,[]);
    end
%currently sets DATA.TemplateScores, DATA.TemplateLabels, and DATA.tmpdips
%if output is requestsed, don't set the figure userdata
    if nargout
        set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));
        out = DATA.TemplateScores;
        return;
    end
    set(DATA.toplevel,'UserData',DATA);
    set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));

    
 function [bs, as] = CalculateTemplateDips(DATA)
     bs = [];
     as = [];
     if ~isfield(DATA,'TemplateScores')
         return;
     end
     for j = 1:size(DATA.tmplots,1)
         if max(DATA.tmplots(j,:)) > size(DATA.TemplateScores,2)
             as(j) = 0;
             bs(j) = 0;
         else
         [as(j),bs(j)] = BestAngle(DATA.TemplateScores(:,DATA.tmplots(j,1)),DATA.TemplateScores(:,DATA.tmplots(j,2)),1);
         end
     end

function DATA = ReplotPCs(DATA,E, varargin)
        
    
meanpos = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'showpos')
        j = j+1;
        meanpos = varargin{j};
        j = j+1;
        dimpos = varargin{j};
    end
    j = j+1;
end
figure(DATA.toplevel)
clusterplot = 0;
if isfield(DATA,'dipvals')
    dipvals = DATA.dipvals;
else
    dipvals = zeros(1,8);
end
if DATA.usebmi == 0
dipvals = zeros(1,8);
end
if DATA.plottype == 2;
    vpts = DATA.vpts;
elseif DATA.plottype == 6;
    vpts = DATA.dvpts;
    DATA.dV= diff(DATA.AllV,1,2);
end
if isempty(E) && isfield(DATA,'cluster') && isfield(DATA.cluster,'space')
    if DATA.cluster.space(1) == DATA.plottype %cluster is in this space
        E = BoundaryFromCluster(E,DATA.cluster);
    end
end

plots = DATA.pcplots(1:8,:);
if DATA.plottype == 3
    if DATA.usestdtemplates
    plots = DATA.tmplots(17:24,:);
    else
    plots = DATA.tmplots(1:8,:);
    end
    dipvals = DATA.tmpdips;
elseif DATA.plottype == 4
    plots = DATA.tmplots(9:16,:);
end

   
if DATA.usegmcid
    if isfield(DATA,'gmcid')
        clid = DATA.gmcid;
    elseif isfield(E,'bestcl')
        clid = E.bestcl;
    else
    clid = DATA.clid;
    DATA.usegmcid = 0; %don't use it if not defined
    end
else
    clid = DATA.clid;
    if isfield(DATA,'clst') 
        clid = DATA.clst;
    end
end
for j = 1:size(plots,1)
    subplot(2,4,j);
    if ismember(DATA.plottype, [3 4 5])
        PlotPCs(DATA.TemplateScores,plots(j,1),plots(j,2),DATA.clplot,clid,DATA.colors,'fixrange');
        title(sprintf('%s vs %s',...
            DATA.TemplateLabels{plots(j,1)},...
            DATA.TemplateLabels{plots(j,2)}));
        if isfield(E,'gmfit') && E.space(1) == 6 && E.space(2) ==4 
            [a,xi] = ismember(plots(j,1),DATA.tmplspace(1,:));
            [b,yi] = ismember(plots(j,2),DATA.tmplspace(1,:));
            if a && b
                plot(E.gmfit.mu(1,xi),E.gmfit.mu(1,yi),'g+','markersize',10,'linewidth',2);
                plot(E.gmfit.mu(2,xi),E.gmfit.mu(2,yi),'g+','markersize',10,'linewidth',2);
            end
            
        end
        set(gca,'UserData',plots(j,:));
    elseif ismember(DATA.plottype,[6])
        PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),DATA.clplot,clid,DATA.colors);
        set(gca,'UserData',vpts(j,:));
        if DATA.usegmcid && sum(ismember([vpts(1,3) vpts(3,4)],DATA.vspace) == 2)
        end
    elseif ismember(DATA.plottype,[2 6])
        PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),DATA.clplot,clid,DATA.colors,'fixrange');
        set(gca,'UserData',vpts(j,:));
        if DATA.usegmcid && sum(ismember([vpts(j,2) vpts(j,4)],DATA.vspace)) == 2
            k = find(DATA.vspace == vpts(j,2));
            m = find(DATA.vspace == vpts(j,4));
            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)
        end
    else
        PlotPCs(DATA.pcs,DATA.pcplots(j,1),DATA.pcplots(j,2),DATA.clplot,clid,DATA.colors, 'fixrange');
        set(gca,'UserData',DATA.pcplots(j,:));
        if DATA.usegmcid && sum(ismember(DATA.pcplots(j,:),[1:4])) == 2 && isfield(DATA.cluster,'gmfit')
            k = find(DATA.pcspace == DATA.pcplots(j,1));
            m = find(DATA.pcspace == DATA.pcplots(j,2));
            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)
        end
    end
    axis('tight');          
    xl = get(gca,'Xlim');
    yl = get(gca,'Ylim');
    text(mean(xl),yl(2),sprintf('%.2f',dipvals(j)),'VerticalAlignment','top','color','k');
    if isempty(E) || isempty(E.pcplot) || E.shape == 2
        clusterplot = 0;
    elseif DATA.plottype == 2 && sum(DATA.vpts(j,:) == E.pcplot) == 4
        clusterplot = j;
    elseif ismember(DATA.plottype,[1 3]) && sum(plots(j,:) == E.pcplot(1:size(plots,2))) == 2
        clusterplot = j;
    end
end
if DATA.plottype == 1
if isfield(DATA,'alldips')
    if size(DATA.alldips,2) > 1
        b = max(DATA.alldips');
    else
        b = DATA.alldips;
    end
   subplot(2,4,1);
   text(-0.2,1.1,sprintf('%s %.2f, dt%.2f, csd%.2f',num2str(DATA.probe),b(1).*100,b(2).*100,b(3).*100),'units','norm');
elseif DATA.dvdt
        text(-0.2,0.1,'dvdt','units','norm');
elseif DATA.csd
    text(-0.2,0.1,'csd','units','norm');
end
end
    
if clusterplot
    subplot(2,4,clusterplot); %need to find right graph
    hold on;
    E.color = 'r';
    DATA.elmousept.h = DrawEllipse(E);
    if E.shape == 1
        tmp = E;
        tmp.pos(1) = mean(E.pos([1 3])) + diff(E.pos([2 4]))/2;
        tmp.pos(3) = mean(E.pos([1 3])) - diff(E.pos([2 4]))/2;
        tmp.pos(2) = mean(E.pos([2 4])) - diff(E.pos([1 3]))/2;
        tmp.pos(4) = mean(E.pos([2 4])) + diff(E.pos([1 3]))/2;
        h = DrawEllipse(tmp,'r');
        set(h,'linestyle',':');
    end
    hold off; 
end
if isfield(DATA,'energy')
GetFigure('VarE');
subplot(1,1,1);
PlotVarE(DATA);
end


function AddMarkToPCs(pos, space, plots, varargin)
    c = 'g';

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'color',4)
            j = j+1;
            c = varargin{j};
        end
        j = j+1;
    end

   for j = 1:size(plots,1)
   [a,xi] = ismember(plots(j,1),space(1,:));
   [b,yi] = ismember(plots(j,2),space(1,:));
   if a && b
       subplot(2,4,j);
       hold on;
       for k = 1:size(pos,1)
           plot(pos(k,xi),pos(k,yi),'+','markersize',10,'linewidth',2,'color',c);
           hold on;
       end
   end            
   end

function PlotVarE(DATA)
        
    colors = 'brgmck';
c = DATA.probe(1);
if DATA.clplot == 1
DensityPlot(DATA.energy(c,DATA.nid),DATA.spkvar(c,:)./DATA.energy(c,:));
else
hold off;
if DATA.usegmcid && isfield(DATA,'gmcid')
    cls = unique(DATA.gmcid);
    for j = length(cls)
        id = find(DATA.gmcid ==cls(j));
        plot(DATA.energy(c,id),DATA.spkvar(c,id)./DATA.energy(c,id),'.','markersize',1,'color',colors(j));
    end
elseif isfield(DATA,'nid')
plot(DATA.energy(c,DATA.nid),DATA.spkvar(c,DATA.nid)./DATA.energy(c,DATA.nid),'.','markersize',1);
hold on;
plot(DATA.energy(c,DATA.clid),DATA.spkvar(c,DATA.clid)./DATA.energy(c,DATA.clid),'r.','markersize',1);
else
plot(DATA.energy(c,:),DATA.spkvar(c,:)./DATA.energy(c,:),'.','markersize',1);
end
end
title(sprintf('%d/%d Spikes',length(DATA.clid),size(DATA.AllV,3)));


function DATA = RestrictTimeRange(DATA, t)
   DATA.uid = find(DATA.t >= t(1) & DATA.t <= t(2));
    DATA.restricttimerange = t;
    if isfield(DATA,'clid')
       DATA.clid = DATA.clid(ismember(DATA.clid,DATA.uid));
       DATA.nid = DATA.nid(ismember(DATA.nid,DATA.uid));
    end
    set(DATA.toplevel,'UserData',DATA);

function PlotSpikes(DATA,spkid, varargin)

dvdt = 1;
fixy = 0;
colors {1} = [0.5 0.5 0.5];
colors {2} = [1 0 0];
colors {3} = [0 1 0];
scale = DATA.plotspk.muscale;
yl = [];
j = 1;
quicktest = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fixy',3)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            yl = varargin{j};
        else
            fixy = 1;
        end
    end
    j = j+1;
end

SetFigure('Spikes');
if fixy
    yl = get(gca,'ylim');
end
subplot(1,1,1);
if ~isfield(DATA,'clid')
    DATA.clid = [];
    DATA.nid = 1:size(DATA.AllV,3);
end

if DATA.plotspk.oneprobe
    chspk = DATA.probe(1);
else
    chspk = DATA.plotspk.probes;
end

if isempty(spkid)
    for c= chspk;
        if c > 0 & c <= size(DATA.AllV,1)
        for j = 1:2
            h = plot([0 30],[0 0],'color',colors{j});
        end
        hold on;
        end
    end
    if fixy
        set(gca,'ylim',yl);
    end
    if DATA.plotspk.bytrial
        nt = DATA.currenttrial;
        title(sprintf('Trial %d: No Spikes  ed%.2f',DATA.Expt.Trials(nt).Trial,DATA.Expt.Trials(nt).ed));
    else
        title(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...
            spkid(1),spkid(end),DATA.t(spkid(1)),DATA.t(spkid(end)),length(id),length(spkid)));
    end
    return;
end


if max(spkid) > size(DATA.AllV,3)
    spkid = spkid(spkid<=size(DATA.AllV,3));
end
spkid = spkid(spkid > 0);

if isempty(spkid)
    ids{2} = DATA.clid;
    ids{1} = DATA.nid;
    id = ids{2};
   nid = ids{1};
elseif DATA.usegmcid
        nc = unique(DATA.gmcid);
        for j = 1:length(nc)
            id = find(DATA.gmcid(spkid) == nc(j));
            ids{j} = spkid(id);
        end
else
    id = find(ismember(spkid,DATA.clid));
    nid = find(ismember(spkid,DATA.nid));
    ids{2} = spkid(id);
    ids{1} = spkid(nid);
end
ispk = DATA.probe;
for j = 1:length(ids)
    V{j} = DATA.AllV(:,:,ids{j});
end

voff = DATA.voffset - DATA.voffset(ispk(1));
if DATA.plotcsd && DATA.csd && isfield(DATA,'AllCSD')
for j = 1:length(ids)
    V{j} = DATA.AllCSD(:,:,ids{j});
end
    if min(ispk) > 1
    ispk = ispk-1;
    end
voff = (DATA.voffset - DATA.voffset(ispk(1)));
elseif DATA.plotdvdt
for j = 1:length(ids)
    V{j} = diff(DATA.AllV(:,:,ids{j}),1,2);
end
end

if DATA.plotspk.submax
    for j = 1:size(V)
        V{j} = V{j} - repmat(max(V{j},[],2),[1 size(V{j},2) 1]);
    end
elseif DATA.plotspk.submin
    for j = 1:size(V)
        V{j} = V{j} - repmat(min(V{j},[],2),[1 size(V{j},2) 1]);
    end
elseif DATA.plotspk.submean
    for j = 1:size(V)
        V{j} = V{j} - repmat(mean(V{j},[],2),[1 size(V{j},2) 1]);
    end
end

V{j} = V{j} .* scale;

if quicktest
    l = size(V{1},2);
    hold off;
    x = [1:l NaN];
for c= chspk;
    if c > 0 & c <= size(DATA.AllV,1)
        for j = 1:length(V)
        nV = squeeze(V{j}(c,:,:)) + voff(c);
        if size(nV,1) > 1
            h = plot([0 30],[0 0],'color',colors{j});
            nV(l+1,:) = NaN;
            set(h,'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));
        else
            plot(1:l,nV,'color',colors{j});
        end
        hold on;
        end
        text(size(V{1},2),voff(c),sprintf('%d',c));
    end

end
else
hold off;


for c= [ispk(1)-1:1:ispk(1)+1];
    if c > 0 & c <= size(DATA.AllV,1)
        for j = 1:size(V)
            plot(squeeze(voff(c)+V{j}(c,:,:)),'color',colors{j});
            hold on;
        end
        text(size(V{1},2),voff(c),sprintf('%d',c));
    end
end
end
if length(yl) == 2
    set(gca,'ylim',yl);
end
hold off;

if DATA.plotspk.bytrial
    nt = DATA.currenttrial;
title(sprintf('Trial %d: %d/%d ed%.2f',DATA.Expt.Trials(nt).Trial, length(id),length(spkid),DATA.Expt.Trials(nt).ed));
else
title(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...
    spkid(1),spkid(end),DATA.t(spkid(1)),DATA.t(spkid(end)),length(id),length(spkid)));
end
        


function PlotTemplateScores(DATA, TemplateScores, probes)

p = 3;
p = probes(2);
slope = mean(squeeze(TemplateScores(p,1,:)))./mean(squeeze(TemplateScores(p,2,:)));
diffs = squeeze(TemplateScores(p,1,:) - slope.* TemplateScores(p,2,:));
hold off;
p = probes(1);
plot(squeeze(TemplateScores(p,1,DATA.nid)),diffs(DATA.nid),'.','markersize',1)
hold on;
plot(squeeze(TemplateScores(p,1,DATA.clid)),diffs(DATA.clid),'r.','markersize',1)
hold off;



function [dp, res] = MaxDprime(x, varargin)

    y = [];
    step = 10;
ndim = 1;
if length(varargin) && isnumeric(varargin{1}) && length(varargin{1}) == length(x)
    y = varargin{1};
    j = 2;
end

if length(y)
    theta = 0:pi/36:pi * 35/36;
    
    for j = 1:length(theta)
        xy = xyrotate(x,y,theta(j));
        rdip(j) = HartigansDipTest(xy(:,1));
        [dp, a] = MaxDprime(xy(:,1));
        dps(j) = a.dp;
        alldp(j,:) = a.dps;
        res.prcs(j) = a.dprc;
    end
    dp = max(dps);
    res.dps = dps;
    res.dips = rdip;
    res.p
else
    pa = sort(x);

    k = 1;
    for j = step:step:length(pa)-step
        dps(k) = (mean(pa(1:j))-mean(pa(j+1:end)))./sqrt(mean([var(pa(1:j)) var(pa(j+1:end))]));
        crit(k) = mean(pa(j:j+1));
        k=k+1;
    end
    id = find(diff(sign(diff(abs(dps)))) < 0); %local maxima=
    res.dps = dps;
    res.crit = crit;
    if isempty(id)
        dp = 0;
        res.dprc = 0;
        res.dp = 0;
        res.maxid =1;
    else
    [res.dp,b] = max(abs(dps(id)));
    res.dprc = (id(b)*step)./length(x);
%may add a check that mean diff(abs(dp)) around here is positive, so that noise on a 
%negative crossing doesn't qualify
    dp = res.dp;
    res.maxid = id(b);
    if res.dprc < 0.01
        dp = 0;
    end
    end
end



    
function [dip, details] = FindDip(values, energy, varargin)
evalcrit = [];
plottype = 0;
domix = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'eval',4)
        j = j+1;
        evalcrit = varargin{j};
    elseif strncmpi(varargin{j},'gmix',4)
        domix = 1;
    elseif strncmpi(varargin{j},'plot',4)
        plottype = 1;
    end
    j = j+1;
end
if diff(size(values)) > 0
    v = sort(values');
else
    v = sort(values);
end
x = diff(v);
w = 10;

if length(v) < 500
first = 50;
stpt=50;
smsd = 20;
else
first = 100;
stpt = 150;
smsd = 40;
end


xs = smooth(x,smsd,'gauss');
[a,b] = min(xs);
id = [];
c = prctile(xs,90).*10;
while isempty(id)
id = find(xs(1:b) > c);
c = c/2;
end
%first = id(end);


if domix
    G = gmdistribution.fit(values,2,'Options',statset('MaxIter',1000)); %2 Gaussians, 1 dimension
    details.gmdprime = abs(diff(G.mu))./sqrt(mean(G.Sigma));
    x = linspace(min(v), max(v),500);
    y = pdf(G, x');
    details.gxy(:,1) = x;
    y = pdf('norm',x,G.mu(1), sqrt(G.Sigma(1))) .* G.PComponents(1);
    details.gxy(:,2) = y;
    z = pdf('norm',x,G.mu(2), sqrt(G.Sigma(2))) .* G.PComponents(2);
%find the dip in the fitted sum between teh two means
    gid = find(x > min(G.mu) & x < max(G.mu));
    if isempty(gid)
        dip(6) = mean(G.mu);
    else
        [aa, peak] = min(z(gid)+y(gid));
        dip(6) = x(peak+gid(1)-1);
    end
    peak = find(diff(sign(z-y)) ~= 0);
    [aa,bb] = max(z(peak)+y(peak));
    details.gxy(:,3) = z;
    details.gmfit = G;
    if isempty(bb)
        dip(5) = 0;
    else
    dip(5) = x(peak(bb));
    end
end
c = prctile(xs,90).*10;
id = find(xs(b:end) > c);
while isempty(id)
    c = c/2;
    id = find(xs(b:end) > c);
end
    
last = id(1)+b-1;
last = length(v)-100;
id = convhull(v([first:last]+1),xs(first:last));
a = find(diff(id) < 0);
iid = unique([1; id(a(end)+1:end)+first; id(2:a(1))+first-1; length(xs)]);
hull = interp1(v(iid+1),xs(iid),v(2:end));
dx = xs - smooth(xs,400);
dx = smooth(xs-hull,100);
dx = xs-hull;
sdx = smooth(dx,100);
    [c,d] = min(dx(b:end));
    d = length(sdx)-stpt;
    if b < d
    [a,dipid] = max(sdx(b:d));
    j = b +dipid-1;
    else
        j = d;
    end
    dip(1) = v(j+1);
    dipsize(1) = (xs(j)-hull(j))./hull(j);
    adipsize(1) = xs(j)-hull(j);
    dipid(1) = j;
    if nargout > 1
        details.dprime(1) = CalcDprime(v(1:j),v(j+1:end));
    end
    [a,dipid] = max(dx(b:d));
    j = dipid+b-1;
    dip(2) = v(j+1);
    dipsize(2) = (xs(j)-hull(j))./hull(j);
    adipsize(2) = xs(j)-hull(j);
    dipid(2) = j;
    if nargout > 1
        details.dprime(2) = CalcDprime(v(1:j),v(j+1:end));
    end
    [c,d] = min(dx(1:b));
    d = stpt;
    [a,dipid] = max(sdx(d:b));
    j = d+dipid-1;
    dip(3) = v(j+1);
    dipsize(3) = (xs(j)-hull(j))./hull(j);
    adipsize(3) = xs(j)-hull(j);
    dipid(3) = j;
    if nargout > 1
        details.dprime(3) = CalcDprime(v(1:j),v(j+1:end));
    end
    [a,dipid] = max(dx(d:b));
    j = d+dipid-1;
    dip(4) = v(j+1);
    dipid(4) = j;
    dipsize(4) = (xs(j)-hull(j))./hull(j);
    adipsize(4) = xs(j)-hull(j);
    if nargout > 1
        details.dprime(4) = CalcDprime(v(1:j),v(j+1:end));
    end

    e(1) = mean(energy(find(values > mean(dip([1 2])))));
    e(2) = mean(energy(find(values < mean(dip([3 4])))));
    sgn = 0;
    if mean(dipsize([3 4])) > 2 * mean(dipsize([1 2])) && mean(adipsize([3 4])) > mean(adipsize([1 2]))
        sgn = -1;
    elseif mean(dipsize([1 2])) > 2 * mean(dipsize([3 4])) && mean(adipsize([1 2])) > mean(adipsize([3 4]))
        sgn = 1;
    end
    if sgn == 0 && e(2) > 2 * e(1)
        sgn = -1;
    elseif sgn == 0 && e(1) > 2 * e(2)
        sgn = 1;
    end
    if sgn == 0 && mean(values) < 0 && skewness(values) < 0
        sgn = -1;
    end
        
    if sgn < 0
        dip(1:4) = dip([3 4 1 2]);
        dipsize = dipsize([3 4 1 2]);
        if nargout > 1
        details.dprime = details.dprime([3 4 1 2]);
        end
        details.sign = -1;
    else
        details.sign = 1;
    end
if plottype == 1
    f = gcf;
    GetFigure('Dips');
    hold off;
    plot(v(2:end),xs);
    hold on;
    plot(v(2:end),sdx,'g');
    plot(v(2:end),hull,'r');
    plot([dip(1) dip(1)],get(gca,'ylim'),'g');
    plot([dip(2) dip(2)],get(gca,'ylim'),'g');
    plot([dip(3) dip(3)],get(gca,'ylim'));
    plot([dip(4) dip(4)],get(gca,'ylim'));
    if domix & ~isnan(dip(5))
    plot([dip(5) dip(5)],get(gca,'ylim'),'r');
    end
    GetFigure(f);
end

details.dipsize = dipsize;
for j = 1:length(evalcrit)
    details.crits(j) = evalcrit(j);
    [a,id] = min(abs(evalcrit(j) - v(2:end)));
    details.cdipsize(j) = (xs(id)-hull(id))./hull(id);
end



function dp = CalcDprime(x, y)

dp = (mean(x)-mean(y))./sqrt(mean([var(x) var(y)]));

function cname = ClusterFile(name,varargin)
filetype = 0;
prefix = [];
    j = 1;
    while j <= length(varargin);
        if isstruct(varargin{j})
            if isfield(varargin{j},'exptno')
                Ex = varargin{j};
                prefix = sprintf('Expt%d',floor(Ex.exptno));
                if Ex.blk > 0
                    prefix = [prefix 'a'];
                end
            end
        elseif strcmp(varargin{j},'auto')
            filetype = 1;
        elseif strcmp(varargin{j},'log')
            filetype = 2;
        end
        j = j+1;
    end
    if isdir(name)
        a = name;
    else
        [a,b] = fileparts(name);
    end
    if filetype == 2
        cname = [a '/AutoLog.txt'];
    elseif filetype == 1
        cname = [a '/' prefix 'AutoClusterTimes.mat'];
    else
        cname = [a '/' prefix 'ClusterTimes.mat'];
    end

        
    
function HistButtonPressed(src, data)
DATA = get(src,'UserData');

start = get(gca,'CurrentPoint');
DATA.elmousept.pos(1) = start(1,1);
DATA.elmousept.pos(3) = start(1,1);
yl = get(gca,'ylim')
if start(1,2) < yl(2) && start(1,2) > yl(1)
    DATA.elmousept.down = 1;
else
    return;
end
DATA.elmousept.axis = gca;
DATA.elmousept.plotargs = {};
DATA.elmousept.h = DrawEllipse(DATA.elmousept);
set(src,'UserData',DATA);

function HistButtonReleased(src, data)
DATA = get(src,'UserData');  % Just for this figure
start = get(gca,'CurrentPoint');
if isfield(DATA,'elmousept') && DATA.elmousept.down
DATA.elmousept.down = 0;
DATA.elmousept.done = 1;
p = DATA.elmousept.pos;
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) diff(p([1 3]))/2 diff(p([2 4]))/2]; 
set(src,'UserData',DATA);
DATA.cluster.crit = p(1);
E = BoundaryFromCluster(DATA.elmousept, DATA.cluster);

%now get main data 
DATA = GetDataFromFig(src);
[cl, DATA.cluster] = ClassifySpikes(DATA,E);
DATA.clid = cl.id;
DATA.nid = cl.nid;
DATA.clst = cl.clst;
DATA.MeanSpike = cl.MeanSpike;
set(DATA.toplevel,'UserData',DATA);
E.h = [];
DATA = ReplotPCs(DATA,E);
end


function HistButtonDragged(src, data)
DATA = get(src,'UserData');

if isfield(DATA,'elmousept') && DATA.elmousept.down > 0
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(1) = start(1,1);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
end
set(src,'UserData',DATA);


function h= oldDrawEllipse(E,varargin)

if E.shape == 1
    h = DrawLine(E,varargin{:});
    return;
end
a = (E.pos(3)-E.pos(1))/2; %x radius
b = (E.pos(4)-E.pos(2))/2;
sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)]+mean(E.pos([1 3]));
y = [y fliplr(-y) -y fliplr(y)]+mean(E.pos([2 4]));
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end

function h= oldDrawLine(E,varargin)

x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end


function h= DrawEllipse(E,varargin)

if E.shape == 1
    h = DrawLine(E,varargin{:});
    return;
end
a = (E.pos(3)-E.pos(1))/2; %x radius
b = (E.pos(4)-E.pos(2))/2;
sn = 0;
cn = 1;
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);

sn = sin(E.angle);
cn = cos(E.angle);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+mean(E.pos([1 3]));
y = yr+mean(E.pos([2 4]));
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y,'color',E.color);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),'color',E.color,varargin{:});
    hold off;
end

function h= DrawLine(E,varargin)

x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
if ishandle(E.h) 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end



function ButtonPressed(src, data)
DATA = GetDataFromFig(src);

start = get(gca,'CurrentPoint');
if InGraph(start,gca)
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})
%    DATA.elmousept.mode = mode;
    distance = DistanceToEllipse(DATA.elmousept,start(1,1:2))
    if  mode  == 2
        DATA.elmousept.down = 3;
        DATA.elmousept.done = 0;
        DATA.elmousept.start = start(1,1:2);
    elseif distance < 1.05 %test fails for NaN
        DATA.elmousept.down = 2;
        DATA.elmousept.done = 0;
        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
        DATA.elmousept.start = start(1,1:2);
        DATA.elmousept.axis = gca;
    else
        DATA.elmousept.down = 1;
        DATA.elmousept.done = 0;
        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
        DATA.elmousept.axis = gca;
    end
set(src,'UserData',DATA);
end

function distance = DistanceToEllipse(E, pos);
   
if isempty(E) | ~isfield(E,'pos');
    distance = NaN;
    return;
end
r(1) = (E.pos(3)-E.pos(1))/2; %x radius
r(2) = (E.pos(4)-E.pos(2))/2;
a(1) = (E.pos(3)+E.pos(1))/2; %x radius
a(2) = (E.pos(4)+E.pos(2))/2;

xy = pos - a;
xy = xy ./r;
cn = cos(-E.angle);
sn = sin(-E.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = sum(p.^2);


    function in = InGraph(pt, ax)
        xl = get(ax,'Xlim');
        yl = get(ax,'Ylim');
      in = pt(1,1) > xl(1) & pt(1,1) < xl(2) & pt(1,2) > yl(1) & pt(1,2) < yl(2);
      
function ButtonReleased(src, data)
DATA = GetDataFromFig(src);
if DATA.elmousept.down == 0 
    return;
end
mode = DATA.elmousept.down;
start = get(gca,'CurrentPoint');
DATA.elmousept.done = 1;
p = DATA.elmousept.pos;
DATA.elmousept.pcplot = get(gca,'UserData');
DATA.elmousept.down = 0;

if mode == 1
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
elseif mode == 2
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
end
set(src,'UserData',DATA);
if mode == 2  %touch inside ellispe to make cut
PCCluster(DATA,DATA.elmousept,1);
end

function ScrollWheel(src, evnt)
DATA = GetDataFromFig(src);
DATA.elmousept.angle = DATA.elmousept.angle+0.02*evnt.VerticalScrollCount;
fprintf('Angle %.2f\n',DATA.elmousept.angle);
DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
set(DATA.toplevel,'UserData',DATA);

 function PCButtonDragged(src, data)
DATA = GetDataFromFig(src);
if isfield(DATA,'elmousept')
if  DATA.elmousept.down == 1
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
elseif  DATA.elmousept.down == 3 %set radius with R mouse button
    start = get(gca,'CurrentPoint');
    r = abs(start(1,1:2) - DATA.elmousept.xyr(1:2));
    sina = sin(DATA.elmousept.angle);
    cosa = cos(DATA.elmousept.angle);
    r = r * [cosa sina; -sina cosa];
    DATA.elmousept.xyr([3 4]) = r;
    DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
elseif  DATA.elmousept.down == 2 %moving ellipse
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(1) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
end
set(src,'UserData',DATA);
end

function SetEllipseDrawing(DATA, shape,varargin)

cluster = 1;
plotargs = {};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cluster',4)
            j = j+1;
            cluster = varargin{j};
        else
            plotargs = {plotargs{:} varargin{j}};
        end
        j = j+1;
    end
    
    F = DATA.toplevel;

DATA.elmousept.h= -1;
DATA.elmousept.shape= shape;
DATA.elmousept.down = 0;
DATA.elmousept.done = 0;
DATA.elmousept.angle = 0;
DATA.elmousept.cluster = cluster;
DATA.elmousept.plotargs = plotargs;
DATA.elmousept.color = DATA.colors{cluster+1};
DATA.elmousept.dragfcn = get(F,'WindowButtonMotionFcn');
%should get old Fcns, then reset them after button release
set(F, 'WindowButtonDownFcn',@ButtonPressed);
set(F, 'WindowButtonMotionFcn',@PCButtonDragged);
set(F, 'WindowButtonUpFcn',@ButtonReleased);
set(F, 'WindowScrollWheelFcn',@ScrollWheel);
set(F,'UserData',DATA);
