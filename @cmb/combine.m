function varargout = combine(name, varargin)
%cmb.combine(file)
%cuts clusters and combines expts for Spike2 generated matlab data files
%
%

%To run combineall +me from commandline after firing up gui
%combine(F,'recombinename','rds.dxXmixacXdd','recombinemucells')
% can't get rid of cluster 2 (M019)
% show cluster colors in time display
% build voltage average for synchronous spikes;
% add checkbox for plotflip
% deleting clusters in AllClusters box
%
% Show spikes for non-stim trials only
%right button on a rotate ellipase sometimes cancels angle
% setting cluster 2 in a new space erates cluster 1
%scaling of PCA often mucked up
%artifact setting second regios sometimes undoes first, lemM102 p3 ex19
%autoscale density plot to peak of cluster?

TOPTAG = 'Combiner'; %default
LFP = [];
setexpts = [];
layout = [];
objver = ObjectVersion(cmb);

%disp(varargin);


nextcell = 0;
newprobe = 0;
DATA.bysuffix = 0;
initialconfig = '';

appendfiles = {};
argon = {};
for j = 1:nargout
    varargout{j} = [];
end;

if length(varargin) & isstruct(varargin{1})
    DATA = varargin{1};
    TOPTAG = DATA.tag.top;
else
    j = 1;
    while j <= length(varargin)
        if isstruct(varargin{j})
        elseif strcmp('Tag',varargin{j})
            TOPTAG = varargin{j+1};
        elseif strncmpi('setprobe',varargin{j},4)
            j = j+1;
            argon = {argon{:} 'setprobe' varargin{j}};
        elseif strncmpi('layout',varargin{j},6)
            j = j+1;
            DATA.layoutfile = varargin{j};
        elseif strncmpi(varargin{j},'nevdir',6)
            j = j+1;
            DATA.nevdir = varargin{j};
        elseif strncmpi('expts',varargin{j},4)
            j = j+1;
            setexpts = varargin{j};
            if exist('DATA','var') & isfield(DATA,'elst')
                set(DATA.elst,'value',setexpts);
                DATA.currentexpt = setexpts;
            end
        end
        j = j+1;
    end
end

if ishandle(name)
    DATA = GetDataFromFig(name);
    it = DATA.toplevel;
    if length(varargin) > 1 && isempty(varargin{1}) % call from menu
        name = varargin{2};
    end
else
    it = findobj('Tag',TOPTAG);
end

if isempty(it)
    DATA = cmb.Init(DATA, name, TOPTAG, layout);
    [DATA.spkvarnames, DATA.spkvarorder] = GetSpikeVals(DATA,NaN,NaN,NaN,[]);
    firstcall = 1;
    recombinename = [];
else
    firstcall = 0;
    if length(varargin) & isstruct(varargin{1})
        DATA = varargin{1};
        TOPTAG = DATA.tag.top;
    else
        DATA = get(it,'UserData');
        if isempty(strmatch(TOPTAG,'Combiner','exact'))
            DATA.forcetag = TOPTAG;
        else
            DATA.forcetag = '';
        end
    end
    if ~isfield(DATA,'toplevel') %%Error
        if strcmp(name,'quit')
            close(it);
        else
            fprintf('Cant Find combiune window - quitting. run cmb.combine(''quit'') to close window\n');
        end
        return;
    end
    set(DATA.toplevel,'Name','Busy......');
    drawnow;
end

if nargin > 0 & strncmpi(name,'nextcell',5)
    dir = splitpath(DATA.datafilename,'dir');
    if dir > 99
        newdir = num2str(str2num(dir)+1);
    else
        newdir = num2str(str2num(dir)+1,'%.3d');
    end
    name = strrep(DATA.datafilename,dir,newdir);
end

reindex = 0;
recombine = 0;
plotall = 0;
preprocess = 0;
plotallspikes = 0;
calltype = 'unknown';
if ischar(name)
    calltype = name;
end
if isfield(DATA,'toplevel') && DATA.state.verbose
    fprintf('Got combine data for %s\n',calltype);
end

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'Trials')
            LFP = varargin{j};
        end
    elseif strncmpi('append',varargin{j},6)
        j = j+1;
        appendfiles = {appendfiles{:} varargin{j}};
    elseif strncmpi(varargin{j},'autoinit',6)
        DATA.state.autocutatstart = 1;
    elseif strncmpi('bysuffix',varargin{j},6)
        DATA.bysuffix = 1;
    elseif strncmpi(varargin{j},'config',5)
        j = j+1;
        DATA.configfile = varargin{j};
        DATA = cmb.DoConfig(DATA,[],'load');
    elseif strncmpi('quicksuffix',varargin{j},9)
        DATA.bysuffix = 1;
        DATA.currentexpt = 1;
    elseif strncmpi(varargin{j},'expand',5)
        j = j+1;
        DATA.state.bigwindow(1) = varargin{j}(1);
        if length(varargin{j}) > 1
            DATA.state.bigwindow(2) = varargin{j}(1);
        end
    elseif strncmpi(varargin{j},'fontsize',6)
        j = j+1;
        DATA.gui.FontSize = varargin{j};
    elseif strncmpi(varargin{j},'fixdrift',6)
        DATA.state.fixdrift = 1;
    elseif strncmpi(varargin{j},'noninteractive',10)
        DATA.state.interactive = 0;
    elseif strncmpi(varargin{j},'newprobe',6)
        newprobe = 1;
        calltype = varargin{j};
    elseif strncmpi(varargin{j},'nospikes',6)
        newprobe = 1;
        argon = [argon varargin(j)];
    elseif strncmpi(varargin{j},'online',6)
        initialconfig = 'onlinespikes';
    elseif strncmpi(varargin{j},'preprocess',6)
        preprocess = 1;
    elseif strncmpi(varargin{j},'plotspikes',6)
        plotallspikes = 1;
    elseif strncmpi(varargin{j},'plotall',6)
        plotall = 1;
    elseif strncmpi(varargin{j},'pretty',4)
        DATA.plot.prettyfigs = 1;
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'recombinemucells',12)
        %do recombine at end = so dont set recombine > 0
        recombine = 0;
    elseif strncmpi(varargin{j},'recombinequit',12)
        recombine = 2;
    elseif strncmpi(varargin{j},'recombineorbw',12)
        recombine = 3;
    elseif strncmpi(varargin{j},'recombinename',12)
        recombine = 4;
        j = j+1;
        recombinename = varargin{j};
    elseif strncmpi(varargin{j},'recombine',5)
        recombine = 1;
    elseif strncmpi(varargin{j},'recount',5)
        DATA.state.recount = 1;
    elseif strncmpi(varargin{j},'relist',4)
        reindex = 1;
    elseif strncmpi(varargin{j},'showonline',5)
        DATA.state.showonlineclusters = 1;
    elseif strncmpi(varargin{j},'showspikes',6)
        DATA.state.showspikes= 1;
    elseif strncmpi(varargin{j},'spool',5)
        cmb.SpoolSpikes(DATA);
    elseif strncmpi(varargin{j},'usenev',5)
        DATA.state.usensx = 2;
    elseif strncmpi(varargin{j},'usensx',5)
        DATA.state.usensx = 1;
    elseif strncmpi(varargin{j},'voff',4)
        j = j+1;
        DATA.plot.voffsets = varargin{j};
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'vmax',4)
        j = j+1;
        DATA.plot.SpikeMaxV = varargin{j};
        set(DATA.toplevel,'UserData',DATA);
    end
    j = j+1;
end




if preprocess
    DATA = cmb.ReadFile(DATA, name);
    DATA = cmb.PrintComments(DATA);
    DATA = cmb.CalcOnlineShapes(DATA);
    return;
end
if firstcall
    DATA = cmb.BuildGUI(DATA);
    if length(setexpts)
        set(DATA.elst,'value',setexpts);
    end
    if ~isempty(initialconfig)
        DATA =cmb.OptionMenu(DATA,initialconfig,'QuickConfig');
    end
end

if iscell(name)  %given a cell arrary of Expts
    Expts = name;
    fname = Expts{1}.Header.expname;
    for j = 1:length(Expts)
        if ~isfield(Expts{j}.Stimvals,'st')
            Expts{j}.Stimvals.st = 2;
        end
        if ~isfield(Expts{j}.Header,'psych')
            if isfield(Expts{j}.Trials,'RespDir')
                Expts{j}.Header.psych = 1;
            else
                Expts{j}.Header.psych = 0;
            end
        end
        Expts{j}.Header.Name = Expts{j}.Header.expname;
        Expts{j}.Header.expname = Expt2Name(Expts{j});
        Expts{j}.gui.clustertype = 5;
        Expts{j}.gui.ncluster = 0;
        Expts{j}.gui.classified = 5;
        Expts{j}.gui.counted = 5;
        Expts{j}.Comments.text = [];
        Expts{j}.Header.trange(1) = Expts{j}.Trials(1).Start(1);
        Expts{j}.Header.trange(2) = Expts{j}.Trials(1).End(end);
    end
    DATA.Expts = Expts;
    DATA.state.plotseq = 5;
    DATA.state.psychonly = 1;
    DATA.state.plotpsych = 1;
    DATA.state.online = 0;
    DATA = cmb.cListExpts(DATA,Expts);
    set(DATA.toplevel,'UserData',DATA);
elseif isfigure(name)
elseif strncmpi(name,'autocut',6)
    eid = get(DATA.elst,'value');
    for j = 1:length(eid)
        cmb.AutoCut(DATA,DATA.expid(eid(j)),eid(j));
    end
elseif strncmpi(name,'autocluster',6)
    eid = get(DATA.elst,'value');
    args = {};
    if DATA.state.optimizeclusters
        args = {args{:} 'minimize'};
    end
    if sum(DATA.plot.useprobe)
        [nr,nc] = Nsubplots(length(DATA.probelist));
        oldprobe = DATA.probe;
        pid = find(DATA.plot.useprobe);
        DATA = cmb.GetAllProbeFig(DATA);
        for j = 1:length(eid)
            for k = pid
                DATA.probe = DATA.probelist(k);
                subplot(nr,nc,k);
                DATA = cmb.AutoCut(DATA,DATA.expid(eid(j)),eid(j),'cluster',args{:});
            end
        end
        DATA.probe = oldprobe;
    else
        for j = 1:length(eid)
            DATA = cmb.AutoCut(DATA,DATA.expid(eid(j)),eid(j),'cluster',args{:});
        end
    end
    DATA.AutoClusters = cmb.SaveAutoClusters(DATA);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'AddProbeList',5)
    cmb.AddProbeList(DATA);
    DATA.plot.showsync = 1;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'mucombine',5)
    DATA.listbycell = 2;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'fixmains',6)
    load(DATA.datafilename);
    tic;
    if isfield(DATA,'AllSpikes')
        for j = 1:length(DATA.AllSpikes)
            [avg(:,j), DATA.AllSPikes{j}] = FixSpkMains(DATA.AllSpikes{j},Ch24.times.*10000,'fix');
        end
        set(DATA.toplevel,'UserData',DATA);
    else
        avg = FixSpkMains(DATA.AllData.Spikes,Ch24.times.*10000);
    end
    toc
    GetFigure('Mains');
    plot(avg);
elseif strncmpi(name,'plotddf',6)
    eid = get(DATA.elst,'value');
    cmb.PlotDDF(DATA);
elseif strncmpi(name,'readlog',6)
    fid = fopen(strrep(DATA.datafilename,'.mat', '.log'),'r');
    a = textscan(fid, '%s','delimiter','\n');
    txt = char(a{1});
    sid = strmatch('Save',txt(:,22:end));
    nid = regexp(a{1}(sid),'c[0-9]\..*\.mat');
    for j = 1:length(sid)
        s = a{1}{sid(j)};
        enames{j} = a{1}{sid(j)}(nid{j}+3:end-4);
        id = strfind(s,'Expts');
        expids{j} = sscanf(s(id+6:end),'%d');
        dates(j) = datenum(s(1:20));
        if strfind(enames{j},'image.ORBW')
            types(j) = 1;
        elseif strfind(enames{j},'ORBW')
            types(j) = 2;
        elseif strfind(enames{j},'image.DCORRC')
            types(j) = 3;
        elseif strfind(enames{j},'DCORRC')
            types(j) = 4;
        elseif strfind(enames{j},'cylinder.ABD')
            types(j) = 5;
        else
            types(j) = 0;
        end
    end
    DATA.log.expnames = enames;
    DATA.log.savedates = dates;
    DATA.log.savexps = expids;
    DATA.log.types = types;
    cmb.PrintLogData(DATA,1);
    cmb.PrintLogData(DATA,2);
    cmb.PrintLogData(DATA,5);
    
elseif strncmpi(name,'relistexps',8)
    DATA = cmb.cListExpts(DATA, DATA.Expts);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'testcounts',8)
    spikelist = cmb.WhichClusters(DATA.toplevel);
    ts = now;
    for j = 1:length(DATA.Expts)
        DATA = cmb.CountSpikes(DATA,j);
    end
    fprintf('CountSpikes Took %.2f\n',mytoc(ts));
    ts = now;
    for j = 1:length(DATA.Expts)
        DATA = cmb.CountSpikesB(DATA,j,DATA.probe,spikelist);
    end
    fprintf('CountSpikes Took %.2f\n',mytoc(ts));
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'listexps',6)
    eid = get(DATA.clst,'value');
    combineids = [];
    if eid(1) > 1
        suff = '';
        it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
        if isempty(it)
            it = strmatch(regexprep(DATA.exptypelist{eid(1)},'RC$',''),DATA.expstrs,'exact');
        end
        stimname = strrep(DATA.explist{eid(1)},['.' DATA.exptypelist{eid(1)}],'');
        if strfind(stimname,'CRC')
            stimname = strrep(stimname,'CRC','');
            suff = 'CRC';
        elseif strfind(stimname,'RC')
            stimname = strrep(stimname,'RC','');
            suff = 'RC';
        end
        if ~isempty(it)
            expname = DATA.expnames{it};
        else
            expname = DATA.exptypelist{eid(1)};
        end
        if 1 %was test for it
            if strmatch(expname,'CO','exact') & strmatch(stimname,'square')
                expname = 'FLSH';
            end
            if regexp(DATA.datafilename,'Expt[0-9]*.mat') %online file
                DATA.outname = [DATA.datafilename '/c' num2str(DATA.spikelist(1)) '.' stimname '.' DATA.expnames{it} suff '.mat'];
            elseif DATA.listbycell
                DATA.outname = cmb.CombinedName(DATA,eid(1),DATA.spikelist(1),'cell',DATA.probe);
            else
                DATA.outname = cmb.CombinedName(DATA,eid(1),DATA.spikelist(1));
            end
            expfile = DATA.outname;
            if DATA.listbycell
                expfile = regexprep(DATA.outname,'\.p[0-9]*c[0-9]\.',sprintf('.cell%d.',DATA.probe));
                args = {};
            else
                args = {'probe'};
            end
            if ~exist(expfile,'file')
                expfile = regexprep(DATA.outname,'\.p[0-9]*c[0-9]\.','.Cells.');
            end
            if exist(expfile,'file')
                cmb.SetFigure(DATA.tag.dataplot,DATA);
                load(expfile);
                if exist('AllExpt','var')
                    Expt = All2Expt(AllExpt,DATA.probe,args{:});
                end
                Expt.Header.Name = cmb.BuildName(Expt.Header.Name);
                DATA.Expt = Expt;
                args = cmb.PlotArgs(DATA, Expt,'combined');
                DATA.plotres = PlotExpt(Expt,args{:});
                csuffs = [];
                if isfield(Expt.Header,'Combineids')
                    combineids = Expt.Header.Combineids;
                    for j = 1:length(Expt.Header.Combineids)
                        if isfield(Expt.Header,'suffixes') && sum(Expt.Header.suffixes >0) == length(Expt.Header.Combineids)
                            si = Expt.Header.suffixes(j);
                            ci = find(ismember(DATA.suffixlist,si)); 
                        elseif isfield(DATA,'suffixlist')
                            ci = Expt.Header.Combineids(j);
                            ci = DATA.suffixlist(ci);
                        else
                            ci = Expt.Header.Combineids(j);
                        end
                        
                        csuffs = [csuffs num2str(ci) ' '];
                        if ci > length(DATA.Expts)
                            fprintf('Combined file refers to too many expts\n');
                        elseif ci ==0
                            fprintf('Invalid combine expt list\n');
                        elseif isfield(DATA,'AllClusters') %Don't Add Expt.Header.Cluster
                        elseif isfield(Expt.Header,'Clusters') & length(Expt.Header.Clusters) >= j & ...
                                (~isfield(DATA.Expts{ci},'Cluster') || isempty(DATA.Expts{ci}.Cluster))
                            DATA.Expts{ci}.Cluster = Expt.Header.Clusters{j};
                        end
                    end
                    title(sprintf('%s Ex %s ed %.2f',splitpath(DATA.outname),csuffs,GetEval(Expt,'ed')));
                else
                    for j = 1:length(Expt.Header.Combined)
                        csuffs = [csuffs num2str(Expt.Header.Combined(j)) ' '];
                    end
                    title(sprintf('%s %s ed %.2f',splitpath(DATA.outname),csuffs,GetEval(Expt,'ed')));
                end
            end
        else
            DATA.outname = [strrep(DATA.datafilename,'.mat','.') '.c' num2str(DATA.spikelist(1)) '.' stimname '.' DATA.exptypelist{eid(1)}  suff '.mat'];
            DATA.outname = cmb.CombinedName(DATA,eid(1),DATA.spikelist(1));
        end
    else
        DATA.outname = 'tmp.mat';
    end
    DATA = cmb.ListSubExpts(DATA,eid);
    if combineids
        id = find(ismember(DATA.expid,combineids));
        if length(id)
            set(DATA.elst,'value',id);
        end
    end
    set(DATA.saveitem,'string',DATA.outname);
    set(DATA.toplevel,'UserData',DATA);
    if nargout
        varargout{1} = DATA;
    end
elseif strncmpi(name,'getexpt',6)
    varargout{1} = DATA.Expt;
    cmb.NotBusy(DATA);
    return;
elseif strncmpi(name,'getstate',6)
    varargout{1} = DATA;
    cmb.NotBusy(DATA);
    return;
elseif strncmpi(name,'ingnoreonline',9)
    DATA.state.useonlineclusters = 0;
elseif strncmpi(name,'checklists',6)
    cmb.CheckLists(DATA);
elseif strncmpi(name,'ClearClusters',7)
    eid = DATA.currentexpt(1);
    nc = cmb.CountClusters(DATA.Expts{eid}.Cluster);
    for k = 1:nc;
        for j = 1:size(DATA.Expts{eid}.Cluster,2);
            if ~isempty(DATA.Expts{eid}.Cluster{k,j}) & isfield(DATA.Expts{eid}.Cluster{k,j},'x')
                DATA.Expts{eid}.Cluster{k,j}.touched = 0;
                DATA.Expts{eid}.Cluster{k,j} = rmfield(DATA.Expts{eid}.Cluster{k,j},{'x' 'y'});
            end
        end
    end
    nc = cmb.CountClusters(DATA.cluster);
    for k = 1:nc;
        for j = 1:size(DATA.cluster,2);
            if ~isempty(DATA.cluster{k,j}) & isfield(DATA.cluster{k,j},'x')
                DATA.cluster{k,j}.touched = 0;
                DATA.cluster{k,j} = rmfield(DATA.cluster{k,j},{'x' 'y'});
            end
        end
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'ClearallClusters',10)
    p = DATA.probe;
    for eid = 1:length(DATA.Expts)
        nc = cmb.CountClusters(DATA.Expts{eid}.Cluster);
        for k = 1:nc;
            if ~isempty(DATA.Expts{eid}.Cluster{k,p}) & isfield(DATA.Expts{eid}.Cluster{k,p},'x')
                DATA.Expts{eid}.Cluster{k,p}.touched = 0;
                DATA.Expts{eid}.Cluster{k,p} = rmfield(DATA.Expts{eid}.Cluster{k,p},{'x' 'y'});
            end
        end
    end
    nc = cmb.CountClusters(DATA.cluster);
    for k = 1:nc;
        for j = 1:size(DATA.cluster,2);
            if ~isempty(DATA.cluster{k,p}) & isfield(DATA.cluster{k,p},'x')
                DATA.cluster{k,p}.touched = 0;
                DATA.cluster{k,p} = rmfield(DATA.cluster{k,p},{'x' 'y'});
            end
        end
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'comments',6)
    cmb.PrintComments(DATA);
    DATA.currentprobe = DATA.probe;
    DATA.Comments = PlotComments(DATA.datadir,'parent',DATA.toplevel);
elseif strncmpi(name,'combine',6)
    DATA.extype = get(DATA.clst,'value');
    if length(DATA.extype) > 1
        DATA.extype = DATA.extype(1);
    end
    [Expt, DATA] = cmb.CombinePlot(DATA,1);
    DATA.Expt = Expt;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'options',5)
    DATA.optionfig = cmb.setoptions(DATA,DATA.tag.options);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'loadclusters',6)
    DATA = cmb.LoadClusters(DATA,cmb.ClusterFile(DATA));
    args = {args{:} varargin{j}};
elseif strncmpi(name,'loadsaccades',6)
    for j = 1:length(DATA.Expts)
        E = LoadEmData(DATA.Expts{j});
        if isfield(E.Trials,'Saccades')
            for k = 1:length(E.Trials)
                DATA.Expts{j}.Trials(k).Saccades = E.Trials(k).Saccades;
            end
        end
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'meanspike',6)
    ts = now;
    if length(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
        DATA = cmb.CalcMeanSpike(DATA,varargin{1});
    else
        DATA = cmb.CalcMeanSpike(DATA,DATA.expid);
    end
    mytoc(ts);
    set(DATA.toplevel,'UserData',DATA);
    if strfind(name,'save')
        cmb.SaveSpikeShape(DATA,DATA.meanspkfile);
    end
elseif strncmpi(name,'plotshape',6)
    GetFigure('ClusterShape');
    PlotSpikeShapes(DATA.MeanSpike,varargin{:});
elseif strncmpi(name,'PlotISI',5)
    GetFigure('ISI');
    [isis, t, s] = CalcISI(DATA.Expts{DATA.currentexpt(1)}.Trials);
    id = find(isis < 1000)
    hist(isis(id),100);
    [a,b] = sort(isis);
    %    DATA.isit = t(b);
    DATA.isis = s(b);
    DATA.ISIpair = 1;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'reclassify',5)
    ts = now;
    DATA = cmb.ReClassifyAll(DATA);
    fprintf('Took %.2f sec\n',(now-ts)*(24*60*60));
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'oclassify',5) %apply online spike classification and calculate mean spike
    DATA = CalcOnlineShaptes(DATA);
elseif strncmpi(name,'revertcelllist',10)
    load(DATA.cellfile);
    DATA.CellList = CellList;
    DATA.CellListCluster = CellListCluster;
    DATA.CellQuality = CellQuality;
    if length(DATA.CellQuality) < length(DATA.CellList)
        lq = length(DATA.CellQuality);
        ll = length(DATA.CellList);
        DATA.CellQuality(:,lq+1:ll) = 0;
    end
    if exist('Templates','var')
        DATA.Templates = Templates;
    end
    if exist('TemplateInfo','var') && ~isempty(DATA.TemplateInfo)
        DATA.TemplateInfo = TemplateInfo;
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'loadconfig',6)
    [outname, pathname] = uigetfile(DATA.configfile);
    if outname
        DATA.configfile = [pathaname outname];
        cmb.DoConfig(DATA,DATA.configfile);
        set(DATA.toplevel,'UserData',DATA);
    end
elseif strncmpi(name,'readconfig',6) %call to readconfig after starting
    DATA  = ReadConfig(DATA);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'showvals',5)
    DATA.showid = cmb.setshow(DATA,DATA.tag.showvals);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'spikes',5)
    cmb.PlotSpikes(DATA,varargin{1});
elseif strncmpi(name,'setexp',6)

    if DATA.playingspk & ishandle(DATA.svfig)
        set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
        set(findobj(DATA.svfig,'Tag','StopSpool'),'UserData',1);
        DATA.playingspk = 0;
        set(DATA.toplevel,'UserData',DATA);
        return;
    end
    if newprobe && DATA.listbycell
        eid = get(DATA.clst,'value');
        DATA = cmb.ListSubExpts(DATA,eid);
        
        
        cmb.SetFigure(DATA.tag.dataplot,DATA);
        cellexpname = regexprep(DATA.outname,'\.[p,0-9]*c[0-9*]\.',sprintf('.cell%d.',DATA.probe));
        if exist(cellexpname,'file')
            load(cellexpname);
            Expt.Header.Name = cmb.BuildName(Expt.Header.Name);
            DATA.Expt = Expt;
            args = cmb.PlotArgs(DATA, Expt,'combined');
            DATA.plotres = PlotExpt(Expt,args{:});
        end
        
        if isempty(DATA.expid)
            varargout{1} = DATA;
            set(DATA.toplevel,'UserData',DATA);
            return;
        end
    end
    DATA = cmb.CheckState(DATA);
    ts  = now;
    id = get(DATA.elst,'value');
    DATA.exabsid = DATA.expid(id);
    cmb.SetFigure(DATA.tag.dataplot,DATA,'noforce');
    if DATA.state.verbose > 1
        fprintf('Created Figure\n');
    end
    plotted = 0;
    if DATA.state.plotcombined || length(DATA.exabsid) == 1
        rc = DATA.plot.condenseRC;
        DATA.plot.condenseRC = 0;
        % why was this plotting the combined data? Want to plot individual ones
        % but using the combined plott style (legends, colors, so that can see what
        % is goin on in individual exps
        %        cmb.PlotCombined(DATA, DATA.Expt);
        for j   = 1:length(DATA.exabsid)
            if sum([DATA.Expts{DATA.exabsid(j)}.Trials.count]) > 0
                if j > 1
                    cmb.PlotCombined(DATA, DATA.Expts{DATA.exabsid(j)},'hold');
                else
                    cmb.PlotCombined(DATA, DATA.Expts{DATA.exabsid(j)});
                end
                plotted = plotted+1;
            else
                fprintf('No Spikes in E%dP%d\n',DATA.exabsid(j),DATA.probe);
            end
        end
        DATA.plot.condenseRC = rc;
        cmb.NotBusy(DATA);
        if nargout
            varargout{1} = DATA;
        end
        if 0 %need next steps when using plot as combined one expt at at time
            %use seomthing else to save time here.
        return;
        end
    end
    
    playargs = {};
    hs = 'nohold';
    colors = mycolors;
    if strncmpi(name,'setexpplot',8)
        playspk = 0;
    else
        playspk = cmb.GetCheck('ShowSpikes',DATA);
        DATA.state.showspikes = playspk;
    end
    if length(varargin) && strncmpi(varargin{1},'flagchange',6)
        playspk = 0; %don't spool spikes every time a flag is changed
    end
    if DATA.state.online && isfield(DATA,'timerobj')
        ton = strcmp(get(DATA.timerobj,'Running'),'on');
        if DATA.expid(id(end)) == length(DATA.Expts) && DATA.state.autolist
            if ~ton
                start(DATA.timerobj);
                DATA = get(DATA.toplevel,'UserData');
            end
        elseif ton
            stop(DATA.timerobj);
        end
    end
    lastexpt = DATA.currentexpt(1);
    if isempty(id)
        DATA.currentexpt = DATA.expid(1);
    else
        DATA.currentexpt = DATA.expid(id);
    end
    DATA.currenttrial = 0;
    DATA.firsttrial = 0; %no partial list
    DATA.lasttrial = 0;
    cmb.CheckSpoolButton(DATA);
    if ~isempty(findobj('Tag',DATA.tag.celllist))  && DATA.state.autoplotcells
        GetFigure(DATA.tag.celllist);
        cmb.PlotCellList(DATA);
    end

    if strncmp(DATA.filetype,'Grid',4) && playspk
        if DATA.state.usensx && DATA.state.online == 0
            DATA = cmb.ReadGridFile(DATA);
            DATA.plot.useprobe = zeros(size(DATA.probelist));
        else
            if ~isfield(DATA,'AllClusters')
                DATA.AllClusters = {};
            end
            if ~isfield(DATA,'grididx')
                DATA.grididx = BuildGridIndex(DATA.nevdir, DATA.Expts);
            elseif max(DATA.currentexpt) > max(DATA.grididx.expt) %% need to reindex
                DATA.grididx = BuildGridIndex(DATA.nevdir, DATA.Expts,'reindex','noerrs');
            end
            
            DATA = cmb.LoadSpikes(DATA, DATA.currentexpt(1));
            DATA = cmb.CheckClusterLoaded(DATA, DATA.currentexpt(1), DATA.probe);
        end
        playargs = [playargs 'Onetrial'];
        %        DATA.plot.useprobe(DATA.probe) = 1;
    elseif DATA.state.online & lastexpt ~= DATA.currentexpt(1)
%DATA.currentexpt(1) is teh key because this is used for spooling spikes, showing XY        
        if DATA.state.usexycache %if showing spikes, to load
            [DATA, ok] = cmb.SpkCache(DATA,DATA.currentexpt(1),DATA.probe,'add');
        else
            ok = 0;
        end
        if ~ok
        DATA = cmb.LoadSpikes(DATA, DATA.currentexpt(1));
        end
        if DATA.plot.quickspks
            cmb.PlaySpikes(DATA,DATA.currentexpt(1),'quickspks');
        end
    elseif isfield(DATA,'AllClusters') && iscell(DATA.AllClusters)
        DATA = cmb.SetProbe(DATA, DATA.probe);
        if DATA.listbycell && newprobe
            DATA =  cmb.ListSubExpts(DATA,0);
        end
        %        for j = 1:length(DATA.expid)
        %                [DATA, counts{j}] = cmb.CountSpikes(DATA,DATA.expid(j)); %% recount in case cluster # changed
        %        end
    elseif DATA.plot.quickspks
        DATA = cmb.LoadSpikes(DATA, DATA.currentexpt(1));
        DATA = cmb.SetCluster(DATA, DATA.currentexpt(1), DATA.probe);
        DATA = cmb.PlaySpikes(DATA,DATA.currentexpt(1),'quickspks');
    end
    
    % if not classified, and autospool is set, show spikes. Unless spkXY is set, in whcih case just
    % show XY plot
    if isempty(DATA.expid)
        varargout{1} = DATA;
        set(DATA.toplevel,'UserData',DATA);
        return;
    elseif DATA.state.online < 2 & DATA.Expts{DATA.expid(id(1))}.gui.classified == 0 && ...
            ~DATA.state.showspkxy && DATA.state.autospool || ...
            (DATA.Expts{DATA.expid(id(1))}.gui.classified == 2 && DATA.state.autospool)
        playspk = 1;
    end
    if DATA.state.showspkxy  %trumps other options
        playspk = 0;
    end
    if DATA.state.autoplotallspikes == 0 && length(id) > 1
        playspk = 0;
    end
    if DATA.state.nospikes && DATA.state.somespikes == 0 %override all the possibilites
        playspk = 0;
    end
    if isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
        DATA.allexp = DATA.currentexpt(1);
    end
    % need to copy Expt.Cluster to DATA.cluster for all probes
    if DATA.state.resetclusters
        for j = 1:size(DATA.cluster,2)
            for k = 1:size(DATA.cluster,1)
                if isfield(DATA.cluster{k,j},'x')
                    DATA.cluster{k,j} = rmfield(DATA.cluster{k,j},{'x' 'y'});
                end
            end
        end
    end
    
    
   %when a new expt is loaded via 'setexp'
   %this checks all probes to see if Cluster is defined, and copies
   %current cluster to all that are absent.
    ndef = 0;
    probe = DATA.probe;
    method = 1;
    if method ==1
        if newprobe
            DATA = cmb.SetCluster(DATA, DATA.currentexpt(1), DATA.probe);
        else
            DATA = cmb.SetCluster(DATA, DATA.currentexpt(1), 1:size(DATA.cluster,2));
        end
    else
    for j = 1:size(DATA.cluster,2)
        pdef(j) = 0;
        ndef(j) = 0;
        for k = 1:size(DATA.cluster,1)
            if isfield(DATA.Expts{DATA.currentexpt(1)},'Cluster') && cmb.iscluster(DATA.Expts{DATA.currentexpt(1)}.Cluster,1,j)
                DATA.cluster{k,j} = DATA.Expts{DATA.currentexpt(1)}.Cluster{1,j};
                DATA.cluster{k,j}.touched = 1;
                ndef(j) = ndef(j)+1;
                pdef(j) = pdef(j)+1;
            else
                DATA.cluster{k,j}.touched = 0;
            end
        end
        if ndef(j) == 0 && cmb.iscluster(DATA.cluster,1,j) & DATA.state.applylastcluster
            for k = 1:size(DATA.cluster,1)
                if cmb.iscluster(DATA.cluster,k,j)
                    DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j} = DATA.cluster{k,j};
                    if isfield(DATA.cluster{k,j},'exptno')
                        DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j}.fromlast = DATA.cluster{k,j}.exptno;
                    else
                        DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j}.fromlast =1;
                    end
                    ndef(j) = ndef(j)+1;
                end
            end
        end
    end
    end
    if DATA.state.classifyallprobes & isfield(DATA,'AllSpikes')
        for j = 1:size(DATA.cluster,2)
            if ndef(j) & ~isempty(DATA.AllSpikes{j})
                DATA.probe = j;
                DATA = SetExptSpikes(DATA,DATA.currentexpt(1),0,'useexp');
                %                DATA.Expts{DATA.currentexpt(1)}.Cluster{1,j}.touched = 1;
            end
        end
        for j = 1:length(DATA.probelist)
            if isfield(DATA.AllSpikes{j},'codes')
                ns(j) = sum(DATA.AllSpikes{j}.codes(:,2) > 0);
            else
                ns(j) = NaN;
            end
        end
        DATA.probe =probe;
        
    end
    
    DATA.nclusters = cmb.CountClusters(DATA.cluster);
    
    
    
    ei = DATA.expid(id(1));
    if playspk & DATA.probe == 100
        DATA = cmb.PlaySpikes(DATA, DATA.expid(id(1)));
    elseif playspk & (~isempty(DATA.AllData.Spikes) || isfield(DATA,'AllSpikes'))
        if DATA.spooling %was looking at a different classification. Re-classif
            for j = id(1:end)
            end
        end
        DATA.spklist = [];
        DATA.spooling = 0;
        
        for j = id(1:end)
            % If this expt last shown with a temporary cut, revert back to the saved
            % one, unless this is a "spool" showing the temporary cut.
            % or if the classifiction has  not been done at all yet.
            if isfield(DATA,'AllSpikes') && DATA.plot.UsePreviousCut
                oldprobe = DATA.probe;
                for k = 1:length(DATA.plot.useprobe)
                    if DATA.plot.useprobe(k)
                        DATA.probe = k;
                        DATA = SetExptSpikes(DATA,DATA.expid(id(1)),0);
                    end
                end
                DATA.probe = oldprobe;
            elseif ~DATA.spooling && ismember(DATA.Expts{DATA.expid(j)}.gui.classified,[0 2])
                DATA = SetExptSpikes(DATA,DATA.expid(id(1)),0);
            else
                [DATA.plot.clusterX , DATA.plot.clusterY, DATA] = cmb.GetClusterSpace(DATA, DATA.Expts{DATA.expid(id(1))});
            end
            DATA = cmb.PlaySpikes(DATA, DATA.expid(j),playargs{:});
        end
        if DATA.state.showspkxy & isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
            DATA =  cmb.PlotAllProbeXY(DATA,0);
        end
        set(DATA.toplevel,'UserData',DATA);
        if ~DATA.state.autoreplotgraph %only show spikes
            cmb.NotBusy(DATA);
            varargout{1} = DATA;
            return; % This stops plotting of results. Why? ? only return if not auto
        end
    elseif DATA.state.showspkxy
        %TODO set nclusters, ed properly. And set recut to 2 if necessary.
        %            or maybe call SetExptSpks here....
        
        nclusters = 2;
        if isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
            if DATA.plot.showallxy
                DATA =  cmb.PlotAllProbeXY(DATA,0);
            end
            [DATA, ispk] = SetExptSpikes(DATA,DATA.expid(id(1)),0,'useexp');
            DATA.spklist = ispk;
            DATA.Expts{DATA.currentexpt(1)}.gui.spks = ispk;
            if DATA.plot.quickspks && DATA.state.showspikes
                cmb.PlaySpikes(DATA,DATA.currentexpt(1),'quickspks');
            else
                cmb.PlaySpikes(DATA,DATA.currentexpt(1),'xyonly');
            end
            
            cmb.SetFigure(DATA.tag.clusterxy,DATA);
            hold off;
            DATA = cmb.DrawXYPlot(DATA, ispk);
            ispk = [];
        else
            if ~isfigure(DATA.svfig) %dont let this spool everything just to create figure
                DATA = cmb.PlaySpikes(DATA,ei,'quickspks');
            end
            cmb.SetSetBox(DATA,ei);
            %           GetFigure(DATA.tag.clusterxy);
            %if use 'useexpall' then no clusters carry over, even if it has
            %neve been touched
            [DATA, ispk] = SetExptSpikes(DATA,DATA.expid(id(1)),0,'useexp');
            DATA.spklist = ispk;
            DATA = CalcClusterVars(DATA, DATA.Expts{DATA.currentexpt(1)}.gui.spks);
            cmb.SetFigure(DATA.tag.clusterxy,DATA);
            hold off;
            DATA = cmb.DrawXYPlot(DATA, ispk);
            nclusters = max(DATA.AllData.Spikes.codes(ispk,2));
        end
        stimes(1) = DATA.Expts{ei}.Trials(end).Start(1);
        stimes(2) = DATA.Expts{ei}.Trials(end).End(end);
        %DATA.currentexpt = ei;
%        tspk = FindSpikes(DATA, stimes, DATA.probe,ei);
%        cmb.DrawSpikeWaves(DATA, tspk, nclusters, 2);
    elseif isfield(DATA.Expts{ei},'Cluster')
        if isfigure(DATA.xyfig)
            set(0,'CurrentFigure',DATA.xyfig);
            hold off;
            DATA.cluster = DATA.Expts{ei}.Cluster;
            cmb.DrawClusters(DATA,DATA.Expts{ei}.Cluster, 0);
        end
    end
    %  clearing spklist causes the problem with online files where the spik
    %  display changes afterthe cluster is cut. If Spklist is empty, then when
    %  the cluster is cut, all spikes (including those not in trials) are
    %  shown, until the spool button is set. What goes wrong if this is not
    %  cleared?
    %    DATA.spklist = []; %once clusters are cut, this should be in them. Don't let first expt set for all
    
    if cmb.isExptCluster(DATA.Expts{ei},DATA.currentcluster,DATA.probe) & ...
            isfield(DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe},'quality');
        a = DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe}.quality;
        if(a > 0) %if its not set in cluster, let it inherit value from last expt
            set(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value',a+1);
        end
    elseif isfield(DATA.Expts{ei},'Cluster') && iscell(DATA.Expts{ei}.Cluster)
        DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe}.quality = 0;
    end
%when online, always plot Expt when probe change brings us here    
%unless we have already plotted becuase of DATA.state,plotcombined
    if (newprobe && DATA.state.somespikes == 0 && DATA.state.online ~= 1) ||  (DATA.state.online == 1 && plotted > 0)
        SetData(DATA);
        varargout{1} = DATA;
        return;
    end
    if newprobe && length(id) == 1
        DATA.Expts{DATA.expid(j)}.gui.counted = 0;
    end
    cmb.SetFigure(DATA.tag.dataplot,DATA,'noforce');
    ClearPlot;
    h = [];
    xargs = {};
    if length(id) > 1 && DATA.Expts{DATA.expid(id(1))}.Header.rc
        xargs = {xargs{:} 'condense'};
    end
    eid = get(DATA.clst,'value');
    DATA.spikelist = cmb.WhichClusters(DATA.toplevel);
    
    for j = id(1:end)
        tc=now;
        fprintf('Counting Spikes for expt %d (%.2f)\n',j,mytoc(ts));
        %        eid = DATA.expid(j); %eid already used for DATA.clst value
        p = GetProbe(DATA, DATA.expid(j), DATA.probe);
        if DATA.state.recount | DATA.Expts{DATA.expid(j)}.gui.counted == 0 & eid(1) > 1
            DATA = cmb.CountSpikes(DATA,DATA.expid(j));
        end
        
        s = cmb.ExptComments(DATA.Expts{DATA.expid(j)});
        
        te = now;
        edepth = GetEval(DATA.Expts{DATA.expid(j)},'ed');
        if length(id)  == 1 || eid(1) > 1 %don't add plots from 'All' list
            args = cmb.PlotArgs(DATA,DATA.Expts{DATA.expid(j)});
            if j > id(1)
                args = {args{:} 'hold'};
            end
            cmb.SetFigure(DATA.tag.dataplot,DATA);
            res = PlotExpt(DATA.Expts{DATA.expid(j)},hs,'forcecolor',colors{j},args{:},'legendpos',7,xargs{:});
            fprintf('%.2f plot (%.2f count)) ',(now-te)*60*60*24,(te-tc)*60*60*24);
            sp  = strfind(DATA.subexplist{j},' ');
            if sp
                names{j} = DATA.subexplist{j}(1:sp(1)-1);
            else
                names{j} = DATA.subexplist{j};
            end
            names{j} = [names{j} cmb.ShowString(DATA,DATA.Expts{DATA.expid(j)})];
            if  ~isempty(res) && isfield(res,'handles') && ishandle(res(1).handles(1))
                h(j) = res(1).handles(1);
            end
            hs = 'Hold';
        end
    end
    if ~isempty(h)
        mylegend(h,names);
    end
    t = get(get(gca,'title'),'String');
    DATA.spikelist = cmb.WhichClusters(DATA.toplevel);
    title([t sprintf(' P%d',DATA.probe) 'Cl' sprintf(' %d',cmb.WhichClusters(DATA.toplevel)) sprintf('ed%.2f',edepth)]);
    %    if DATA.state.verbose
    if 1
        fprintf('Took %.2f\n',(now-ts)*60*60*24);
    end
    varargout{1} = DATA;
    set(DATA.toplevel,'UserData',DATA);
    
elseif strncmpi(name,'Close',4)
    cmb.CloseCombine(DATA);
    return;
elseif strncmpi(name,'postperiod',8)
    j = j+1;
    DATA.state.postperiod = varargin{1};
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'newoptions',8)
    newoptions(DATA,DATA.tag.options);
elseif strncmpi(name,'preperiod',8)
    j = j+1;
    DATA.state.preperiod = varargin{1};
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'profile',6)
    DATA.profiling = ~DATA.profiling;
    fprintf('Profiling %d\n',DATA.profiling);
    set(DATA.toplevel,'UserData',DATA);
    
elseif strncmpi(name,'showdprime',9)
    DATA.plot.showdprime = 1;
    set(DATA.toplevel,'UserData',DATA);
    fprintf('Showing additional dprime plots\n');
    dp = cmb.CalcDprime(DATA,1);
    fprintf('Dp = %.3f\n',dp);
    dp = cmb.CalcDprime(DATA,1,'yparam',8); %spk var
    fprintf('Dp = %.3f\n',dp);
    dp = cmb.CalcDprime(DATA,1,'yparam',44); %spk var/sqrt(energy)
    fprintf('Dp = %.3f\n',dp);
    dp = cmb.CalcDprime(DATA,1,'yparam',45); %spk var/sqrt(energy)
    fprintf('Dp = %.3f\n',dp);
    dp = cmb.CalcDprime(DATA,1,'xparam',46,'yparam',44); %spk var/sqrt(energy)
    fprintf('Dp = %.3f\n',dp);
elseif strncmpi(name,'savelfp',7) || strncmpi(name,'makelfp',7)
    outname = get(DATA.saveitem,'string');
    Expt = DATA.Expt;
    Expt.Header.CombineDate = now;
    if DATA.datafilename(2) == ':'
        drive = DATA.datafilename(1:2);
    else
        drive = '';
    end
    tic;
    useraw = 0;
    outname = regexprep(outname,'\.c[0-9]\.','.lfp.');
    outname = regexprep(outname,'\.p[0-9]*c[0-9]\.','.lfp.');
    if strncmp(DATA.filetype,'Grid',4)
        Expt = LoadLFP(Expt,'zeropad','ft','verbose');
        name = 'makelfp'; %don't save these outputs
    elseif DATA.state.online
        for j = 1:length(DATA.combineids)
            if max([DATA.probes.source]) >= 1 %if >1 data file, LFP is in #2
                lfpfile = strrep(DATA.Expts{DATA.combineids(j)}.Header.Name,'/Expt','A/LFPS');
                lfpfile = regexprep(lfpfile,'(\.[0-9]*).mat','A$1.lfp.mat');
            else
                lfpfile = strrep(DATA.Expts{DATA.combineids(j)}.Header.Name,'/Expt','/LFPS');
            end
            lfpfiles{j} = [drive lfpfile];
        end
        Expt= LoadSpike2LFP(Expt,'reload','lfpfile',lfpfiles,'drive',drive,'fixshort',10);
        % above reps should fail with online files, but remove .lfp first to be sure
        if ~isfield(Expt.Trials,'LFP')
            fprintf('No LFP Data in %s\n',lfpfiles{:});
            return;
        end
        if DATA.state.scalelfp
            Expt = ScaleLFP(Expt, 'scale', 0.05);
        end
        outname = strrep(outname,'.lfp.mat','.mat');
        outname = strrep(outname,'.mat','.lfp.mat');
    elseif ~isempty(LFP)
        Expt= LoadSpike2LFP(Expt,LFP);
    elseif useraw
        lfpfile = [drive strrep(Expt.Header.Name,'.mat','A.lfp.mat')];
        Expt= LoadSpike2LFP(Expt,'reload','lfpfile',lfpfile,'fixshort',5);
    elseif strncmp(DATA.filetype,'Grid',4)
        Expt = LoadLFP(Expt,'zeropad','ft','verbose');
        name = 'makelfp'; %don't save these outputs
    else
        Expt= LoadSpike2LFP(Expt,'reload','drive',drive,'fixshort',5);
    end
    ltime = toc;
    if isfield(Expt.Trials,'LFP') %successfully got LFP Data
        ck = CheckLFP(Expt,'verbose');
        ns = mode(ck.lens);
        Expt.Trials = rmfields(Expt.Trials,{'Spikes' 'OSpikes' 'Ocodes' 'count'});
        LFP = Expt;
        tic;
        if strncmpi(name,'savelfp',7)
            if sum(ck.lens .* ck.nch) > 100e6
                mycprintf('blue','Removing FTs from LFP before saving\n');
                LFP.Trials = rmfield(LFP.Trials,'FTlfp');
            end
            save(outname,'LFP');
            fprintf('Saved %d samples/trial. %s Loading LFP tookd %.2f, saving %.2f\n',ns,outname,ltime,toc)
        end
        LFP = FixExpt(LFP,'auto')
        setappdata(DATA.toplevel,'LFPExpt',LFP)
        if DATA.plot.lfpplot > 0
            if DATA.state.online
                gains = std(cat(1,LFP.Trials.LFP));
                for j = 1:length(LFP.Trials)
                    for k = 1:length(gains);
                        LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k)./gains(k);
                    end
                end
            end
            cmb.ShowLFPPlot(DATA);
        end
        set(DATA.toplevel,'UserData',DATA);
    else
        fprintf('No LFP Data for %s\n',outname);
    end
    if DATA.logfid
        fprintf(DATA.logfid, '%s, LFP file %s (Took %.1f %.1f)\n',datestr(now),outname,ltime,toc);
    end
    
elseif strncmpi(name,'save',4) & isfield(DATA,'Expt')
    DATA.outname = get(DATA.saveitem,'string');
    Expt = DATA.Expt;
    Expt.Header.CombineDate = now;
    if isfield(DATA,'cluster')
        Expt.Header.Cluster = DATA.cluster;
    end
    Expt.Header.clist = DATA.spikelist;
    id = regexp(DATA.outname,'\.c[0-9]\.');
    if ~isempty(id)
        co = str2num(DATA.outname(id(1)+2:id(1)+3));
        if ~ismember(co,DATA.spikelist)
            a = questdlg(sprintf('Cluster %d Not in SpikeList',co),'Cluster-Name mismatch','Cancel','Proceed','Proceed');
            if strmatch(a,'Cancel');
                return;
            end
        end
    end
    
    if isfield(Expt.Trials,'LFP')
        Expt.Trials = rmfield(Expt.Trials,{'LFP' 'lfptime' 'lfpo' 'FTlfp'});
    end
    
    save(DATA.outname,'Expt');
    nspk = sum([Expt.Trials.count]);
    fprintf('Saved %d spikes (Expts%s) to %s\n',nspk,sprintf(' %d',DATA.combineids),DATA.outname);
    if isfield(Expt.Header,'suffixes')
        estr = sprintf('%s Suffs %s',sprintf(' %d',Expt.Header.Combineids),sprintf(' %d',Expt.Header.suffixes));
    else
        estr = sprintf(' %d',Expt.Header.combineids);
    end
    if DATA.logfid
        fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s) to %s (%s)\n',datestr(now),nspk,estr,DATA.outname,DATA.user);
    end
    cfile = cmb.CombinerLst(DATA);
    if exist(cfile,'file') == 2
        load(cfile);
    else
        combines.name = {};
    end
    id = strmatch(splitpath(DATA.outname),combines.name,'exact');
    if isempty(id)
        id = length(combines.name)+1;
    end
    combines.name{id} = splitpath(DATA.outname);
    combines.lastdate(id) = now;
    combines.list{id} = Expt.Header.Combined;
    combines.starts{id} = Expt.Header.BlockStart;
    combines.trials(id,:) = [Expt.Trials(1).Trial Expt.Trials(end).Trial];
    DATA.combines = combines;
    %tic;
    DATA = cmb.cListExpts(DATA,DATA.Expts,'nowarn');
    % toc
    set(DATA.toplevel,'UserData',DATA);
    %  cmb.combine('listexpts');
    %  save(cfile,'combines');
elseif strncmpi(name,'show',4)
    if iscell(varargin{1})
        for j = 1:length(varargin{1})
            DATA.show.(varargin{1}{j}) = 1;
        end
    else
        DATA.show.(varargin{1}) = 1;
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'setadcid',6)
    DATA.adcid = varargin{1};
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'setvplot',6)
    DATA.plot.voltxy = varargin{1};
    set(DATA.toplevel,'UserData',DATA);
    if DATA.plot.voltxy == 4
        cmb.PlotTrodeXcorr(DATA,0);
    end
elseif strncmpi(name,'pcaplot',6)
    if length(varargin) & isnumeric(varargin{1})
        eid = varargin{1};
    else
        eid = DATA.currentexpt(1);
    end
    cmb.PlotPCAs(DATA,eid);
elseif strncmpi(name,'trodecorr',6)
    if length(varargin) & isnumeric(varargin{1})
        eid = varargin{1};
    else
        eid = DATA.currentexpt(1);
    end
    out= cmb.CalcTrodeCorrs(DATA,eid);
elseif strncmpi(name,'tetrodemov',6)
    ispk = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
    for j = 1:32
        hold off;
        plot(DATA.AllData.Spikes.values(ispk,j),DATA.AllData.Spikes.values(ispk,j+32),'.','markersize',1);
        hold on;
        plot(DATA.AllData.Spikes.values(ispk,j+64),DATA.AllData.Spikes.values(ispk,j+96),'r.','markersize',1);
        
        drawnow;
        pause
    end
elseif strncmpi(name,'setxycluster',6)
    if isnumeric(varargin{1}) & length(varargin{1}) > 1
        if length(varargin{1}) > 2
            cmb.SetXYCluster(DATA,0,5,2,varargin{1});
        else
            cmb.SetXYCluster(DATA,0,4,2,varargin{1});
        end
    elseif ischar(varargin{1})
        a = strmatch(varargin{1},DATA.spkvarnames,'exact')
        b = strmatch(varargin{2},DATA.spkvarnames,'exact')
        if length(a) && length(b)
            cmb.SetXYCluster(DATA,0,4,2,a,b);
        end
    end
elseif strncmpi(name,'sortexpts',5)
    
    
elseif strncmpi(name,'l',7)
    cmb.lProbeClusters(DATA);
elseif strncmpi(name,'store',5)
    set(DATA.toplevel,'UserData',varargin{1});
elseif strncmpi(name,'trackcluster',8)
    id = get(DATA.elst,'value');
    sb = findobj(DATA.xyfig,'tag','Set+Next');
    for j = id:length(DATA.expid)
        DATA.exabsid = DATA.expid(j);
        set(DATA.elst,'value',j);
        cmb.combine('setexpt');
        DATA = get(DATA.toplevel,'UserData');
        C = cmb.OptimizeDprime(DATA);
        DATA.cluster{DATA.currentcluster,DATA.probe} = C;
        [DATA, newd, nc] = SetSpkCodes(DATA,DATA.spklist,DATA.probe,2);
        cmb.DrawClusters(DATA,DATA.cluster, 0);
        set(sb,'backgroundcolor','r');
        drawnow;
    end
    
elseif strncmpi(name,'relist',6)
    if DATA.state.online == 0
        return;
    end
    DATA = cmb.ReadDir(DATA, DATA.name,'setprobe',DATA.probe);
    if  strncmp(DATA.filetype,'Grid',4)
        DATA.grididx = BuildGridIndex(DATA.nevdir, DATA.Expts,'reindex','noerrs');
        DATA.fullvidx = ListExpts(DATA.name,'fullv');
    end
    DATA = cmb.combine('listexps',DATA,'Tag',DATA.tag.top);
    eid = get(DATA.clst,'value');
    DATA = cmb.ListSubExpts(DATA,eid);
    
    set(DATA.toplevel,'UserData',DATA);
    if DATA.state.autoupdatelist
        names = get(DATA.clst,'String');
        eid = strmatch(DATA.Expts{end}.Header.expname, names,'exact');
        if length(eid) == 1
            set(DATA.clst,'value', eid);
            DATA = cmb.ListSubExpts(DATA,eid);
            nex = length(get(DATA.elst,'String'));
            set(DATA.elst,'value',nex);
            set(DATA.toplevel,'UserData',DATA);
        end
        WinFront(DATA.tag);
    end
elseif exist(name,'dir') % do dir before file, since dirs pass exist(name,'file')
    if isfield(DATA,'name') && strcmp(DATA.name,name)
        fprintf('%s already loaded\n',DATA.name);
    else
        if DATA.bysuffix
            DATA.AllClusters = {};
            DATA.state.showspikes = 0;
        end

    DATA.Expts = {};
    args = argon;
    if reindex
        args = {args{:}, 'relist'};
    end
    if DATA.state.online
        DATA.state.online = 1;
        DATA.state.autospool = 0;
        DATA.state.autoupdatelist = 1;
    end
    DATA = cmb.ReadDir(DATA, name, varargin{:});
    if isfield(DATA,'AllClusters')
        for j = 1:length(DATA.AllClusters)
            for k = 1:length(DATA.AllClusters{j})
                if isfield(DATA.AllClusters,'mahal')
                    DATA.Expts{j}.Cluster{k}.mahal = DATA.AllClusters{j}(k).mahal;
                end
            end
        end
    elseif strncmp(DATA.filetype,'Grid',4)
        DATA.grididx = BuildGridIndex(DATA.nevdir, DATA.Expts,'noerrs');
        DATA.fullvidx = ListExpts(DATA.name,'fullv');
        DATA = cmb.LoadClusters(DATA,cmb.ClusterFile(DATA),'allprobes');
    else
        DATA = cmb.LoadClusters(DATA,cmb.ClusterFile(DATA));
    end
    DATA.lastread = 0;
    DATA.lastsize = 0;
    if ~isfield(DATA,'timerobj')  && DATA.state.autolist
        DATA.timerobj = timer('TimerFcn',@cmb.timerfna, 'Period', 1.0,...
            'Tag',DATA.tag.top, 'ExecutionMode','FixedSpacing');
        DATA.lastread = 0;
    end
    DATA.dirname = name;
    set(DATA.toplevel,'UserData',DATA);
    cmb.SetGui(DATA);
    end
elseif strncmpi(name,'newfile',6)
    DATA.datafilename = get(findobj(DATA.toplevel, 'Tag','FileName'),'string');
    cmb.combine(DATA.datafilename);
elseif strncmpi(name,'recombineorbw',12)
    ReClassify
    cmb.ReCombineAll(DATA,0,'orbw');
elseif strncmpi(name,'recombine',6)
    cmb.ReCombineAll(DATA,0);
elseif strncmpi(name,'reptrial',6)
    cmb.PlayOneTrial(DATA,DATA.currenttrial,0);
elseif strncmpi(name,'rfupdate',6)
    if isfield(DATA.fitvals,'Pp') & isfield(DATA.fitvals,'Op')
        if abs(DATA.fitvals.OpRo-DATA.fitvals.PpRo) > 2
            msgbox('Ro mismatch');
        end
        xy = op2xy([DATA.fitvals.Op DATA.fitvals.Pp],DATA.fitvals.OpRo);
        uflname = strrep(DATA.datafilename,'.mat','.ufl');
        fid = fopen(uflname,'a');
        if fid > 0
            fprintf(fid,'FitRF %.2f %.2f %.2f %.2f\n',xy(1),xy(2),DATA.fitvals.Opw,DATA.fitvals.Ppw);
            fclose(fid);
        end
    end
elseif strncmpi(name,'TrackTemplates',7)
    if ~isfield(DATA,'TemplateScores') | strcmpi(name,'TrackTemplateForce')
        cmb.TrackTemplates(DATA);
        DATA = get(DATA.toplevel,'UserData');
    end
    if length(varargin) & isnumeric(varargin{1})
        cl = varargin{1};
        first = 2;
    else
        cl = 2;
        first = 1;
    end
    cmb.PlotTemplateScores(DATA,cl,varargin{first:end});
elseif strncmpi(name,'winfront',6)
    WinFront(DATA.tag);
elseif strncmpi(name,'xcorr',4)
    type = strmatch(name,{'xcorrdc' 'xcorrv' 'xcorr'},'exact');
    probes = varargin{1};
    
elseif exist(name,'file')
    args = argon;
    if reindex
        args = {args{:}, 'relist'};
    end
    if setexpts
        DATA.exabsid =setexpts;
    else
        DATA.exabsid = [];
    end
    
    DATA.appending = 0;
    DATA = cmb.ReadFile(DATA, name, args{:});
    X = [];
    X = CopyFields(X,DATA,{'datadir' 'datafilename'});
    for j = 1:length(appendfiles)
        DATA.appending = 1;
        DATA.appendexptid(j) = length(DATA.Expts)+1;
        DATA = cmb.ReadFile(DATA, appendfiles{j}, args{:});
        DATA.appenddir{j} = DATA.datadir;
    end
    DATA = CopyFields(DATA,X,{'datadir' 'datafilename'});
    set(DATA.toplevel,'UserData',DATA);
    if plotallspikes
        cmb.PlotAllExptsAndProbes(DATA);
    end
    if recombine == 3
        cmb.ReCombineAll(DATA,0,'orbw');
    elseif recombine == 4
        cmb.ReCombineByName(DATA,recombinename,1);
    elseif recombine > 0
        cmb.ReCombineAll(DATA,0);
    end
    if plotall
        d = fileparts(DATA.datafilename);
        PlotAllProbes(d,'save','sptrig')
    end
    DATA = cmb.PrintComments(DATA);
    if strfind(name,'ic-169')
        cmb.CheckSaccades(DATA.Expts,'Z:/smr/icarus/169/ic169.em.mat');
    end
    if ismember(recombine,[2 3 4])
        cmb.combine('close');
    end
elseif iscell(name)
    args = argon;
    if reindex
        args = {args{:}, 'relist'};
    end
    names = name;
    for j = 1:length(names)
        name = names{j};
        if exist(name,'file')
            DATA = cmb.ReadFile(DATA, name, args{:});
            set(DATA.toplevel,'UserData',DATA);
            if recombine == 1
                cmb.ReCombineAll(DATA,0);
                if plotall
                    d = fileparts(DATA.datafilename);
                    PlotAllProbes(d,'save','sptrig')
                end
            end
        end
    end
    
elseif regexp(name,'ruf[0-9][0-9][0-9]') %if its not a file name, try building path
    dnum = strrep(name,'ruf','');
    set(DATA.toplevel,'UserData',DATA);
    name = ['/data/rufus/' dnum '/' name '.mat'];
    if exist(name,'file')
        cmb.combine(name,'Tag',DATA.tag.top);
    end
elseif firstcall
    fprintf('No file or directory %s\n',name);
    close(DATA.toplevel);
    return;
end


%
%
%
%END OF PARSING name
%
%
%






%some varaings have to be processed after the file has been opened
recombinenames = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'readconfig',8)
        if length(varargin) > j && ischar(varargin{j+1}) & exist(varargin{j+1},'file')
            j = j+1;
            DATA  = ReadConfig(DATA,'file', varargin{j});
        else
            DATA  = ReadConfig(DATA);
        end
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'loadclusters',8)
        ts= now;
        for j = 1:length(DATA.Expts)
            fprintf('Loading Expt %d',j);
            [DATA, D] = cmb.CheckClusterLoaded(DATA,j,DATA.probe);
            loaddur(j) = D.loaddur(end);
            fprintf(' took%s\n',sprintf(' %.2f',D.loaddur));
        end
        C = DATA.AllClusters;
        x = whos('C');
        byterate = x.bytes/(1024 * 1024 * sum(loaddur));
        fprintf('Took %.2f sec (%.2f in LoadClusters %.2f Mb/sec)\n',mytoc(ts),sum(loaddur),byterate);
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(varargin{j},'combinemucells',12)
        if length(varargin) > j
            args = varargin(j+1:end);
        else
            args = {};
        end
        varargout{1} = cmb.CombineAllCells(DATA,'mucells', args{:});
    elseif strncmpi(varargin{j},'nospikes',6)
        DATA.state.showspikes = 0;
    elseif strncmpi(varargin{j},'recombinemucells',12)
        if length(varargin) > j
            args = varargin(j+1:end);
        else
            args = {};
        end
        if strncmpi(varargin{j},'recombinemucellsquick',16)
            args = {'quick' 'saveonlyall' args{:}};
        end
        varargout{1} = cmb.ReCombineAll(DATA,0,'mucells',recombinenames, 'reapply', args{:});
        varargout{2}.toplevel = double(DATA.toplevel);
    elseif strncmpi(varargin{j},'recombinename',9)
        j = j+1;
        if iscell(varargin{j})
            recombinenames = varargin{j};
        else
            recombinenames = {recombinenames{:} varargin{j}};
        end
    elseif strncmpi(varargin{j},'recombine',6)
        varargout{1} = cmb.ReCombineAll(DATA,0);
    elseif strncmpi(varargin{j},'quit',4)
        cmb.CloseCombine(DATA);
        return;
    elseif strncmpi(varargin{j},'plotargs',8)
        j = j+1;
        DATA.plotargs = {};
        while j <= length(varargin)
            DATA.plotargs = {DATA.plotargs{:} varargin{j}};
            j = j+1;
        end
        set(DATA.toplevel,'UserData',DATA);
    end
    j= j+1;
end
cmb.SetGui(DATA);
cmb.NotBusy(DATA);
if DATA.state.verbose
    fprintf('Leaving Combine (%s)\n',calltype);
end

if nargout && isempty(varargout)
    varargout{1} = DATA;
end


