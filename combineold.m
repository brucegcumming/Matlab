function out = combine(name, varargin)
%combine(file)
%cuts clusters and combines expts for Spike2 generated matlab data files
%
%

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
nextcell = 0;
argon = {};
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
        elseif strncmpi('expts',varargin{j},4)
            j = j+1;
            setexpts = varargin{j};
            if exist('DATA','var') & isfield(DATA,'elst')
                set(DATA.elst,'value',setexpts);
                DATA.currentexpt = setexpts(1);
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
DATA.expstrs = {'dx' 'xo' 'dp' 'dO'  'dxP' 'dxXceC' 'dxXce' 'sf' 'dxXId' 'dxXIdP' 'or' 'Op' 'Pp' 'sz' 'jv'...
    'orXob' 'orXobP' 'orXme' 'szXob' 'DcXorP' 'DcXor' 'DcXorXUsPRC' 'DcXorXUsP' 'OpXme' 'PpXme' 'co' 'CtfXip' 'tfXip' 'IB' 'backMov' 'backMovXannTyp' 'IBXannTyp' 'dirXme' ...
    'OpXip' 'PpXip' 'me' 'OlXOr' 'ce' 'dxXUs' 'dxXUsP' 'dxXpi' 'dxXor' 'orXdx' 'TwoCylDispXhxPRB' 'TwoCylDispXhxP' 'TwoCylDispXhx' 'dxXIdPD' 'dxXIdD'  ...
    'orXobXUs' 'orXobXUsP' 'dxXceXFr' 'yo' 'orXobXfyP' 'orXobXfy'...
    'orXobXfxP' 'orXobXfx' 'orXId' 'dirXId' 'orXIdD' 'dirXIdD' ...
    'szXobXme' 'orXip' 'dxXceXrC' 'meXFr' 'dxXces' 'dxXcesRC' 'dxXmixac' 'dxXcesXFr' ...
    'dxPRB' 'dxRB' 'bdXIdP' 'dxXIdPDB' 'sMXId' 'sMXdx' ...
    'sMXPd'};
DATA.expnames = {'DT', 'XO', 'DP' 'ODX' 'DTP' 'AC' 'AC' 'SF' 'DTID' 'ABD' 'OT' 'OP' 'PP' 'SZ' 'VE' ...
    'ORBW' 'ORBW' 'OXM' 'SZOB' 'DCOR' 'DCOR' 'UDCOR' 'UDCOR' 'OPM' 'PPM' 'CO' 'CTFIP' 'TFIP' 'IB' 'Movie' 'MovieANN' 'IB' ....
    'OD' 'OPIP' 'PPIP' 'OD' 'BOP' 'CE' 'UADT' 'UADT' 'DPI' 'DXO' 'DXO' 'TWO' 'TWO' 'TWO' 'DID' 'DID' 'UORBW' 'UORBW' 'FAC' 'YO' ...
    'ORBW' 'ORBW' 'ORBW' 'ORBW' 'SRID' 'SRID' 'DRID' 'DRID' ...
    'SFMOB' 'ORP' 'AC' 'FOD' 'XAC' 'XACRC' 'MIXAC' 'FXAC' 'DXRB' 'DXRB' 'BDID' 'DIDB' 'FLID'...
    'DTPL' 'PDPL'};
DATA.expstrs = {DATA.expstrs{:} 'sfXob' 'sfXme' 'szXme' 'dxXdx' 'dOXdO' 'dOXceXFr' 'dOXce'};
DATA.expnames = {DATA.expnames{:} 'SFOB' 'SFM' 'SZM' 'DT' 'ODX' 'FOXAC' 'OXAC'};
    DATA.subprobe = 0;
    DATA.subprobes = 0;
    DATa.xprobes = [];



[DATA.spkvarnames, DATA.spkvarorder] = GetSpikeVals(DATA,NaN,NaN,NaN,[]);
    DATA.suffs = 'abcdeghijklmnopqrstuvwxyz';
    DATA.state.recut = 1; %default
    DATA.state.resetclusters = 0;
    DATA.state.recount = 0;
    DATA.state.plotpsych = 0;
    DATA.state.psychonly = 0;
    DATA.state.plotcombined = 0;
    DATA.state.autoplot = 0;
    DATA.state.autolist = 0;
    DATA.state.autoplotcells = 0;
    DATA.state.autoplotnewprobe = 1;
    DATA.state.autoreplotgraph = 1;
    DATA.state.autosetlist = 1;
    DATA.state.NeedAllProbeFile = 0;
    DATA.state.plotseq  = 0;
    DATA.state.uselfp = 0;
    DATA.state.listlen = 6;
    DATA.state.bigwindow(1) =1.1;
    DATA.state.bigwindow(2) =1.1;
    DATA.state.includeprobename = 0;
    DATA.state.applylastcluster = 0;
    DATA.state.classifyallprobes = 0;
    DATA.state.verbose = 0;
    DATA.tag.top = TOPTAG;
    DATA.tag.options = [TOPTAG 'Options'];
    DATA.tag.spikev = [TOPTAG 'SpikeV'];
    DATA.tag.clusterxy = [TOPTAG 'ClusterPlot'];
    DATA.tag.emplot = [TOPTAG 'EM'];
    DATA.tag.lfpplot = [TOPTAG 'LFP'];
    DATA.tag.rcfiga = [TOPTAG 'RCa'];
    DATA.tag.rcfigb= [TOPTAG 'RCb'];
    DATA.tag.showvals= [TOPTAG 'ShowVal'];
    DATA.tag.celllist= [TOPTAG 'CellList'];
    DATA.tag.psych = [TOPTAG 'Psych'];
    DATA.tag.cptag = [TOPTAG 'CombineCP'];
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.tag.dataplot = [TOPTAG 'Combiner'];

    DATA.wsc= [1 1];
    DATA.defaults.fz = 96;
    DATA.state.fixrange = 0;
    DATA.spooling = 0;
    DATA.densityplot = 0;
    DATA.alldensityplot = 0;
    DATA.plot.showspikeshapes = 0;
    DATA.plot.clusterXrange = [0 10];
    DATA.plot.clusterYrange = [0 1.5];
    DATA.plot.clusterZrange = [0 10];
    DATA.cluster{1,1}.Arange = [6:8];
    DATA.cluster{1,1}.Brange = [11:20];
    DATA.cluster{1,1}.Erange = [1:100];
    DATA.clusterArange = [6:8];
    DATA.clusterBrange = [11:20];
    DATA.clusterErange = [1:100];
    DATA.cluster{1,1}.params = [1 2];
    DATA.cluster{1,1}.autocut = 0;
    DATA.plot.autoscale = 1;
    DATA.plot.autoscalemode = 1;
    DATA.state.autofit = 0;
    DATA.state.autonext = 0;
    DATA.state.autospool = 0;
    DATA.state.showspkxy = 0;
    DATA.state.forcebuild = 0;
    DATA.state.preperiod = 2000;
    DATA.state.postperiod = 2000;
    DATA.state.reclassify = 0;
    DATA.state.redoautocut = 0;
    DATA.plot.nodc = 1;
    DATA.plot.autoclustermode = 1;
    DATA.test.fastplot = 0;
    DATA.state.showonlineclusters = 0;
    DATA.state.savedvdt = 0;
    DATA.show.sz = 0;
    DATA.show.or = 0;
    DATA.show.jv = 0;
    DATA.show.jx = 0;
    DATA.show.sf = 0;
    DATA.show.me = 0;
    DATA.show.tf = 0;
    DATA.show.dw = 0;
    DATA.show.Dw=0;
    DATA.show.Dm=0;
    DATA.show.ed = 1;
    DATA.show.Ri = 0;
    DATA.show.ob = 0;
    DATA.show.Ro = 0;
    DATA.show.Fr = 0;
    DATA.show.ei = 0;
    DATA.show.dx = 0;
    DATA.show.js = 0;
    DATA.show.times = 0;
    DATA.show.backxo = 0;
    DATA.show.backyo = 0;
    DATA.show.co = 0;
    DATA.show.Bc = 0;
    DATA.show.bo = 0;
    DATA.xyfig = 0;
    DATA.optionfig = 0;
    DATA.svfig = 0;
    DATA.plot.showem = 0;
    DATA.plot.showcp = 0;
    DATA.plot.showsync = 0;
    DATA.spikelist = 1;
    DATA.plot.dvdt = 0;
    DATA.plot.voltxy = 0;
    DATA.adcid = [8 40 72 100];
    DATA.minplottime = 0.1;
    DATA.plot.clusterX = 1;
    DATA.plot.clusterY = 2;
    DATA.plot.clusterZ = 0;
    DATA.plot.nmin = 0;
    DATA.plot.nminrc = 0;
    DATA.plot.sdfw = 166;
    DATA.plot.showISI = 0;
    DATA.plot.acov = 0;
    DATA.plot.collapse = 0;
    DATA.plot.flip = 0;
    DATA.plot.bsdelaymax = 5000;
    DATA.plot.setptsize = 0;
    DATA.plot.condenseRC = 0;
    DATA.plot.SpikeMaxV = 5;
    DATA.plot.SpikeVsep = 3;
    DATA.plot.lfpplot = 0;
    DATA.plot.plotmod = 0;
    DATA.plot.xcorr = 0;
    DATA.plot.addhash = 0;
    DATA.probenames = {'1' '2'};
    DATA.probelist = [5];
    DATA.probe = 5;
    DATA.linesread = 0;
    DATA.currentcluster = 1;
    DATA.forceclusterid = 0;
    DATA.state.online = 0;
    DATA.outname = 'tmp.mat';
    DATA.badnames = {};
    DATA.allcombineids = {};
    DATA.fitvals = [];
    DATA.logfid = 0;
    DATA.plot.lfptype = 0;
    DATA.state.scalelfp = 0;
    DATA.state.autoplotallspikes = 0;
    DATA.state.useonlineclusters = 1;
    DATA.filetype = 'Spike2';
    if ischar(name)
        DATA.datafilename = name;
        if exist(name,'dir') % do dir before file, since dirs pass exist(name,'file')
            DATA.state.online = 1;
            DATA.state.autospool = 1;
            DATA.state.autosetlist = 1;
        end
    else
        DATA.datafilename = 'List';
    end
    DATA.meanspkfile = strrep(DATA.datafilename,'.mat','spks.mat');
    
    
    DATA.syncsign = 2;
    DATA.TriggerSign = 0;
    DATA.plot.syncoverlay = 1;
    DATA.plot.synccluster = 0;
    DATA.plot.showwave = 0;
    DATA.plot.timebyspikeprobe = 0;
    DATA.plot.centerRFmeasures = 0;
    load('StdTemplate');
    DATA.Templates = Templates;
    DATA.TemplateInfo = [];
    DATA.tag.allprobes = 'AllProbeSpikes';
    DATA.currenttrial = 0;
    DATA.firsttrial = 0;
    DATA.lasttrial = 0;
    DATA.playingspk = 0;
    DATA.state.nospikes = 0;
    DATA.state.optimizeclusters = 0;
    DATA.plot.UsePreviousCut = 1;
    DATA.plot.showartifacts = 0; %display events classified as artifacts
    DATA.plot.showN = 1;
    DATA.plot.showdprime = 0;
    DATA.plot.type = 'Mean';
    DATA.plot.autoVrange = 0;
    DATA.state.showspikes = 0;
    
    DATA.state.currentclusterquality = 0;
    DATA.figs.allprobes = 0;
    DATA.plot.emskip = 0;
    DATA.em.stattime = [1500 21000]; %time window to look for saccades
    DATA.em.maxdrift = 0;
    DATA.AllData.pcs = [];
    DATA.savedclusters = 1; %no changes made yet
    DATA.plot.prettyfigs = 0;
    DATA.plot.voffsets = [];
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
j = 1;
while j <= nargin -1
    if isstruct(varargin{j})
        if isfield(varargin{j},'Trials')
            LFP = varargin{j};
        end
    elseif strncmpi(varargin{j},'expand',5)
        j = j+1;
        DATA.state.bigwindow(1) = varargin{j}(1);
        if length(varargin{j}) > 1
            DATA.state.bigwindow(2) = varargin{j}(1);
        end
    elseif strncmpi(varargin{j},'preprocess',6)
        preprocess = 1;
    elseif strncmpi(varargin{j},'plotspikes',6)
        plotallspikes = 1;
    elseif strncmpi(varargin{j},'plotall',6)
        plotall = 1;
    elseif strncmpi(varargin{j},'pretty',4)
        DATA.plot.prettyfigs = 1;
        set(DATA.toplevel,'UserData',DATA);
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
        SpoolSpikes(DATA);
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
    DATA = ReadFile(DATA, name);
    DATA = PrintComments(DATA);
    DATA = CalcOnlineShapes(DATA);
    return;
end
if firstcall
        DATA = BuildGUI(DATA);
        if length(setexpts)
            set(DATA.elst,'value',setexpts);
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
    DATA = ListExpts(DATA,Expts);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'autocut',6)
    eid = get(DATA.elst,'value');
    for j = 1:length(eid)
        AutoCut(DATA,DATA.expid(eid(j)),eid(j));
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
        DATA = GetAllProbeFig(DATA);
        for j = 1:length(eid)
            for k = pid
                DATA.probe = DATA.probelist(k);
                subplot(nr,nc,k);
                DATA = AutoCut(DATA,DATA.expid(eid(j)),eid(j),'cluster',args{:});
            end
        end
        DATA.probe = oldprobe;
    else
    for j = 1:length(eid)
        DATA = AutoCut(DATA,DATA.expid(eid(j)),eid(j),'cluster',args{:});
    end
    end
    DATA.AutoClusters = SaveAutoClusters(DATA);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'AddProbeList',5)
    AddProbeList(DATA);
    DATA.plot.showsync = 1;
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
    PlotDDF(DATA);
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
       PrintLogData(DATA,1);
       PrintLogData(DATA,2);
       PrintLogData(DATA,5);
       
elseif strncmpi(name,'listexps',6)
    eid = get(DATA.clst,'value');
    if eid(1) > 1
        suff = '';
        it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
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
            if strmatch(expname,'CO','exact') & strmatch(stimname,'square')
                expname = 'FLSH';
            end
            if regexp(DATA.datafilename,'Expt[0-9]*.mat') %online file
                DATA.outname = [DATA.datafilename '/c' num2str(DATA.spikelist(1)) '.' stimname '.' DATA.expnames{it} suff '.mat'];
            else
                DATA.outname = CombinedName(DATA,eid(1),DATA.spikelist(1));
            end
            if exist(DATA.outname,'file')
                GetFigure(DATA.tag.dataplot);
                load(DATA.outname);
                Expt.Header.Name = BuildName(Expt.Header.Name);
                DATA.Expt = Expt;
                args = PlotArgs(DATA, Expt,'combined');
                DATA.plotres = PlotExpt(Expt,args{:});
                csuffs = [];
                if isfield(Expt.Header,'Combineids')
                    for j = 1:length(Expt.Header.Combineids)
                        ci = Expt.Header.Combineids(j);
                        csuffs = [csuffs num2str(ci) ' '];
                        if isfield(Expt.Header,'Clusters') & length(Expt.Header.Clusters) >= j & ...
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
            DATA.outname = [strrep(DATA.datafilename,'.mat','.') 'c' num2str(DATA.spikelist(1)) '.' stimname '.TMP' suff '.mat'];
        end
    else
            DATA.outname = 'tmp.mat';
    end
    DATA = ListSubExpts(DATA,eid);
    set(DATA.saveitem,'string',DATA.outname);
    set(DATA.toplevel,'UserData',DATA);
    if nargout
        out = DATA;
    end
elseif strncmpi(name,'getexpt',6)
    out = DATA.Expt;
    NotBusy(DATA);
    return;
elseif strncmpi(name,'getstate',6)
    out = DATA;
    NotBusy(DATA);
    return;
elseif strncmpi(name,'ingnoreonline',9)
    DATA.state.useonlineclusters = 0;
elseif strncmpi(name,'checklists',6)
    CheckLists(DATA);
elseif strncmpi(name,'ClearClusters',7)
    eid = DATA.currentexpt;
    nc = CountClusters(DATA.Expts{eid}.Cluster);
    for k = 1:nc;
        for j = 1:size(DATA.Expts{eid}.Cluster,2);
            if ~isempty(DATA.Expts{eid}.Cluster{k,j}) & isfield(DATA.Expts{eid}.Cluster{k,j},'x')
            DATA.Expts{eid}.Cluster{k,j}.touched = 0;
            DATA.Expts{eid}.Cluster{k,j} = rmfield(DATA.Expts{eid}.Cluster{k,j},{'x' 'y'});
            end
        end
    end
    nc = CountClusters(DATA.cluster);
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
    nc = CountClusters(DATA.Expts{eid}.Cluster);
    for k = 1:nc;
        if ~isempty(DATA.Expts{eid}.Cluster{k,p}) & isfield(DATA.Expts{eid}.Cluster{k,p},'x')
            DATA.Expts{eid}.Cluster{k,p}.touched = 0;
            DATA.Expts{eid}.Cluster{k,p} = rmfield(DATA.Expts{eid}.Cluster{k,p},{'x' 'y'});
        end
    end
    end
    nc = CountClusters(DATA.cluster);
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
    PrintComments(DATA);
elseif strncmpi(name,'combine',6)
    DATA.extype = get(DATA.clst,'value');
    [Expt, DATA] = CombinePlot(DATA,1);
    DATA.Expt = Expt;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'options',5)
    DATA.optionfig = setoptions(DATA,DATA.tag.options);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'loadclusters',6)
    DATA = LoadClusters(DATA,ClusterFile(DATA));
elseif strncmpi(name,'meanspike',6)
    ts = now;
    if length(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
        DATA = CalcMeanSpike(DATA,varargin{1});
    else
        DATA = CalcMeanSpike(DATA,DATA.expid);
    end
    mytoc(ts);
    set(DATA.toplevel,'UserData',DATA);
    if strfind(name,'save')
        SaveSpikeShape(DATA,DATA.meanspkfile);
    end
elseif strncmpi(name,'plotshape',6)
    GetFigure('ClusterShape');
    PlotSpikeShapes(DATA.MeanSpike,varargin{:});
elseif strncmpi(name,'PlotISI',5)
    GetFigure('ISI');
    [isis, t, s] = CalcISI(DATA.Expts{DATA.currentexpt}.Trials);
    id = find(isis < 1000)
    hist(isis(id),100);
    [a,b] = sort(isis);
%    DATA.isit = t(b);
    DATA.isis = s(b);
    DATA.ISIpair = 1;
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'reclassify',5)
    ts = now;
    DATA = ReClassifyAll(DATA);
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

    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'showvals',5)
    DATA.showid = setshow(DATA,DATA.tag.showvals);
    set(DATA.toplevel,'UserData',DATA);
elseif strncmpi(name,'spikes',5)
    PlotSpikes(DATA,varargin{1});
elseif strncmpi(name,'setexp',6)
    if DATA.playingspk & ishandle(DATA.svfig)
        set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
        set(findobj(DATA.svfig,'Tag','StopSpool'),'UserData',1);
        DATA.playingspk = 0;
        set(DATA.toplevel,'UserData',DATA);
        return;
    end
    ts  = now;
    id = get(DATA.elst,'value');
    DATA.exabsid = DATA.expid(id);
    GetFigure(DATA.tag.dataplot,'noforce');
    if DATA.state.plotcombined
        rc = DATA.plot.condenseRC;
        DATA.plot.condenseRC = 0;
% why was this plotting the combined data? Want to plot individual ones
% but using the combined plott style (legends, colors, so that can see what
% is goin on in individual exps
%        PlotCombined(DATA, DATA.Expt);
       for j   = 1:length(DATA.exabsid)
           if j > 1
               PlotCombined(DATA, DATA.Expts{DATA.exabsid(j)},'hold');
           else
               PlotCombined(DATA, DATA.Expts{DATA.exabsid(j)});
           end
       end
        DATA.plot.condenseRC = rc;
        NotBusy(DATA);
        if nargout 
            out = DATA;
        end
        return;
    end

    hs = 'nohold';
    colors = mycolors;
    if strncmpi(name,'setexpplot',8)
        playspk = 0;
    else
        playspk = get(findobj(DATA.toplevel,'Tag','ShowSpikes'),'value');
        DATA.state.showspikes = playspk;
    end
    if length(varargin) && strncmpi(varargin{1},'flagchange',6)
        playspk = 0; %don't spool spikes every time a flag is changed
    end
    if DATA.state.online 
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
    lastexpt = DATA.currentexpt;
    DATA.currentexpt = DATA.expid(id(1));
    DATA.currenttrial = 0; 
    DATA.firsttrial = 0; %no partial list
    DATA.lasttrial = 0;
    CheckSpoolButton(DATA);
    if ~isempty(findobj('Tag',DATA.tag.celllist))
        GetFigure(DATA.tag.celllist);
        PlotCellList(DATA);
    end
    if strncmp(DATA.filetype,'Grid',4) & playspk
        DATA = ReadGridFile(DATA);
        DATA.plot.useprobe = zeros(size(DATA.probelist));
%        DATA.plot.useprobe(DATA.probe) = 1;
    elseif DATA.state.online & lastexpt ~= DATA.currentexpt
        DATA = LoadSpikes(DATA, DATA.currentexpt);
    end

% if not classified, and autospool is set, show spikes. Unless spkXY is set, in whcih case just
% show XY plot
    if DATA.state.online < 2 & DATA.Expts{DATA.expid(id(1))}.gui.classified == 0 && ...
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
    if DATA.state.nospikes %override all the possibilites
        playspk = 0;
    end
    if isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
        DATA.allexp = DATA.currentexpt;
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


    ndef = 0;
    probe = DATA.probe;
    for j = 1:size(DATA.cluster,2)
        pdef(j) = 0;
        ndef(j) = 0;
        for k = 1:size(DATA.cluster,1)
        if isfield(DATA.Expts{DATA.currentexpt},'Cluster') && iscluster(DATA.Expts{DATA.currentexpt}.Cluster,1,j)
            DATA.cluster{k,j} = DATA.Expts{DATA.currentexpt}.Cluster{1,j};
            DATA.cluster{k,j}.touched = 1;
            ndef(j) = ndef(j)+1;
            pdef(j) = pdef(j)+1;
        else
            DATA.cluster{k,j}.touched = 0;
        end
        end
        if ndef(j) == 0 && iscluster(DATA.cluster,1,j) & DATA.state.applylastcluster
            for k = 1:size(DATA.cluster,1)
                if iscluster(DATA.cluster,k,j)
                DATA.Expts{DATA.currentexpt}.Cluster{k,j} = DATA.cluster{k,j};
                ndef(j) = ndef(j)+1;
                end
            end
        end
    end
    if DATA.state.classifyallprobes & isfield(DATA,'AllSpikes')
        for j = 1:size(DATA.cluster,2)
            if ndef(j) & ~isempty(DATA.AllSpikes{j})
                DATA.probe = j;
                DATA = SetExptSpikes(DATA,DATA.currentexpt,0,'useexp');
%                DATA.Expts{DATA.currentexpt}.Cluster{1,j}.touched = 1;
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
    
    DATA.nclusters = CountClusters(DATA.cluster);



    ei = DATA.expid(id(1));
    if playspk & DATA.probe == 100
        DATA = PlaySpikes(DATA, DATA.expid(id(1)));
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
                [DATA.plot.clusterX , DATA.plot.clusterY, DATA] = GetClusterSpace(DATA, DATA.Expts{DATA.expid(id(1))});        
            end
            DATA = PlaySpikes(DATA, DATA.expid(j));
        end
        if DATA.state.showspkxy & isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
           DATA =  PlotAllProbeXY(DATA,0);
        end
        set(DATA.toplevel,'UserData',DATA);
        if ~DATA.state.autoreplotgraph %only show spikes
            NotBusy(DATA);
            out = DATA;
            return; % This stops plotting of results. Why? ? only return if not auto
        end
    elseif DATA.state.showspkxy
%TODO set nclusters, ed properly. And set recut to 2 if necessary.
%            or maybe call SetExptSpks here....

        nclusters = 2;
        if isfield(DATA,'AllClusters') | isfield(DATA,'AllSpikes')
            DATA =  PlotAllProbeXY(DATA,0);
            [DATA, ispk] = SetExptSpikes(DATA,DATA.expid(id(1)),0,'useexp');
            DATA.spklist = ispk;
            DATA.Expts{DATA.currentexpt}.gui.spks = ispk;
            GetFigure(DATA.tag.clusterxy);
            hold off;
            DATA = DrawXYPlot(DATA, ispk);
            ispk = [];
        else
            if ~isfigure(DATA.svfig)
                DATA = PlaySpikes(DATA,ei);
            end
            SetSetBox(DATA,ei);
            GetFigure(DATA.tag.clusterxy);
            %if use 'useexpall' then no clusters carry over, even if it has
            %neve been touched
            [DATA, ispk] = SetExptSpikes(DATA,DATA.expid(id(1)),0,'useexp');
            DATA.spklist = ispk;
            hold off;
            DATA = DrawXYPlot(DATA, ispk);
            nclusters = max(DATA.AllData.Spikes.codes(ispk,2));
        end
        stimes(1) = DATA.Expts{ei}.Trials(end).Start(1);
        stimes(2) = DATA.Expts{ei}.Trials(end).End(end);
        tspk = FindSpikes(DATA, stimes, DATA.probe,ispk);
        DrawSpikeWaves(DATA, tspk, nclusters, 2);
    elseif isfield(DATA.Expts{ei},'Cluster')
        if isfigure(DATA.xyfig)
        set(0,'CurrentFigure',DATA.xyfig);
        hold off; 
        DATA.cluster = DATA.Expts{ei}.Cluster;
        DrawClusters(DATA,DATA.Expts{ei}.Cluster, 0);
        end
    end
%  clearing spklist causes the problem with online files where the spik
%  display changes afterthe cluster is cut. If Spklist is empty, then when
%  the cluster is cut, all spikes (including those not in trials) are
%  shown, until the spool button is set. What goes wrong if this is not
%  cleared?
%    DATA.spklist = []; %once clusters are cut, this should be in them. Don't let first expt set for all
    
    if isExptCluster(DATA.Expts{ei},DATA.currentcluster,DATA.probe) & ...
        isfield(DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe},'quality');
        a = DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe}.quality;
        if(a > 0) %if its not set in cluster, let it inherit value from last expt
        set(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value',a+1);
        end
    else
        DATA.Expts{ei}.Cluster{DATA.currentcluster,DATA.probe}.quality = 0;
    end
    GetFigure(DATA.tag.dataplot,'noforce');
    ClearPlot;
    h = [];
    xargs = {};
    if length(id) > 1 && DATA.Expts{DATA.expid(id(1))}.Header.rc
        xargs = {xargs{:} 'condense'};
    end
    eid = get(DATA.clst,'value');
 
    for j = id(1:end)
        if DATA.state.recount | DATA.Expts{DATA.expid(j)}.gui.counted == 0 & eid(1) > 1
            DATA = CountSpikes(DATA,DATA.expid(j));
        end
        
        s = ExptComments(DATA.Expts{DATA.expid(j)});
 
        te = now;
        edepth = GetEval(DATA.Expts{DATA.expid(j)},'ed');
       if length(id)  == 1 || eid(1) > 1 %don't add plots from 'All' list
            args = PlotArgs(DATA,DATA.Expts{DATA.expid(j)});
            if j > id(1)
                args = {args{:} 'hold'};
            end
        GetFigure(DATA.tag.dataplot);
        res = PlotExpt(DATA.Expts{DATA.expid(j)},hs,'forcecolor',colors{j},args{:},'legendpos',7,xargs{:});
         fprintf('%.2f ',(now-te)*60*60*24);
        sp  = strfind(DATA.subexplist{j},' ');
        if sp
            names{j} = DATA.subexplist{j}(1:sp(1)-1);
        else
            names{j} = DATA.subexplist{j};
        end
        names{j} = [names{j} ShowString(DATA,DATA.Expts{DATA.expid(j)})];
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
    DATA.spikelist = WhichClusters(DATA.toplevel);
    title([t sprintf(' P%d',DATA.probe) 'Cl' sprintf(' %d',WhichClusters(DATA.toplevel)) sprintf('ed%.2f',edepth)]);
%    if DATA.state.verbose
    if 1
    fprintf('Took %.2f\n',(now-ts)*60*60*24);
    end
    out = DATA;
    set(DATA.toplevel,'UserData',DATA);
    
elseif strncmpi(name,'Close',4)
    if isfield(DATA,'timerobj')
        stop(DATA.timerobj);
    end
    if isfield(DATA,'tag')
        names = fieldnames(DATA.tag);
        for j = 1:length(names)
            CloseTag(DATA.tag.(names{j}));
        end
    end
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
elseif strncmpi(name,'showdprime',9)
    DATA.plot.showdprime = 1;
    set(DATA.toplevel,'UserData',DATA);
    fprintf('Showing additional dprime plots\n');
    dp = CalcDprime(DATA,1);
    fprintf('Dp = %.3f\n',dp);
    dp = CalcDprime(DATA,1,'yparam',8); %spk var
    fprintf('Dp = %.3f\n',dp);
    dp = CalcDprime(DATA,1,'yparam',44); %spk var/sqrt(energy)
    fprintf('Dp = %.3f\n',dp);
    dp = CalcDprime(DATA,1,'yparam',45); %spk var/sqrt(energy)
    fprintf('Dp = %.3f\n',dp);
    dp = CalcDprime(DATA,1,'xparam',46,'yparam',44); %spk var/sqrt(energy)
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
    if DATA.state.online
        for j = 1:length(DATA.combineids)
            if max([DATA.probes.source]) >= 1 %if >1 data file, LFP is in #2
            lfpfile = strrep(DATA.Expts{DATA.combineids(j)}.Header.Name,'/Expt','A/LFPS');
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
    else   
        Expt= LoadSpike2LFP(Expt,'reload','drive',drive,'fixshort',5);
    end
    ltime = toc;
    if isfield(Expt.Trials,'LFP') %successfully got LFP Data
    ck = CheckLFP(Expt,'verbose');
    ns = mode(ck.lens);
    Expt.Trials = rmfield(Expt.Trials,{'Spikes' 'OSpikes' 'Ocodes' 'count'});
    LFP = Expt;
    if sum(ck.lens .* ck.nch) > 100e6
        LFP.Trials = rmfield(LFP.Trials,'FTlfp');
    end
    tic;
    if strncmpi(name,'savelfp',7)
    save(outname,'LFP');
    end
    if DATA.plot.lfpplot > 0
        DATA.LFP = LFP;
        if DATA.state.online
            gains = std(cat(1,DATA.LFP.Trials.LFP));
            for j = 1:length(DATA.LFP.Trials)
                for k = 1:length(gains);
                    DATA.LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k)./gains(k);
                end
            end
        end
        set(DATA.toplevel,'UserData',DATA);
        ShowLFPPlot(DATA);
    end
    fprintf('Saved %d samples/trial. %s Loading LFP tookd %.2f, saving %.2f\n',ns,outname,ltime,toc)
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
    if DATA.logfid
        fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s) to %s\n',datestr(now),nspk,sprintf(' %d',DATA.combineids),DATA.outname);
    end
    cfile = CombinerLst(DATA);
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
    DATA = ListExpts(DATA,DATA.Expts);
   % toc
    set(DATA.toplevel,'UserData',DATA);
  %  combine('listexpts');
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
        PlotTrodeXcorr(DATA,0);
    end
elseif strncmpi(name,'pcaplot',6)
    if length(varargin) & isnumeric(varargin{1})
        eid = varargin{1};
    else
        eid = DATA.currentexpt;
    end
    PlotPCAs(DATA,eid);
elseif strncmpi(name,'trodecorr',6)
    if length(varargin) & isnumeric(varargin{1})
        eid = varargin{1};
    else
        eid = DATA.currentexpt;
    end
    out= CalcTrodeCorrs(DATA,eid);
elseif strncmpi(name,'tetrodemov',6)
    ispk = DATA.Expts{DATA.currentexpt}.gui.spks;
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
        SetXYCluster(DATA,0,5,2,varargin{1});
        else
        SetXYCluster(DATA,0,4,2,varargin{1});
        end
    elseif ischar(varargin{1})
        a = strmatch(varargin{1},DATA.spkvarnames,'exact')
        b = strmatch(varargin{2},DATA.spkvarnames,'exact')
        if length(a) && length(b)
            SetXYCluster(DATA,0,4,2,a,b);
        end
    end
elseif strncmpi(name,'sortexpts',5)
    

elseif strncmpi(name,'l',7)
    lProbeClusters(DATA);
elseif strncmpi(name,'store',5)
    set(DATA.toplevel,'UserData',varargin{1});
elseif strncmpi(name,'trackcluster',8)
    id = get(DATA.elst,'value');
    sb = findobj(DATA.xyfig,'tag','Set+Next');
    for j = id:length(DATA.expid)
        DATA.exabsid = DATA.expid(j);
        set(DATA.elst,'value',j);
        combine('setexpt');
        DATA = get(DATA.toplevel,'UserData');
        C = OptimizeDprime(DATA);
        DATA.cluster{DATA.currentcluster,DATA.probe} = C;
        [DATA, newd, nc] = SetSpkCodes(DATA,DATA.spklist,DATA.probe,2);
        DrawClusters(DATA,DATA.cluster, 0);   
        set(sb,'backgroundcolor','r');
        drawnow;
    end
    
elseif strncmpi(name,'relist',6)
    DATA = ReadDir(DATA, DATA.name,'setprobe',DATA.probe);
    DATA = combine('listexps',DATA,'Tag',DATA.tag.top);
    eid = get(DATA.clst,'value');
    DATA = ListSubExpts(DATA,eid);

    set(DATA.toplevel,'UserData',DATA);
    if DATA.state.autoupdatelist
        names = get(DATA.clst,'String');
        eid = strmatch(DATA.Expts{end}.Header.expname, names,'exact');
        if length(eid) == 1
           set(DATA.clst,'value', eid);
           DATA = ListSubExpts(DATA,eid);
            nex = length(get(DATA.elst,'String'));
            set(DATA.elst,'value',nex);
            set(DATA.toplevel,'UserData',DATA);
        end
        WinFront(DATA.tag);
    end
elseif exist(name,'dir') % do dir before file, since dirs pass exist(name,'file')
    DATA.Expts = {};
    args = argon;
    if reindex
        args = {args{:}, 'relist'};
    end
    DATA.state.online = 1;
    DATA.state.autospool = 1;
    DATA.state.autoupdatelist = 1;
    DATA.state.showspikes = 1;
    DATA = ReadDir(DATA, name, varargin{:});
    if strncmp(DATA.filetype,'Grid',4)
    DATA = LoadClusters(DATA,ClusterFile(DATA),'allprobes');
    else
    DATA = LoadClusters(DATA,ClusterFile(DATA));
    end
    DATA.lastread = 0;
    DATA.lastsize = 0;
    if ~isfield(DATA,'timerobj')
            DATA.timerobj = timer('TimerFcn',@timerfna, 'Period', 1.0,...
                'Tag',DATA.tag.top, 'ExecutionMode','FixedSpacing');
            DATA.lastread = 0;
    end
    set(DATA.toplevel,'UserData',DATA);
    SetGui(DATA);
elseif strncmpi(name,'newfile',6)
    DATA.datafilename = get(findobj(DATA.toplevel, 'Tag','FileName'),'string');
    combine(DATA.datafilename);
elseif strncmpi(name,'recombineorbw',12)
    ReClassify
    ReCombineAll(DATA,0,'orbw');
elseif strncmpi(name,'recombine',6)
    ReCombineAll(DATA,0);
elseif strncmpi(name,'reptrial',6)
    PlayOneTrial(DATA,DATA.currenttrial,0);
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
        TrackTemplates(DATA);
        DATA = get(DATA.toplevel,'UserData');
    end
    if length(varargin) & isnumeric(varargin{1})
        cl = varargin{1};
        first = 2;
    else
        cl = 2;
        first = 1;
    end
    PlotTemplateScores(DATA,cl,varargin{first:end});
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

    DATA = ReadFile(DATA, name, args{:});
    set(DATA.toplevel,'UserData',DATA);
    if plotallspikes
        PlotAllExptsAndProbes(DATA);
    end
    if recombine == 3
        ReCombineAll(DATA,0,'orbw');
    elseif recombine == 4
        ReCombineByName(DATA,recombinename,1);
    elseif recombine > 0
        ReCombineAll(DATA,0);
    end
    if plotall
        d = fileparts(DATA.datafilename);
        PlotAllProbes(d,'save','sptrig')
    end
    DATA = PrintComments(DATA);
    if strfind(name,'ic-169')
        CheckSaccades(DATA.Expts,'Z:/smr/icarus/169/ic169.em.mat');
    end
    if ismember(recombine,[2 3 4])
        combine('close');
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
            DATA = ReadFile(DATA, name, args{:});
            set(DATA.toplevel,'UserData',DATA);
            if recombine == 1
                ReCombineAll(DATA,0);
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
        combine(name,'Tag',DATA.tag.top);
    end
elseif firstcall
    fprintf('No file or directory %s\n',name);
    close(DATA.toplevel);
    return;
end

SetGui(DATA);
NotBusy(DATA);


function mytoc(tstart)
   took = now-tstart;
   fprintf('Took %.3f sec\n',took * 24 * 60 * 60);
       
   
   
function DATA = CalcOnlineShapes(DATA)
    ts = now;
    DATA = ClearExptClusters(DATA);
    DATA = LoadOnlineClusters(DATA, ClusterFile(DATA,'getonline'));
    DATA = ReClassifyAll(DATA,'mkmean');
    SaveSpikeShape(DATA, strrep(DATA.meanspkfile,'spks','ospk'));
    fprintf('Took %.2f sec\n',(now-ts)*(24*60*60));


function DATA = PlotXcorr(DATA, probes, type)
        nc = length(probes);
    GetFigure('CrossCorrelation');
    if nc > 2
    for j = 1:nc
        subplot(nc,nc,((j-1)*nc)+1);
        for k = 1:j
            subplot(nc,nc,((j-1)*nc) + k);
            if type == 3
                CalcXcorr(DATA, DATA.currentexpt,probes(j),probes(k));
            elseif type == 1
                CalcXcorrDC(DATA, DATA.currentexpt,probes(j),probes(k));
            else
                CalcXcorrV(DATA, DATA.currentexpt,probes(j),probes(k));
            end
            set(gca,'ytick',[],'xtick',[]);
            if k == 1
                ylabel(sprintf('P%d',probes(j)));
            end
            if j == nc
                xlabel(sprintf('P%d',probes(k)));
            end
        end
    end
    else
        if type == 3
            CalcXcorr(DATA, DATA.currentexpt,probes(1),probes(2));
        elseif type == 1
            CalcXcorrDC(DATA, DATA.currentexpt,probes(1),probes(2));
        else
            CalcXcorrV(DATA, DATA.currentexpt,probes(1),probes(2));
        end
    end

    
function DATA = ReadGridFile(DATA)
    
    idx = [];
   datdir = 'F:/Utah/dufus/';
   datdir = fileparts(DATA.name);
   [a,b,c,d] = fileparts(DATA.name);
   if DATA.state.online
       id = strfind(DATA.name,'/');
       if ~isempty(id)
           dname = DATA.name(id(end)+1:end);
       else
           dname = test;
       end
       if strfind(a,'Spike')
           datdir = ['C:/Utah/' dname;]
       else
           datdir = DATA.name;
       end
   end
   idxfile = [datdir '/FileIdx.mat'];
   relist = 0;
   if exist(idxfile,'file') & ~relist
       a = load(idxfile);
       idx = a.idx;
   end
   nidx = length(idx);
   for j = 1:length(DATA.Expts)
       starts(j) = DATA.Expts{j}.Header.CreationDate + DATA.Expts{j}.Header.Start./(10000 * 60 *60 *24);
       ebstimes{j} = DATA.Expts{j}.bstimes;
   end
   d = dir(datdir);
   newf = 0;
   for j = 1:length(d)
       if strfind(d(j).name,'.nev') %this has starttime and Dig Events
           nevfile = [datdir '/' d(j).name];
           if relist || isempty(idx) || isempty(strmatch(d(j).name,idx.names))
               nev = openNEV('read','noparse',nevfile);
               ts = nev.MetaTags.DateTimeRaw;
               t  = datenum(ts(1),ts(2),ts(4),ts(5),ts(6),ts(7));
               [a, b] = min(abs(t-starts));
               onoff = diff(int16(bitand(4,nev.Data.SerialDigitalIO.UnparsedData)));
               id = find(onoff > 0);
               bstimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
               id = find(onoff < 0);
               estimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
               ns = length(bstimes);
               if ns && length(estimes)
               nidx = nidx+1;
               newf = newf+1;
                              idx.names{nidx} = d(j).name;

               idx.starttime(nidx) = t;
               idx.expt(nidx) = b;
               idx.tdiff(nidx) = a;
               idx.bsdiff(nidx,1) = sum(abs(diff(bstimes)-diff(ebstimes{b}(1:ns))));
               idx.bsdiff(nidx,2) = length(ebstimes{b})-ns;
               idx.toff(nidx) = ebstimes{b}(1)-bstimes(1);
               end
           
           end
       end
   end
   if newf
       save(idxfile,'idx');
   end
   id = find(idx.expt == DATA.currentexpt);
   if length(id) > 1
       [a,b] = min(idx.toff(id));
       id = id(b);
   end
   nfile = [];
   if length(id) == 1
       nfile = [datdir '/' idx.names{id}];
   end
   if ~isempty(nfile)     
       DATA = ReadGridSpikes(DATA,nfile,'toff',idx.toff(id));
   end
        
        
function DATA = ReadFile(DATA, name, varargin)
            
args = {};
j = 1;
while j <= length(varargin)
    args = {args{:} varargin{j}};
    j =j+1;
end
    DATA.datafilename = name;
    DATA.state.online = 0;
    DATA.logfid = fopen(strrep(name,'.mat', '.log'),'a');
    if isfield(DATA,'toplevel')
    set(findobj(DATA.toplevel, 'Tag','FileName'),'string',name);
    end
    if strfind(name,'.txt') %% an online text file
        Expts = ReadOnlineTxt(name, args{:});
        Expts = CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
        DATA.state.online = 2;
        DATA.lastread = 0;
        DATA.lastsize = 0;
        if ~isfield(DATA,'timerobj')
            DATA.timerobj = timer('TimerFcn',@timerfn, 'Period', 1.0, ...
                'Tag',DATA.tag.top, 'ExecutionMode','FixedSpacing');
            DATA.lastread = 0;
        end
    else
        [Trials, Expts, All] = APlaySpkFile(name, args{:});
        if isfield(Trials,'DataType')
            DATA.filetype = Trials.DataType;
        end
        if strncmp(DATA.filetype,'Grid',4)
            np = sscanf(DATA.filetype,'GridData %d');
            if isempty(np)
                np = 96; % default
            end
            DATA.probelist = 1:np;
            DATA.probesource = 100;
            DATA.probenames = cellstr(int2str(DATA.probelist'));
            DATA.state.applylastcluster = 1;
            DATA.state.classifyallprobes = 1;
            DATA.probe = DATA.probelist(1);

            for j = 1:np
                DATA.probes(np).probe = j;
            end
            if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
                DATA.probe = DATA.probelist(1);
            end
            if np > 1
                DATA = AddMultiProbeGUI(DATA);
            end
            
            if isfield(Trials,'bstimes')
                for nexp = 1:length(Expts)
                id = find(Trials.bstimes > Expts{nexp}.Header.Start & Trials.bstimes < Expts{nexp}.Header.End);
                Expts{nexp}.bstimes = Trials.bstimes(id);
                end
            end

        end

        if isfield(Trials,'Probes') &  ~isempty(Trials.Probes)
%  If alreayd have multi-probe setut, do not change
            if length(DATA.probelist) <= 1 & length(Trials.Probes) > 1
                if sum(unique([Trials.Probes.probe]) < 100) > 1 % count # of recording chans
                    DATA.state.includeprobename = 1;
                end

                DATA.probes = Trials.Probes;
                DATA = AddMultiProbeGUI(DATA);
                if isfield(Trials,'Comments') &  isfield(Trials.Comments,'Peninfo')
                    id = strfind(Trials.Comments.Peninfo.trode,'Contact');
                    if length(id)
                        x = id(1);
                        id = strfind(Trials.Comments.Peninfo.trode(id:end),' ');
                        x = sscanf(Trials.Comments.Peninfo.trode(id+x:end),'%d');
                        Trials.Comments.Peninfo.probesep = x;
                    end
                end
            end
            if isfield(Trials,'Probes') & ~isempty(Trials.Probes)
                        
            DATA.probelist = unique([Trials.Probes.probe]);
            if isfield(Trials.Probes,'var')
            DATA.probevars = {Trials.Probes.var};
            end
            DATA.probenames = cellstr(int2str([DATA.probelist]'));
            DATA.probes = Trials.Probes;
            DATA.probesource = FindProbeSources(DATA);
            else
                DATA.probelist = 1;
                DATS.probesource = 1;
            end
            if isfield(Trials,'setprobe')
                id = find(DATA.probelist == Trials.setprobe(1));
                if isempty(id)
                    id = 1;
                end
            else
                id = 1;
            end
            DATA.probe = DATA.probelist(id(1));
            if isfield(DATA, 'toplevel')
            it = findobj(DATA.toplevel,'Tag','ProbeId');
            set(it,'string',DATA.probenames,'value',id(1));
            end
        else
            DATA.state.savedvdt = 0; %don't need this if not switching probes
        end
        if DATA.state.online
            DATA.cellfile = [name, '.cells.mat'];
        else
        DATA.cellfile = strrep(name,'.mat','.cells.mat');
        end
        if exist(DATA.cellfile,'file')
            load(DATA.cellfile);
            DATA.CellList = CellList;
            if exist('CellQuality','var')
                if size(CellQuality,2) < size(CellList,2)
                   nq = size(CellQuality,2);
                   nl = size(CellList,2);
                   CellQuality(:,nq+1:nl) = 0;
                end
                DATA.CellQuality = CellQuality;
            else
                DATA.CellQuality = zeros(size(DATA.CellList));
            end
            if length(DATA.CellQuality) < length(DATA.CellList)
                lq = length(DATA.CellQuality);
                ll = length(DATA.CellList);
                DATA.CellQuality(:,lq+1:ll) = 0;
            end
            if exist('CellListCluster','var')
                DATA.CellListCluster = CellListCluster;
            else
                DATA.CellListCluster = zeros(size(CellList));
                id = find(CellList > 0)
                DATA.CellListCluster(id) = 1;
            end
            if exist('Templates','var')
                DATA.Templates = Templates;
            end
            if exist('TemplateInfo','var')
                DATA.TemplateInfo = TemplateInfo;
            end
        end
        
    end
    DATA = ListExpts(DATA,Expts);
    Trialids = [];
    TrialTimes = [];
    for nexp = 1:length(Expts)
        newt = [Expts{nexp}.Trials.Trial];
        news = [];
        for j = 1:length(newt)
            news(j) = Expts{nexp}.Trials(j).Start(1);
        end
        Trialids = [Trialids newt];
        TrialTimes = [TrialTimes news];
        Expts{nexp}.gui.classified = 0;
        Expts{nexp}.gui.counted = 0;
        Expts{nexp}.gui.clustertype = 0;
        if isfield(Expts{nexp}.Header,'Nam1e')
            Expts{nexp}.Header.Name = Expts{nexp}.Header.Nam1e;
        end
    end
    if isempty(Expts)
        NotBusy(DATA);
        return;
    end
    DATA.Expts = Expts;
    if isempty(DATA.exabsid)
        DATA.exabsid = 1:length(DATA.Expts);
    end
    if DATA.state.online
        DATA.AllData = [];
    else
    if isfield(All,'AllClusters')
        DATA.AllClusters = All.AllClusters;
        All = rmfield(All,'AllClusters');
        DATA.state.nospikes = 1;
    end
    DATA.AllData = All;
    DATA.AllData.SpikeIdx = Trials.Spkid;
    DATA.AllData.Trialids = Trialids;
    DATA.AllData.TrialTimes = TrialTimes;
    DATA.AllData.pcs = [];
    if isfield(Trials,'Comments')
        DATA.Comments = Trials.Comments;
    end

    DATA = SetSpkLists(DATA);
% load saved cluster params
    
    if DATA.state.useonlineclusters
        DATA = LoadOnlineClusters(DATA, ClusterFile(DATA,'getonline'));
    end
    DATA = LoadClusters(DATA, ClusterFile(DATA));
    DATA.name = name;

    if isfield(DATA,'toplevel')
    DATA = ListSubExpts(DATA,1,'relist');
    end
    if length(DATA.probelist) > 1 && DATA.state.NeedAllProbeFile
    DATA = LoadClusters(DATA, ClusterFile(DATA,'allprobes'),'allprobes');
    end
    cfile = CombinerLst(DATA);
    if exist(cfile,'file')
        load(cfile);
        DATA.combines = combines;
    end
    if isfield(Trials,'setprobe') & length(Trials.setprobe) > 1
        DATA.plot.useprobe = zeros(size(DATA.probelist));
        DATA.plot.useprobe(Trials.setprobe) = 1;
            DATA.currentexpt = 1;
            DATA = LoadAllProbeSpikes(DATA, 0, 'select');
    else
%        DATA.plot.useprobe = zeros(size(DATA.probelist));
%        DATA.plot.useprobe(DATA.probe) = 1;
    end
    
    
    end

    function probesource = FindProbeSources(DATA)
        if ~isfield(DATA.probes,'source')
            probesource = ones(size(DATA.probelist));
        else
            for j = 1:length(DATA.probelist)
                id = find([DATA.probes.probe] == DATA.probelist(j));
                probesource(j) = DATA.probes(id(1)).source;
            end
        end

function PrintLogData(DATA, type)
       id = find(DATA.log.types == type);
       if isempty(id)
           return;
       end
       allexp = unique(cat(1,DATA.log.savexps{id}));
       [done, b] = max(DATA.log.savedates(id));
       last = id(b);
       label = DATA.log.expnames{last};
       fprintf('%s %s: Expts ',label,datestr(done));
       fprintf('%d ',DATA.log.savexps{last});
       x = setdiff(allexp, DATA.log.savexps{last});
       if length(x)
           fprintf('Other '); 
           fprintf('%d ',x);
       end
       fprintf('\n');
       

function Expt = ParseComment(Expt,s)

if strmatch('cm=back',s)
    sscanf(s,'cm=back=%s,wi=%f,hi=%f');
    c = findstr(s,',');
    Expt.Stimvals.backstim = s(9:c(1)-1);
    a = sscanf(s(c(1):end),',wi=%f,hi=%f,ce=%f,xo=%f,yo=%f,dx=%f,co=%f');
    Expt.Stimvals.backwi = a(1);
    Expt.Stimvals.backhi = a(2);
    Expt.Stimvals.backxo = a(4);
    Expt.Stimvals.backyo = a(5);
    Expt.Stimvals.backdx = a(6);
    Expt.Stimvals.backco = a(7);
end
        
       
function s = ExptComments(Expt)
    if isfield(Expt,'Comments') & ~isempty(Expt.Comments.text);
        s = Expt.Comments.text;
        bid = strmatch('cm=back=',s);
        gid = setdiff(1:length(s),bid);
        for j = gid
            fprintf('%s\n',s{j});
        end
        for j = 1:length(bid)
%        ParseComment(Expt, s{bid(j)});
        end
    else
        s = {};
    end

function DATA = PrintComments(DATA)
    
    if isfield(DATA.Expts{1}.Header,'rfstr');
        fprintf('%s\n',DATA.Expts{1}.Header.rfstr);
    end
    if isfield(DATA,'Comments')
        if isfield(DATA.Comments,'Peninfo')
            fprintf('%s\n',DATA.Comments.Peninfo.trode);
            a= InterpretLine(DATA.Comments.Peninfo.trode);
            DATA.Header.trode = DATA.Comments.Peninfo.trode;
            CopyFields(DATA.Header,a);
        end
        for j = 1:length(DATA.Comments.text)
            if isempty(strmatch('cm=back=',DATA.Comments.text{j})) & ...
                isempty(strfind(DATA.Comments.text{j},'Tube Protrudes'))
                id = find(DATA.AllData.TrialTimes < DATA.Comments.times(j));
                if isempty(id)
                    id = 1;
                end
            fprintf('%.1f(Trial%.0f): %s\n',DATA.Comments.times(j)./10000,...
                DATA.AllData.Trialids(id(end)),...
                DATA.Comments.text{j});
            end
        end
    end
    cmfile = strrep(DATA.datafilename,'.mat','.txt');
    fid = fopen(cmfile,'r');
    if fid > 0 
        fprintf('Offline Comments:\n');
        a = textscan(fid,'%s','delimiter','\n');
        for j = 1:length(a{1})
            fprintf('%s\n',a{1}{j});
        end
        fclose(fid);
    end


function NotBusy(DATA)

[a,b,c] = fileparts(DATA.datafilename);
if isfield(DATA,'toplevel') && ishandle(DATA.toplevel)
set(DATA.toplevel,'Name',[DATA.tag.top ':' b]);
end

function ShowLFPPlot(DATA)
    b = [];

    GetFigure(DATA.tag.lfpplot);
    a = [];
    if ismember(DATA.plot.lfpplot,[4 5 6 7 8 9 10 11 12 13 14])
        if DATA.LFP.Header.rc
            if DATA.LFP.Stimvals.n2 > 4 & DATA.LFP.Stimvals.nt > 4
            [a,b] = PlotRevCorAny(DATA.LFP,'lfp',DATA.probe,'nmin',20,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb,'twoslice');
            else
            [a,b] = PlotRevCorAny(DATA.LFP,'lfp',DATA.probe,'nmin',20,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb);
            end
            a.Header = b.Header;
            GetFigure(DATA.tag.lfpplot);
        else
            DATA.probelist=1:24;
            a = PlotLFP(DATA.LFP);
            a.Header = DATA.LFP.Header;
            a.Header.exptype = regexprep(a.Header.expname,'.*\.','');
        end
    end
    hold off;
    probelist = DATA.probelist;
    probelist=1:24;
    if DATA.plot.lfpplot == 1
        if DATA.LFP.Header.rc > 0
            if isempty(b) 
                [a,b] = PlotRevCorAny(DATA.LFP,'lfp',DATA.probe,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb,'nmin',DATA.plot.nminrc);
                a.Header = b.Header;
                GetFigure(DATA.tag.lfpplot);
            end
            PlotAllProbes(a);
        elseif strmatch('square.co',DATA.LFP.Header.expname)
            a = PlotMLFP(DATA.LFP,'image');
        else
            a = PlotLFP(DATA.LFP);
            a.Header = DATA.LFP.Header;
            a.Header.exptype = regexprep(a.Header.expname,'.*\.','');
            subplot(1,1,1);
            PlotAllProbes(a,'LFPTrial','probes',probelist);
        end     
    end
    issdf = isfield(a,'sdfs');
    if DATA.plot.lfpplot == 3
        PlotMLFP(DATA.LFP,'image','probes',probelist);
    elseif DATA.plot.lfpplot == 2
        PlotMLFP(DATA.LFP,'stack','probes',probelist);
    elseif DATA.plot.lfpplot == 4
        PlotAllProbes(a,'blank','probes',probelist);
    elseif DATA.plot.lfpplot == 5

    elseif DATA.plot.lfpplot == 6
        PlotAllProbes(a,'probes',probelist);
    elseif DATA.plot.lfpplot == 7 & issdf
        PlotAllProbes(a,'onestim','xvals',1:length(a.sdfs.x(:,1)),'yvals',1:length(a.sdfs.x(1,:)),'probes',probelist);
    elseif DATA.plot.lfpplot == 8
        PlotAllProbes(a,'monoc','probes',probelist);
    elseif DATA.plot.lfpplot == 9
        PlotAllProbes(a,'LFPEig','probes',probelist);
    elseif DATA.plot.lfpplot == 10
        PlotAllProbes(a,'stimvar','probes',probelist);
    elseif DATA.plot.lfpplot == 11
        PlotAllProbes(a,'varblank','probes',probelist);
    elseif DATA.plot.lfpplot == 12
        PlotAllProbes(a,'frameresp','probes',probelist);
    elseif DATA.plot.lfpplot == 13
        PlotAllProbes(a,'LFPTrial','probes',probelist);
    elseif DATA.plot.lfpplot == 14
        PlotAllProbes(a,'CSD','probes',probelist);
    end

function str =ShowString(DATA, Expt)
    str = '';
    if isfield(Expt.Stimvals,'rb') & Expt.Stimvals.rb ~= 0
        str = sprintf('rb%.0f',Expt.Stimvals.rb);
    end
    fn = fields(DATA.show);
    for k = 1:length(fn)
        if isempty(strmatch(fn{k},'times') )
            x = GetEval(Expt,fn{k});
        if strmatch(fn{k},'Fr') & x > 1
         str = [str sprintf(' %s=%.1f',fn{k},x)];
        elseif DATA.show.(fn{k})
            if strmatch(fn{k},{'jx' 'pi'})
                str = [str sprintf(' %s=%.3f',fn{k},x)];
            elseif strmatch(fn{k},{'or' 'sl' 'se' })
                str = [str sprintf(' %s=%.0f',fn{k},x)];
            else
                str = [str sprintf(' %s=%.2f',fn{k},x)];
            end
        end
        end
    end

        
function Expts = CountTxtSpikes(Expts, probe, cl)

if length(cl) > 2
    cl = cl(1:2);
end
for j = 1:length(Expts)
    if isfield(Expts{j},'Trials')
    for k = 1:length(Expts{j}.Trials)

        if size(Expts{j}.Trials(k).AllSpikes,2) < cl(1)+1
            Expts{j}.Trials(k).Spikes = [];
        elseif length(cl) > 1 && size(Expts{j}.Trials(k).AllSpikes,2) > cl(end)
            Expts{j}.Trials(k).Spikes = union(Expts{j}.Trials(k).AllSpikes{probe,cl+1});
        else
            Expts{j}.Trials(k).Spikes = Expts{j}.Trials(k).AllSpikes{probe,cl(1)+1};
        end
    end
    else
        fprintf('Expt %d no trials\n');
    end
    Expts{j}.gui.counted = 1;
end

function Spks=GetSpikeStruct(DATA)
        
    if isempty(DATA.AllData.Spikes)
        Spks = DATA.AllSpikes{DATA.probe};
    else
        Spks = DATA.AllData.Spikes;
    end
    
function DATA = LoadClusters(DATA, cfile, varargin)
    clearold = 1;
    allprobes = 0;
    nprobes = 1;
    probe = DATA.probe;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{1},'allprobes',4)
            allprobes = 1;
        elseif strncmpi(varargin{j},'noclear',5)
            clearold = 0;
        elseif strncmpi(varargin{j},'probe',5)
            j =j+1;
            probe = varargin{j};
        end
        j = j+1;
    end
    if ~isfield(DATA.AllData,'Spikes')
        return;
    end
    if isempty(DATA.AllData.Spikes)
        Spks = DATA.AllSpikes{probe};
    else
        Spks = DATA.AllData.Spikes;
    end
    if isempty(Spks)
        return;
    end
    if exist(cfile,'file')
        excludelist = [];
        clid = [];
        classified = 0;
        load(cfile);
        DATA.state.recut = 1;
        firstspk = length(Spks.times);
        lastspk = 0;
        if size(clid,1) == size(Spks.codes,1);
            Spks.codes(:,2) = clid;
            last = max(find(clid > 0));
            maxclasst = Spks.times(last);
            classified = 1;
        end
        for j = 1:min([length(excludelist) length(DATA.Expts)]);
            if max(excludelist{j}) <= length(DATA.Expts{j}.Trials)
            lst = -abs([DATA.Expts{j}.Trials(excludelist{j}).Trial]);
            if lst
                for k = 1:length(lst)
                    [DATA.Expts{j}.Trials(excludelist{j}(k)).Trial] = lst(k);
                end
            end

            end
        end
        p = probe;
%should chnage this to use AllClusters, then find correct expt
%especially for loading clusters define online, where Expt order can
%change.
        nx=0;
        for j = 1:length(DATA.Expts)
            eid(j,1) = DATA.Expts{j}.Trials(1).id;
            eid(j,2) = DATA.Expts{j}.Trials(end).id;
        end
        nset = 0;
        for j = 1:length(AllClusters)
            %only load the data for the current probe if this is a single probe file
            if ~isempty(AllClusters{j})
                if isfield(AllClusters{j},'ids')
                    exid = find(eid(:,1) < AllClusters{j}.ids(2) & ...
                        eid(:,2) > AllClusters{j}.ids(1));
                    if length(exid) == 1
                        nx = exid;
                    elseif length(exid) > 1
                        overlap = [];
                        for k = 1:length(exid)
                            overlap(k) = sum(ismember([eid(exid(k),1):eid(exid(k),2)],[AllClusters{j}.ids(1):AllClusters{j}.ids(2)]));
                        end
                        [a,k] = max(overlap);
                        nx = exid(k);                        
                    end
                elseif j <= length(DATA.Expts)
                    nx = j;
                end
                if nx > 0 && isfield(AllClusters{j},'Cluster')
                sz = size(AllClusters{j}.Cluster);
                if allprobes
                    DATA.Expts{nx}.Cluster = AllClusters{j}.Cluster;
                    % some old single probe files have AllClusters the other way around
                elseif sz(1) == 1 && sz(2) > 1 && length(DATA.probes) == 1
                    for k = 1:sz(2)
                        DATA.Expts{nx}.Cluster{k,1} = AllClusters{j}.Cluster{k};
                    end
                elseif sz(1) >= 0 && sz(2) >= 1 && length(DATA.probes) == 1
                    DATA.Expts{nx}.Cluster(1:sz(1),1) = AllClusters{j}.Cluster(:,sz(2));
                elseif sz(1) > 0 && sz(2) >= p
                    DATA.Expts{nx}.Cluster(1:sz(1),p) = AllClusters{j}.Cluster(:,p);
                    for k = (sz(1)+1):size(DATA.Expts{nx}.Cluster,1)
                    DATA.Expts{nx}.Cluster{k,p} = [];
                    end
%If Allclsuters{j).Cluster is empty, make sure DATA.Expts is too.
                elseif sz(1) == 0 && sz(2) >= p  & isfield(DATA.Expts{nx},'Cluster')
                    for k = 1:size(DATA.Expts{nx}.Cluster,1)
                    DATA.Expts{nx}.Cluster{k,p} = [];
                    end
                end
                if exist('clustertypes','var') %%saved type
                    DATA.Expts{nx}.gui.clustertype = clustertypes(j);
                else
                    DATA.Expts{nx}.gui.clustertype = 1;
                end
                firstspk = length(Spks.times);
                lastspk = 0;
                nset = 0;
                if ~isfield(DATA.cluster{1,1},'Arange')
                    DATA.cluster{1,1}.Arange = [6:8];
                    DATA.cluster{1,1}.Brange = [11:20];
                    DATA.cluster{1,1}.Erange = [1:100];
                end
                end

                %only do this for one probe, as AllData.Spikes only has one probes data
                p = probe;
                if nx > 0 && isfield(DATA.Expts{nx},'Cluster')
                    nprobes = size(DATA.Expts{nx}.Cluster,2);

                    for k = 1:size(DATA.Expts{nx}.Cluster,1)
                        if p > nprobes
                            DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                            DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                            nprobes = p;
                        end
                        if isempty(DATA.Expts{nx}.Cluster{k,p}) | ~isfield(DATA.Expts{nx}.Cluster{k,p},'firstspk')
                            DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                            DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                        else
                            nset = nset + 1;
                        end
                        if DATA.Expts{nx}.Cluster{k,p}.firstspk == 0
                            DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                        end
                        if DATA.Expts{nx}.Cluster{k,p}.lastspk == 0
                            DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                        end
                        if DATA.Expts{nx}.Cluster{k,p}.firstspk < firstspk
                            ts = Spks.times(DATA.Expts{nx}.Cluster{k,p}.firstspk);
                            if ts  > DATA.Expts{nx}.Header.trange(1) && ts < DATA.Expts{nx}.Header.trange(2)
                                firstspk = DATA.Expts{nx}.Cluster{k,p}.firstspk;
                            end
                        end
                        if DATA.Expts{nx}.Cluster{k,p}.lastspk > lastspk && ...
                                DATA.Expts{nx}.Cluster{k,p}.lastspk <= length(Spks.times)
                            ls = Spks.times(DATA.Expts{nx}.Cluster{k,p}.lastspk);
                            if ls  > DATA.Expts{nx}.Header.trange(1) && ls < DATA.Expts{nx}.Header.trange(2)
                                lastspk = DATA.Expts{nx}.Cluster{k,p}.lastspk;
                            end
                        end
                        if DATA.Expts{nx}.gui.clustertype == 2 % online cluster - can't use spk coun
                            DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                            DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                        end
                        if ~isfield(DATA.Expts{nx}.Cluster{k,p},'params')
                            DATA.Expts{nx}.Cluster{k,p}.params = [1 2];
                        end
                        if ~isfield(DATA.Expts{nx}.Cluster{k,p},'Arange') & isfield(DATA.cluster{1,1},'Arange')
                            DATA.Expts{nx}.Cluster{k,p}.Arange = DATA.cluster{1,1}.Arange;
                            DATA.Expts{nx}.Cluster{k,p}.Brange = DATA.cluster{1,1}.Brange;
                            DATA.Expts{nx}.Cluster{k,p}.Erange = DATA.cluster{1,1}.Erange;
                        end

                    end
                end

                if lastspk < firstspk %nevet found these in cluster file
                    firstspk = NaN;
                    lastspk = NaN;
                end
                if lastspk > size(Spks.codes,1)
                    lastspk = size(Spks.codes,1);
                end
                % with multiple probes, empty clsuters can arise because clusters have beed set
                % on other probes.
                if nset == 0 && nprobes == 1 % actively set no clusters. Only works if one probe
                    DATA.Expts{nx}.gui.classified = 1;
                    DATA.Expts{nx}.gui.ncluster = 0;
                elseif isnan(firstspk)
                    DATA.Expts{nx}.gui.classified = 0;
                    DATA.Expts{nx}.gui.ncluster = 0;
                elseif classified && lastspk > firstspk && max(Spks.codes(firstspk:lastspk,2)) > 0
                    DATA.Expts{nx}.gui.classified = classified;
                    DATA.Expts{nx}.gui.ncluster = nset;
                else
                    DATA.Expts{nx}.gui.classified = 0;
                    DATA.Expts{nx}.gui.ncluster = 0;
                end
            end
        end
    elseif clearold
        p = probe;
        for j = 1:length(DATA.Expts);
            if isfield(DATA.Expts{j},'Cluster')
            for k = 1:size(DATA.Expts{j}.Cluster,1)
            DATA.Expts{j}.Cluster{k,p} = {};
            end
            end
        end
    end

    if isempty(DATA.AllData.Spikes)
        DATA.AllSpikes{probe}.codes = Spks.codes;
    else
        DATA.AllData.Spikes.codes = Spks.codes;
    end
function DATA = ClearExptClusters(DATA)

    for j = 1:length(DATA.Expts)
        DATA.Expts{j}.Cluster = {};
    end
    
    
function DATA = LoadOnlineClusters(DATA, cfile)
    if length(DATA.probelist) > 4
        minq = 5;
    else
        minq = 0;
    end
    if exist(cfile,'file')
        load(cfile);
        DATA.state.recut = 1;
        p = DATA.probe;
        for j = 1:length(DATA.Expts)
            eid(j,1) = DATA.Expts{j}.Trials(1).id;
            eid(j,2) = DATA.Expts{j}.Trials(end).id;
        end
        for j = 1:length(AllClusters);
            if ~isempty(AllClusters{j})
                if isfield(AllClusters{j},'ids')
                    exid = find(eid(:,1) < AllClusters{j}.ids(2) & ...
                        eid(:,2) > AllClusters{j}.ids(1));
                    if length(exid) == 1 & isfield(AllClusters{j},'Cluster')
                        k = exid;
                        DATA.Expts{k}.OnlineCluster = AllClusters{j}.Cluster;
                        for ic = 1:length(AllClusters{j}.Cluster(:))
%only copy online cuts that were not automatic
%Feb 20011 Change to apply test to the quality assigned online
%don't copy assigned quality - it needs resetting when offline clsuteres
%are cut
                            C= AllClusters{j}.Cluster{ic};
                            if ~isfield(C,'autocut') || C.autocut == 0 
                                if (isfield(C,'quality') && C.quality >= minq) || minq == 0
                                    C.quality = 0;
                                    DATA.Expts{k}.Cluster{ic} = C;
                                    DATA.Expts{k}.Cluster{ic}.touched = -1;
                                else
                                    DATA.Expts{k}.OnlineCluster{ic}.touched = -1;
                                end
                            end
                            DATA.Expts{k}.OnlineCluster{ic}.touched = -1;
                        end
                        DATA.Expts{k}.gui.clustertype = 2;
                        if ~isfield(DATA.Expts{k},'Cluster')
                            DATA.Expts{k}.Cluster = {};
                        end
%
%Fill in fields that might be missing from old cluster files
%remove firstpk, lastspk refs from the online file - these spike id #s will
%not necessarily match the saved fil
                        for m = 1:size(DATA.Expts{k}.Cluster,1)
                            DATA.Expts{k}.Cluster{m,p}.firstspk = NaN;
                            DATA.Expts{k}.Cluster{m,p}.lastspk = NaN;
                            if ~isfield(DATA.Expts{k}.Cluster{m,p},'params')
                                DATA.Expts{k}.Cluster{m,p}.params = [1 2];
                            end
                        end
                        DATA.Expts{k}.gui.classified = 0; %loaded def, but not classified
                    end
                end
            end
        end
    end


function [DATA, Stimvals, needv, retval] = CheckCombine(DATA, interactive)
    
    if nargin < 2
        interactive = 1;
    end
    retval = 1;
    err = '';
    id = get(DATA.elst,'value');
    for j = 1:length(id)
        c.stims(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'st');
        c.bws(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'ob');
        c.mes(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'me');
        c.tfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'tf');
        c.sfs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'sf');
        c.ors(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'or');
        c.sls(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'sl');
        c.ces(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'ce');
        c.xos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'xo');
        c.yos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'yo');
        c.fxs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'fx');
        c.fys(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'fy');
        c.jxs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'jx');
        c.cos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'co');
        c.Bcs(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'Bc');
        c.bos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'bo');
        c.backxos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'backxo');
        c.backyos(j) = GetEval(DATA.Expts{DATA.expid(id(j))},'backyo');
        ets{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.et;
        e2s{j} = DATA.Expts{DATA.expid(id(j))}.Stimvals.e2;
        Stimvals = DATA.Expts{DATA.expid(id(j))}.Stimvals;
    end

    
    
 % make c.(xx) a list of unigue values of xx
 % if there are > 1 of these, it needs to be a field in Trials
    f = fields(c);
    needv = {};
    n = 1;
    for j = 1:length(f)
        vals = unique(c.(f{j}));
        vals = vals(find(~isnan(vals)));
        c.(f{j}) = vals;
        if length(vals) > 1
            needv{n} = f{j}(1:end-1); %remove the 's'
            n = n+1;
        end
    end

    Stimvals.BlockedStim = 0;
    for j = 1:length(e2s)
            if isempty(e2s{j}) 
            e2s{j} = 'e0';
            end
            if isempty(ets{j}) 
            ets{j} = 'e0';
            end
    end
    blank = strmatch('none',c.stims);
    nstims = length(c.stims) - length(blank);
    blank = strmatch('e0',unique(e2s));
    e2lst = unique(e2s);
    e1lst = unique(ets);
    ne2 = length(e2lst) - length(blank);
    bws = unique(c.bws);
    bws = bws(find(~isnan(bws)));
    xos = unique(c.xos);
    yos = unique(c.yos);
%If there are blocks with different values for some paramerter, combine them
%and be sufe to fill trials, then set to expt 2. 
%Careful - this forces combine to generate a new name. 
    if length(bws) > 1 & ne2 == 0 & Stimvals.st== 22 %bw only means anything for stim image
        err = [err sprintf('%d Bandwidths in %s',length(unique(bws)),DATA.outname) sprintf('%.2f ',unique(bws))];
        for j = id(1:end)
            DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'ob');
            Stimvals.e2 = 'ob';
            Stimvals.BlockedStim = 1;
        end
    elseif ne2 == 1 & length(blank) > 0 &  strmatch('me',e2lst) %Expt with me varying + one without
        for j = id(1:end)
            DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},'me');
            Stimvals.e2 = 'me';
            Stimvals.BlockedStim = 1;
        end        
    end
    if length(c.sfs) >1 & isempty(strmatch('sf',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d SFS in %s',length(c.sfs),DATA.outname) sprintf('%.2f ',c.sfs)];
    end
    if length(c.ors) > 1 & isempty(strmatch('or',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d Oris in %s',length(c.ors),DATA.outname) sprintf('%.2f ',c.ors)];
    end
    if length(c.sls) > 1 & isempty(strmatch('sl',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d Seed Loops in %s',length(c.sls),DATA.outname) sprintf('%.2f ',c.sls)];
    end
    if length(c.jxs) > 1 & isempty(strmatch('jx',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d jx vals in %s',length(c.jxs),DATA.outname) sprintf('%.3f ',c.jxs)];
    end
    if range(xos) > 0.5 & isempty(strmatch('Op',{e1lst{:} e2lst{:}})) & isempty(strmatch('Pp',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d X pos (%.1f deg) in %s',length(c.xos),range(xos),DATA.outname) sprintf('%.2f ',c.xos)];
    end
    
    if length(c.bos) > 1 && (c.stims ~= 2 || sum(c.sls) > 0) && isempty(strmatch('b',{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d BackOr in %s',length(c.bos),DATA.outname) sprintf('%.3f ',c.bos)];
    end
    f = {'cos' 'backxos' 'backyos' 'Bcs'};
    ex = {'co' 'backxo' 'backyo' 'Bc'};
    names = {'Contrast' 'BackGr Ori' 'Back Xpos' 'Back Ypos' 'Back Contrast'};
    for j = 1:length(f)
    if length(c.(f{j})) > 1 & isempty(strmatch(ex{j},{e1lst{:} e2lst{:}}))
      err = [err sprintf('%d %s in %s',length(c.(f{j})),names{j},DATA.outname) sprintf('%.3f ',c.(f{j}))];
    end
    end
    if length(err) > 1
        if interactive 
        a = questdlg(err,'Combine these Files','Cancel','OK','OK');
        if strmatch(a,'Cancel')
            retval = -1;
        end
        else
            fprintf('%s\n',err);
        end
    end
 
    
 function AutoClusters = SaveAutoClusters(DATA)
 
     cfile = ClusterFile(DATA,'auto');
     if exist(cfile,'file')
         load(cfile);
     end
     for j = 1:length(DATA.Expts)
         for k = 1:length(DATA.probelist)
             p = DATA.probelist(k);
             if isfield(DATA.Expts{j}.Cluster{1,p},'autocut') && DATA.Expts{j}.Cluster{1,p}.autocut
                 AutoClusters{j}{1,p} = DATA.Expts{j}.Cluster{1,p};
             end
         end
     end
     save(cfile,'AutoClusters');
     
 function DATA = AutoCut(DATA,eid, expno, varargin)
 
     minimize = 0;
     plotmode = 1;
     show = 2;
     j = 1;
     mode = DATA.plot.autoclustermode;
     while j <= length(varargin)
         if strncmpi(varargin{j},'cluster',4)
             mode = 1;
         elseif strncmpi(varargin{j},'minimize',4)
             minimize = 1;
         elseif strncmpi(varargin{j},'noplot',4)
             plotmode = 0;
             show = 0;
         end
         j = j+1;
     end
     
     
     if DATA.state.online
         DATA = LoadSpikes(DATA,eid);
     end
     [DATA, ispk] = SetExptSpikes(DATA, eid, 'setrange');

     if length(ispk) < 4
         dprime = 0;
     end
     if length(ispk) > 3
         DATA = CalcClusterVars(DATA, DATA.Expts{eid}.gui.spkrange);

         if mode == 1
             if isfield(DATA,'AllClusters')
                 Cy = DATA.AllClusters(DATA.probe).cy;
                 Cx = DATA.AllClusters(DATA.probe).cx;
             else
                 Cy = DATA.Spikes.cy;
                 Cx = DATA.Spikes.cx;
             end
             [x,y,z] = CalcDensity(DATA, ispk, 3);
             xv = x(1,:);
             yv = y(:,1);
             % find local maxima in smoothed 2-D density
             sd = 5;
             np = 10;

             nl = 1;
             stepscale = 1;
             while (np > 2 || np == 1) & stepscale > 0.1 & sd < 25 & sd > 0.25
                 if sd > 15
                     [a,b,K] = gauss2d(sd,[-60:60]);
                 else
                     [a,b,K] = gauss2d(sd,[-30:30]);
                 end
                 K = K./sum(K(:));
                 Z = conv2(z,K,'same');
                 dy = diff(Z,[],2);
                 dx = diff(Z);
                 yz = diff(sign(dy),[],2);
                 xz = diff(sign(dx));
                 lastnp = np;
                 pid = find(yz(2:199,:) < 0 & xz(1:198,2:199) < 0);
                 if length(ispk) > 5000
                     zthresh = max(Z(:))./40;
                 elseif length(pid) > 4
                     zthresh = max(Z(:))/10;
                 elseif length(pid) > 2
                     zthresh = max(Z(:))./20;
                 else
                     zthresh = 0;
                 end
                 pid = find(yz(2:199,:) < 0 & xz(1:198,2:199) < 0 & Z(2:199,2:199) > zthresh);

                 np = length(pid);
                 if np == 2
                     [py,  px] = ind2sub([198 198], pid);
                     d = abs((px(2)-px(1) + i * (py(2)-py(1))));
                     if d < 30
                         np = 1;
                     end
                 end
                 nps(nl) = np;
                 sds(nl) = sd;
                 zts(nl) = zthresh;
                 if nl ==1 && np > 10
                     stepscale = 4;
                 end
                 if np == 1 && sd < 50 && sd > 1
                     sd = sd - 0.5 * stepscale;
                 elseif np == 1 && sd  < 1
                     sd = sd * 0.8;
                 elseif np == 1
                     sd = sd * 0.8;
                 elseif nl > 1 && nps(nl-1) == 1
                     stepscale = stepscale  * 0.25;
                     sd = sd + stepscale;
                 else
                     sd = sd +stepscale;
                 end
                 nl = nl+1;
             end
             if np == 0
                 [a,pid] = max(Z(2:199,2:199));
                 np = 1;
             end
             [py,  px] = ind2sub([198 198], pid);
             px = px+1;
             py = py+1;
             [pd, peaki] = max(abs(px +i*py)); %find most distant local maximum
             peaka = [px(peaki) py(peaki)];

             % now find highest peak that is not this one
             [id, idx] = setdiff(pid,pid(peaki));
             zid = sub2ind([200 200],py(idx),px(idx));
             [a,b] = max(Z(zid));
             if np > 1
                 % This didn't help...
                 %             [a,c,K] = gauss2d(sd/2,[-30:30]);
                 %             K = K./sum(K(:));
                 %             Zf = conv2(z,K,'same');
                 peakb = [px(idx(b)) py(idx(b))];
                 xdi = linspace(x(1,peaka(1)),x(1,peakb(1)));
                 ydi = linspace(y(peaka(2),1),y(peakb(2),1));
                 dd = interp2(x,y,Z,xdi,ydi);
                 [dmin,b] = min(dd);
                 minpt = [xdi(b) ydi(b)];
                 zmax = Z(peaka(2),peaka(1));

                 d = sqrt(sum((peaka-peakb).^2));
             else
                 zmax = max(Z(:));
                 peakb = [];
                 d = 0;
             end
             if length(ispk) < 500
                 minzt = 0.1;
             elseif length(ispk) < 2000
                 minzt = 0.05;
             else
                 minzt = 0.02;
             end
             if np == 0
             elseif np == 1 || d < 30 || dmin/zmax > 0.8 || zmax/max(Z(:)) < minzt % too close treat as a single blob
                 if np > 1 %% peak is not for real
                     [a,id] = max(Z(:));
                     [peaka(2),  peaka(1)] = ind2sub([200 200], id);
                     if dmin/zmax > 0.8
                         dmin = zmax/5;
                     else
                         dmin = zmax/20;
                     end
                 else
                     dmin = zmax/50;
                 end
                 twopeaks = 0;
                 xu = max(Cx(ispk));
                 xl = peaka(1);
                 xc = (xv(xl)+xu)/2;
                 xr = (xu-xv(xl))/2;
                 id = find(Z(peaka(2):end,peaka(1)) < dmin);
                 if isempty(id)
                     id = find(Cx(ispk) > xv(xl));
                     yu = max(Cy(ispk(id)));
                 else
                     yu = yv(id(1)+peaka(2)-1);
                 end
                 yr = yu-yv(peaka(2));
                 yc = yv(peaka(2));
             else
                 twopeaks = 1;
                 %dmin is the low density value between teh cluster peak and the rest
                 %should find lower values in the outward directions, but may not - if the
                 %cluster is very good, dmin can be low.
                 id = find(Z(peaka(2),1:peaka(1)) > dmin);
                 if isempty(id)
                     id = find(Z(peaka(2),1:peaka(1)) < zmax/20);
                 end
                 if isempty(id)
                     xl = 1;
                 else
                     xl = id(1);
                 end
                 id = find(Z(peaka(2),peaka(1):end) < dmin);
                 if isempty(id)
                     xu = length(xv);
                 else
                     xu = id(1)+peaka(1)-1;
                 end
                 id = find(Z(1:peaka(2),peaka(1)) > dmin);
                 if isempty(id) % check for local minima
                     yl = 1;
                 else
                     yl = id(1);
                 end
                 id = find(Z(peaka(2):end,peaka(1)) < dmin);
                 if isempty(id)
                     yu = length(yv);
                 else
                     yu = id(1)+peaka(2)-1;
                 end
                 k = (xv(xu)-xv(xl))/(yv(yu)-yv(yl)); %ratio of radii
                 yr = sqrt((ydi(b)-ydi(1)).^2 + (xdi(b)-xdi(1)).^2/k^2);
                 xr = yr * k;

                 id = find(Z(peaka(2), peaka(1):end) < zmax/20);
                 if isempty(id)
                     xu = length(xv);
                 else
                     xu = peaka(1)+id(1);
                     if xu > length(xv)
                         xu = length(xv);
                     end
                 end
                 id = find(Z(peaka(2):end,peaka(1)) < zmax/20);
                 if isempty(id)
                     yu = max(Cy(ispk));
                 else
                     yu = yv(peaka(2)+id(1)-1);
                 end
                 k =  (xv(xu)-xdi(1))/xr;
                 xc = xdi(1) + xr * (k-1)/2;
                 xr = xr + xr *(k-1)/2;
                 k = (yu-ydi(1))/yr;
                 yc = ydi(1) + yr * (k-1)/2;
                 yr = yr + yr *(k-1)/2;
             end
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             C.angle = 0;
             C.firstspk = NaN;
             C.lastspk = NaN;
             C.autocut = 2;
             C.params(1) = DATA.plot.clusterX;
             C.params(2) = DATA.plot.clusterY;
             C.Arange = DATA.clusterArange;
             C.Brange = DATA.clusterBrange;
             C.Erange = DATA.clusterErange;
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA.Expts{eid}.gui.classified = 3;
             if plotmode
                 GetFigure('AutoCut');
                 subplot(2,1,1);
                 hold off;
                 pcolor(x,y,Z); shading('interp');
                 hold on;
                 plot(xv(peaka(1)),yv(peaka(2)),'ok');
                 if ~isempty(peakb)
                     plot(xv(peakb(1)),yv(peakb(2)),'ow');
                 end
                 plot(xc,yc,'+k');
                 if np > 1
                     plot(minpt(1),minpt(2),'+w');
                 end
                 x = linspace(0,xr);
                 y =  sqrt(yr.^2 - (x.*yr/xr).^2);
                 x = [x fliplr(x) -x fliplr(-x)];
                 y = [y fliplr(-y) -y fliplr(y)];
                 plot(x+xc,y+yc,'r');
                 title(sprintf('%s%d:P%d %d peaks',DATA.Expts{eid}.Header.expname,eid,DATA.probe,np));
                 subplot(2,1,2);
                 hold off;
                 plot(Z(peaka(1),:));
                 hold on;
                 plot(Z(:,peaka(2)),'r');
                 if isfigure(DATA.figs.allprobes)
                     figure(DATA.figs.allprobes);
                     [nr,nc] = Nsubplots(length(DATA.probelist));
                     subplot(nr,nc,DATA.probe);
                 else
                     GetFigure(DATA.tag.clusterxy);
                 end
                 [DATA, ispk, dprime, details] = SetExptSpikes(DATA, eid, show,'useexpall');
             else
                 [DATA, ispk, dprime, details] = SetExptSpikes(DATA, eid, 0,'useexpall');

             end
             DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
             nspk = sum(details.nc);

             % IF there are two peaks, optimize cluster size
             if twopeaks & nspk > 0
                 oldC = C;
                 if minimize
                     DATA.ispk = ispk;
                     C = OptimizeDprime(DATA);
                     DATA.cluster{1,DATA.probe} = C;
                     [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe,0);
                 else
                     newd = dprime * 1.001;
                     dps(1) = dprime;
                     nl = 1;
                     frac = 0;
                     while newd > dprime
                         dprime = newd;
                         xr = xr.*1.1;
                         %                 yr = yr.*1.1;
                         C.x = [xc xr 1];
                         C.y = [yc yr 1];
                         DATA.cluster{1,DATA.probe} = C;
                         [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe, 0);
                         frac = d.nc(1)./length(ispk);
                         nl = nl+1;
                         dps(nl) = newd;
                     end
                     C.x = [xc xr./1.1 1];
                     dprime = dps(nl-1);
                     newd = dprime*1.001;
                     while newd > dprime
                         dprime = newd;
                         %                 xr = xr.*1.1;
                         yr = yr.*1.1;
                         %                C.x = [xc xr 1];
                         C.y = [yc yr 1];
                         DATA.cluster{1,DATA.probe} = C;
                         [DATA, newd, d] = SetSpkCodes(DATA,ispk,DATA.probe,0);
                         nl = nl+1;
                         dps(nl) = newd;
                     end
                     C.y = [yc yr./1.1 1];
                 end
                 frac = d.nc(1)./length(ispk);
                 if frac > 0.99  %didn't work
                     C = oldC;
                 end

                 DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
                 DATA.cluster{1,DATA.probe} = C;
                 [DATA, dprime, details] = SetSpkCodes(DATA,ispk,DATA.probe,show);
                 DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
             end
             if np == 1 | dprime < 1
                 C.dprime = dprime;
                 DATA = CheckAutoRate(DATA, eid, C, 10);
             end
             if show
                 DrawClusters(DATA,DATA.Expts{eid}.Cluster,0);
             end
             title(sprintf('%s%d %d/%d spikes: dprime %.1f',DATA.Expts{eid}.Header.expname,eid,details.nc(1),length(ispk),dprime));

             drawnow;
             return;





         elseif mode == 2
             [a,b] = max(z);
             py = smooth(y(b,1),5,'gauss');
             dens = smooth(a,5,'gauss');
             [a,maxi] = max(dens);
             maxpt = [b(maxi) maxi];
             for j = 1:length(dens)
                 dip(j) = dens(j)./max(dens(1:j));
             end
             id = find(diff(sign(diff(dip))) >0);
             GetFigure('AutoCut');
             subplot(2,1,1);
             hold off;
             plot(x(1,:),dip);
             hold on;
             for j = 1:length(id)
                 plot([x(1,id(j)+1) x(1,id(j)+1)],[0 1],':');
             end
             id = find(dip(2:end) > 0.5 & diff(dip) > 0.005);
             [a,b] = min(dip(id));
             if isempty(id)
                 id = find(dip < max(dip)/2);
                 ll = id(1);
                 ul = length(dens);
             else
                 ll = id(b);
                 %     ll = id(end);
                 ul= id(end)+10;
                 id = find(diff(dip(1:ll)) <= 0);


                 ll = id(end); %first positive ddip
                 peakdip = max(dip(ll:end))
                 id = find(diff(dip(ul:end)) > 0 & dip(ul+1:end) > peakdip);
                 if id
                     ul = ul+id(1);
                 else
                     ul = length(dens);
                 end
             end

             if x(1,ul) > max(Cx(ispk)) | ul == size(x,1)
                 xu = max(Cx(ispk)).*1.05;
             else
                 xu = x(1,ul);
             end
             plot([x(1,ll) x(1,ll)],[0 1],'r:');
             plot([xu xu],[0 1],'r:');
             subplot(2,1,2);
             hold off;
             pcolor(x,y,z);
             hold on;
             shading('flat');
             plot(x(1,:),py);

             xc = mean([xu x(1,ll)]);
             xr = (xu-x(1,ll))/2;
             yd = sum(z(:,ll:ul),2);
             yc = sum(y(:,1).*yd)./sum(yd);
             sigma = sqrt(sum(yd.*(y(:,1)-yc).^2)./sum(yd));
             yr = 6.* sigma;
             id = find(Cx > x(1,ll) & Cx <= xu);
             sigma = std(Cy(id));

             C.x = [xc xr 1];
             C.y = [yc yr 1];
             C.angle = 0;
             C.firstspk = NaN;
             C.lastspk = NaN;
             C.autocut = 2;
             C.params(1) = DATA.plot.clusterX;
             C.params(2) = DATA.plot.clusterY;
             C.Arange = DATA.clusterArange;
             C.Brange = DATA.clusterBrange;
             C.Erange = DATA.clusterErange;

             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA.Expts{eid}.gui.classified = 3;
             DATA.Expts{eid}.gui.clustertype = 3;

             GetFigure(DATA.tag.clusterxy);
             [DATA, ispk, dprime]  = SetExptSpikes(DATA, eid, show,'useexpall');
             if np == 1 | dprime < 2
                 DATA = CheckAutoRate(DATA, eid, C, 10);
             end
             DrawClusters(DATA,DATA.Expts{eid}.Cluster,0);
             return;
         end
         xid = find(DATA.Spikes.cx(ispk) > median(DATA.Spikes.cx(ispk)));
         yid = find(DATA.Spikes.cy(ispk) < median(DATA.Spikes.cy(ispk)));
         xc = prctile(DATA.Spikes.cx(ispk(yid)),90);
         xmax = prctile(DATA.Spikes.cx(ispk(yid)),99);
         yc = prctile(DATA.Spikes.cy(ispk(xid)),90);
         ymax = prctile(DATA.Spikes.cy(ispk(xid)),99);
         xr = (xmax-xc) * 2;
         if xr == 0
             xr = xc/3;
         end
         yr = (ymax-yc) * 1.5;
         if yr == 0
             yr = yc/3;
         end
         xc = xc.*1.4;
         yc = yc.*0.7;
         C.x = [xc xr 1];
         C.y = [yc yr 1];
         C.angle = 0;
         C.firstspk = NaN;
         C.lastspk = NaN;
         C.autocut = 1;
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
         if show & isfigure(DATA.xyfig);
             GetFigure(DATA.xyfig);
             hold off;
         end
         if ~isfield(DATA,'ptsize')
             DATA.ptsize = 4;
         end
         DATA = SetExptSpikes(DATA, eid, show,'useexpall');
         dprime = DATA.cluster{1,DATA.probe}.dprime;
         [rate, frac] = SumExptSpikes(DATA, eid, 1);
         if isinf(xr) | isinf(yr)
             fprintf('xc,yc Inf problem %f %f',xc,yc);
             return;
         end
         while rate < 10 & frac < 0.9
             xr = xr .*1.2;
             yr = yr.*1.2;
             if ~(xr < 1000) | ~(yr < 1000) %catches inf also..
                 fprintf('Big raduis',xc,yc);
             end
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             fprintf('Only %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xr,yr);
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA = SetExptSpikes(DATA, eid, show,'useexpall');
             dprime = DATA.cluster{1,DATA.probe}.dprime;
             [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
         while rate > 100 & frac > 0.1
             xc = xc .* 1.2;
             yc = yc .* 1.2;
             C.x = [xc xr 1];
             C.y = [yc yr 1];
             fprintf('Rate %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xc,yc);
             DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
             DATA = SetExptSpikes(DATA, eid, show,'useexpall');
             dprime = DATA.cluster{1,DATA.probe}.dprime;
             [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
         DATA.cluster{1,DATA.probe}.autocut = 1;
         DATA.cluster{1,DATA.probe}.dprime = dprime;
         DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
         fprintf('Probe %d Expt %d Cluster at %.3f+- %.3f, %.3f +- %.3f %.2fHz\n',DATA.probe,expno, C.x(1), C.x(2), C.y(1), C.y(2),rate);
         if show
             title(sprintf('probe %d, Expt %d %.f Hz Dp %.2f',DATA.probe,expno,rate,dprime));
         end
     else
         fprintf('Probe %d Expt %d No Spikes\n',DATA.probe,expno);
         C.firstspk = NaN;
         C.lastspk = NaN;
         C.autocut = 1;
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
     end

function OptimizeDprimeHit(a,b)
DATA = GetDataFromFig(a);
if isfield(DATA.plot,'useprobe') & sum(DATA.plot.useprobe)
    DATA = GetAllProbeFig(DATA);
    [nr, nc] = Nsubplots(length(DATA.probelist));
    probes = find(DATA.plot.useprobe);
    for j = probes;
        subplot(nr,nc,j);
        DATA.probe = DATA.probelist(j);
        DATA.ispk = DATA.AllClusters(j).spklist;
        C = OptimizeDprime(DATA);
        DATA.cluster{DATA.currentcluster,DATA.probe} = C;
        [DATA, newd, nc] = SetSpkCodes(DATA,DATA.ispk,DATA.probe,2);
        DrawClusters(DATA,DATA.cluster, 0);
    end
else
       
DATA.ispk = DATA.Expts{DATA.currentexpt}.gui.spks;
C = OptimizeDprime(DATA);
DATA.cluster{DATA.currentcluster,DATA.probe} = C;
[DATA, newd, nc] = SetSpkCodes(DATA,DATA.ispk,DATA.probe,2);
DrawClusters(DATA,DATA.cluster, 0);
end
set(DATA.toplevel,'UserData',DATA);

function C = OptimizeDprime(DATA)
    global ncall;
    global ycs;
    ncall = 0;
    ycs = [];
    eid = DATA.currentexpt;
    cl = DATA.currentcluster;
    C = DATA.cluster{cl,DATA.probe};
  
    byforce = 1;
%    options = optimset();
      if byforce
          x = [C.x(1) C.x(2) C.y(1) C.y(2)];
          for k= 1:3
     scales = [0.9:0.2:1.1];
     n = 0;
     for a = scales
         for b = scales
             for c = scales
                 for d = scales;
                     n = n+1;
                     xs(n,:) = x .* [a b c d];
                 end
             end
         end
     end
     for j = 1:n
    DATA.cluster{1,DATA.probe}.x = [xs(j,1) xs(j,2) C.x(3)];
    DATA.cluster{1,DATA.probe}.y = [xs(j,3) xs(j,4) C.x(3)];
    if isfield(DATA,'ispk')
    [DATA, dprimes(j), details] = SetSpkCodes(DATA,DATA.ispk,DATA.probe, 0);
    else
    [DATA, dprimes(j), details] = SetSpkCodes(DATA,DATA.spklist,DATA.probe, 0);
    end
    [dprime, best] = max(dprimes);
    x = xs(best,:);
     end
          end
      else 
          
%for some reason this seems to get suck in very local minima. I
%gueess because dprime is not a smoth function of the variables on
%a fine scale. see lemM0600 probe 4 expt 39
    options = optimset('TolFun',0.1,'TolX',0.01,'Largescale','off','simplex','on');
    [x, dp, flag, output] = fminsearch(@ClusterDprime, [C.x(1) C.x(2) C.y(1) C.y(2)],options, DATA, C.nspk);
%    [x, dp, flag] = fminsearch(@ClusterDprime, x,options, DATA, C.nspk);
    dprime = -dp;
      end
    C.x = [x(1) x(2) 1];
    C.y = [x(3) x(4) 1];
    C.dprime = dprime;
             
         
             
function dprime = ClusterDprime(x, DATA, nspks)
    global dprimes;
    global ncall;
    global ycs;

    eid = DATA.currentexpt;
    C = DATA.cluster{DATA.currentcluster,DATA.probe};
    if length(x) == 2
    C.x(2) = x(1);
    C.y(2) = x(2);

    else
    C.x(1) = x(1);
    C.x(2) = x(2);
    C.y(1) = x(3);
    C.y(2) = x(4);
    end
    ispk = DATA.ispk;

    DATA.cluster{1,DATA.probe} = C;
    [DATA, dprime, details] = SetSpkCodes(DATA,ispk,DATA.probe, 0);

    ns = details.nc(1);
    ncall = ncall+1;
    dprimes(ncall) = dprime;
    ycs(ncall,:) = x;
    if ns < nspks * 0.5 | ns > nspks * 1.5 | C.x(1) > details.maxx | C.y(1) > details.maxy 
        dprime = 1000;
    else
    dprime = -dprime;
    end
    
    
function DATA = CheckAutoRate(DATA, eid, C, minrate)
     xr = C.x(2);
     yr = C.y(2);
     xc = C.x(1);
     yc = C.y(1);
     show = 2;
     [rate, frac] = SumExptSpikes(DATA, eid, 1);
     ispk = DATA.Expts{eid}.gui.spks;
     xm = max(DATA.Spikes.cx(ispk));
     while rate < minrate & frac < 0.9
         xr = xr .*1.2;
         yr = yr.*1.2;
         if ~(xr < 1000) | ~(yr < 1000) %catches inf also..
             fprintf('Big raduis',xc,yc);
             rate = NaN;
         end
         C.x = [xc xr 1];
         C.y = [yc yr 1];
         fprintf('Only %.2f Hz (%.2f) Trying %.2f,%.2f\n',rate,frac,xr,yr);
         DATA.Expts{eid}.Cluster{1,DATA.probe} = C;
         [DATA, ispk, dprime]  = SetExptSpikes(DATA, eid, show,'useexpall');
         DATA.Expts{eid}.Cluster{1,DATA.probe}.dprime = dprime;
         DATA.Expts{eid}.Cluster{1,DATA.probe}.nspk = DATA.cluster{1,DATA.probe}.nspk;
         if ~isnan(rate)
         [rate, frac] = SumExptSpikes(DATA, eid, 1);
         end
     end
         
function [Expt, DATA] = CombinePlot(DATA, dlgs, varargin)
    
    spikelist = DATA.spikelist;
    ids = get(DATA.elst,'value');
    j = 1;
    while  j <= length(varargin)
        if strncmpi(varargin{j},'cluster',5)
            j = j+1;
            spikelist = varargin{j};
        elseif strncmpi(varargin{j},'ids',3)
            j = j+1;
            ids = varargin{j};
        end
        j = j+1;
    end
    
    needfields = {};
    if DATA.state.online < 2
        [DATA, Stimvals, needfields, ok] = CheckCombine(DATA, dlgs);
        if ok < 0
            Expt = [];
            return;
        end
    end

    
%
% when combining expts with certain differences, make sure these are
% recorded for each trial, so that they can be treared differently
% add to this list any parameters shown manually with DATA.show
    CheckF = {'sx' 'Fr' 'rb' 'ap' 'jx' 'ns' 'sf' 'or' 'backxo' 'bo' 'backyo' 'co' 'Bc' 'bh' 'SpikeGain' 'dd'};
    f = fields(DATA.show);
    for j = 1:length(f)
      if isempty(strmatch(f{j},{'ed'},'exact')) && DATA.show.(f{j})
            CheckF = {CheckF{:} f{j}};
        end
    end
    for f = 1:length(CheckF)
        vals = [];
        for j = ids(1:end)
            if isfield(DATA.Expts{DATA.expid(j)}.Stimvals,CheckF{f})
                vals(j) = DATA.Expts{DATA.expid(j)}.Stimvals.(CheckF{f});
            else 
                vals(j) = Inf;
            end
        end
        if length(unique(vals(ids))) > 1
            for j = ids(1:end)
            DATA.Expts{DATA.expid(j)}= FillTrials(DATA.Expts{DATA.expid(j)},CheckF{f});
            end
        end
    end
    tfields = {};
    Frs = [];
    Comments = [];
    for j = ids(1:end)
        tfields = union(tfields,fields(DATA.Expts{DATA.expid(j)}.Trials));
        if isfield(DATA.Expts{DATA.expid(j)}.Stimvals,'Fr')
        Frs(j) = DATA.Expts{DATA.expid(j)}.Stimvals.Fr;
        else
            Frs(j) = 0;
        end
        Comments = [Comments  DATA.Expts{DATA.expid(j)}.Comments];
    end
    tfields = union(tfields,needfields);
    if length(unique(Frs(ids))) > 1
        for j = ids(1:end)
            DATA.Expts{DATA.expid(j)}= FillTrials(DATA.Expts{DATA.expid(j)},'Fr');
        end
    end
    j = DATA.expid(ids(end));
    DATA.combineids = DATA.expid(ids);
    if length(DATA.extype) == 1
    DATA.allcombineids{DATA.extype} = DATA.combineids;
    end
    Expt.Header = DATA.Expts{j}.Header;
    Trials = [];
    Pulses = [];
    nb = 1;
    Clusters = {};
    if DATA.plot.autoclustermode < 3
        DATA.autocut = 1;
    else
        DATA.autocut = 0;
    end
    counts = {}; nspk = 0;
    for j = ids(1:end)
        eid = DATA.expid(j);
        BlockStart(nb) = DATA.Expts{DATA.expid(j)}.Trials(1).Trial;
        exped(nb) = GetEval(DATA.Expts{DATA.expid(j)},'ed');
        Header = DATA.Expts{DATA.expid(j)}.Header;
        if isfield(DATA.Expts{DATA.expid(j)},'Cluster')
                Clusters{nb} = DATA.Expts{DATA.expid(j)}.Cluster;
                recount = 0;
        end
        if DATA.state.nospikes
            recount = 1; %% for chagning probes
        elseif DATA.autocut & (size(Clusters) < nb | size(Clusters{nb},2) < DATA.probe ...
                | isempty(Clusters{nb}) | isempty(Clusters{nb}{1,DATA.probe})...
                | ~isfield(Clusters{nb}{1,DATA.probe},'x')...
                | (isfield(Clusters{nb}{1,DATA.probe},'autocut') & Clusters{nb}{1,DATA.probe}.autocut > 0  &DATA.state.redoautocut))
            DATA = AutoCut(DATA,DATA.expid(j),j);
            Clusters{nb} = DATA.Expts{DATA.expid(j)}.Cluster;
            recount = 1;
        else
            recount = 0;
        end
        if DATA.state.psychonly
            recount = 0;
        elseif DATA.state.online == 0 | recount
%            if iscluster(DATA.Expts{eid}.Cluster,DATA.currentcluster,DATA.probe)
%                DATA.cluster = DATA.Expts{eid}.Cluster;
%            end
            [DATA, counts{j}] = CountSpikes(DATA,DATA.expid(j)); %% recount in case cluster # changed
        end
%tfields is the combined fields of all expts (above). Sao any differences
%mean it is missing from the current Expt;
        newf = setdiff(tfields, fields(DATA.Expts{DATA.expid(j)}.Trials));
        for k = 1:length(newf)
            if strcmp(newf{k},'FalseStart')
                [DATA.Expts{DATA.expid(j)}.Trials.(newf{k})] = deal(0);
            else
                DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},newf{k});
            end
        end
        newtrials = DATA.Expts{DATA.expid(j)}.Trials;
        Trials = [Trials newtrials];
        if isfield(DATA.Expts{DATA.expid(j)},'ExcludeId')
        end
        if isfield(DATA.Expts{DATA.expid(j)},'Pulses')
            Pulses = [Pulses DATA.Expts{DATA.expid(j)}.Pulses];
        end
        if Header.trange(1) < Expt.Header.trange(1)
            Expt.Header.trange(1) = Header.trange(1);
        end
        if Header.trange(2) > Expt.Header.trange(2)
            Expt.Header.trange(2) = Header.trange(2);
        end
        nb = nb+1;
    end
    for j = 1:length(counts)
        for k = 1:length(counts{j})
            nspk(j,k) = counts{j}(k);
        end
    end
    eid = DATA.expid(ids(end));
    Expt.Stimvals = Stimvals;
    Expt.Trials = Trials;
    Expt.Header.Name = BuildName(DATA.Expts{eid}.Header.Name);
    Expt.Header.BlockStart = BlockStart;
    Expt.Header.depths = exped;
    Expt.Header.Clusters = Clusters;
    Expt.Header.Combined = ids;
    Expt.Header.Combineids = DATA.expid(ids);
    Expt.Header.Spikelist = DATA.spikelist;
    Expt.Header.CombineDate = now;
    Expt.Header.nspk = sum(nspk);
    if isfield(DATA,'AllClusters')
        Expt.Header.probe = DATA.probe;
        Expt.Header.dips = DATA.AllClusters(DATA.probe).dips;
    end

    if ~isempty(Comments)
    Expt.Comments = Comments;
    end
    if isfield(DATA,'Comments') & isfield(DATA.Comments,'Peninfo')
        Expt.Header.Peninfo = DATA.Comments.Peninfo;
                    id = strfind(DATA.Comments.Peninfo.trode,'Contact');
            if length(id)
                x = id(1);
                id = strfind(DATA.Comments.Peninfo.trode(id:end),' ');
                x = sscanf(DATA.Comments.Peninfo.trode(id+x:end),'%d');
                Expt.Header.probesep = x;
            end
    end

    if ~DATA.state.online & ~isfield(DATA,'AllClusters') & ~DATA.state.psychonly
    Expt.Header.SpkStats = GetSpkStats(DATA);
    end
    if isfield(Expt.Trials,'FalseStart')
        id = find([Expt.Trials.FalseStart] > 0);
        if ~isempty(id)
            for j = 1:length(id)
                Expt.Trials(id(j)).Start = Expt.Trials(id(j)).Start - Expt.Trials(id(j)).FalseStart;
            end
            msgbox(sprintf('%s,%s: %d Trials with long delays',Expt.Header.Name,Expt.Header.expname,length(id)));
        end
    end
    if ~isempty(Pulses)
        Expt.Pulses = Pulses;
    end
    if ~isfield(Expt.Header,'Spike2Version')
        Expt.Header.Spike2Version = 1.0;
    end
    
    if DATA.state.uselfp & length(ids) > 0  %reload LFP to match lengths etc
        Expt = LoadSpike2LFP(Expt,'reload');
    end
    id = strmatch(DATA.Expts{eid}.Header.expname,DATA.expstrs,'exact');
    if Stimvals.BlockedStim
        DATA.outname = Expt2FileName(DATA, Expt, DATA.spikelist(1));
        set(DATA.saveitem,'string',DATA.outname);
    end
    [Expt, plotres] = PlotCombined(DATA, Expt);
    if DATA.state.autofit
        fit = FitExpt(plotres,'plotfit');
        DATA = AddFitToData(DATA, plotres,fit);
    end
  
function [Expt, res] = PlotCombined(DATA, Expt, varargin)
res = [];
xargs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'hold',4)
        xargs = {xargs{:} varargin{j}};
    end
    j = j+1;
end
 if isfield(Expt,'ExcludeCluster')
    for j = 1:length(Expt.ExcludeCluster{1})
        id = find([Expt.Trials.Trial] == Expt.ExcludeCluster{1}(j));
        if id
            Expt.Trials(id).Trial = -Expt.Trials(id).Trial;
        end
    end
 end
 if ismember(DATA.plot.showcp,[5 6]) & isfield(Expt.Trials,'RespDir')
     psfargs = {};
     on = get(findobj('Tag','ShowN','Parent',DATA.toplevel),'value');
     if on
         psfargs = {psfargs{:} 'shown'};
     end
     GetFigure(DATA.tag.psych);
     hold off;
     if DATA.plot.showcp == 6
         ExptPsych(Expt, 'smooth',10,psfargs{:});
     else
         ExptPsych(Expt,psfargs{:});
     end
     return;
 end
    GetFigure(DATA.tag.dataplot);
    args = PlotArgs(DATA,Expt,'combined');
    args = {xargs{:} args{:}};
    res = PlotExpt(Expt,args{:});
    if DATA.state.plotseq == 4 & isfield(res,'fig')
        set(res(1).fig, 'WindowButtonDownFcn',@FigButtonPressed);
        set(res(1).fig, 'WindowButtonUpFcn',@FigButtonReleased);
        dat = get(res(1).fig,'UserData');
        dat.parentfigtag = DATA.tag.top;
        set(res(1).fig,'UserData',dat);
    elseif DATA.state.plotseq == 10  %show acloop plot
        PlotRC(res,'acloop');
    elseif DATA.state.plotseq == 11  %show acloop plot
        PlotRC(res,'acresp');
    end
    t = get(get(gca,'title'),'String');
    spikelist = WhichClusters(DATA.toplevel);
    edepth = GetEval(Expt,'ed');
%    title([t 'Cl' sprintf(' %d',WhichClusters(DATA.toplevel))]);
    title([t sprintf(' P%d',DATA.probe) sprintf(' Cl %d',spikelist) sprintf('ed %.2f',edepth)]);
    if DATA.plot.showem
        GetFigure(DATA.tag.emplot);
        Expt = LoadEmData(Expt);
        PlotExptEM(Expt);
    end
    if DATA.state.uselfp
        GetFigure('LFP');
        hold off;        
        CalcLFPPulse(Expt,DATA.AllData,'plot');
        if size(Expt.Trials(1).LFP,2) > 5
            if DATA.plot.lfpplot == 1
                PlotMLFP(Expt,'stack',0);
            elseif DATA.plot.lfpplot ==2
                PlotMLFP(Expt,'image');
            end
        end
    end
        
%    save(DATA.outname,'Expt'); %use save to do the saving....

    
   
function DATA = ReadDir(DATA, name, varargin)  %% Online Plots

d = dir(name);
[a,b] =sort([d.datenum]);
d = d(b);
reindex = 0;
args = {};
j = 1;
nskip = 0;
DATA.cellfile = [name, '/Cells.mat'];

while j <= nargin-2
        if strncmpi(varargin{j},'relist',3)
            reindex =1;
            args = {args{:} 'relist'};
        elseif strncmpi(varargin{j},'skip',3)
            j = j+1;
            nskip = varargin{j};
        else
            args = {args{:} varargin{j}};
        end
    j = j+1;
end
    expnames = {};
    if DATA.logfid == 0
        DATA.logfid = fopen([name,'/log.txt'],'a');
    end

    if isempty(DATA.Expts)
        nexp = 1;
        SpkId = [];
        Spikes = [];
        Trialids = [];
        TrialTimes = [];
    else
        for j=1:length(DATA.Expts)
            expnames{j} = splitpath(DATA.Expts{j}.Header.Name);
 %           DATA.Expts{j}.gui.classified = 0;
 % don't reset this here - resets clusters for files already read. Reset
 % it below only when a file is read.
        end
        nexp = j+1;
        SpkId = DATA.AllData.SpikeIdx;
        Spikes = DATA.AllData.Spikes;
        Trialids = DATA.AllData.Trialids;
        TrialTimes = DATA.AllData.TrialTimes;
        All = DATA.AllData;
        Expts = DATA.Expts;
    end
    DATA.defaults.starttrial = 1;
    lastn = nexp;
    nbad = length(DATA.badnames)+1;
    badidx = 1:nbad-1;
    nex = 0;
    nexpts = 0;
    for j = 1:length(d)
        if regexp(d(j).name,'Expt[0-9]*.mat') & ...
                d(j).bytes > 128 & .....
                isempty(strfind(d(j).name,'idx.mat')) & ...    %exclude the .idx file
                d(j).datenum < now-(1/(24 * 60 * 60)) & ... %at least  1 sec old
                isempty(strmatch(d(j).name,{expnames{:}}))%dont read if we already have
            nexpts = nexpts+1;
        end
    end
    fprintf('%d Epxts in %s\n',nexpts,name);
    for j = 1:length(d)
        if regexp(d(j).name,'Expt[0-9]*.mat') & ...
                d(j).bytes > 128 & .....
                isempty(strfind(d(j).name,'idx.mat')) & ...    %exclude the .idx file
                d(j).datenum < now-(1/(24 * 60 * 60)) & ... %at least  1 sec old
                isempty(strmatch(d(j).name,{expnames{:}}))%dont read if we already have
            nex = nex+1;
            if nex <= nskip & nex < nexpts
                id = NaN;
            elseif length(DATA.badnames)
                id = strmatch((d(j).name),DATA.badnames(badidx));
            else
                id = [];
            end
            if isempty(id) || (~isnan(id) && d(j).datenum > DATA.badtimes(id(end)))
            [trls, exps, All] = APlaySpkFile([name '/' d(j).name],'Defaults',DATA.defaults,'online',args{:});
            if isfield(trls,'DataType')
                DATA.filetype = trls.DataType;
            end
                
            if isempty(exps)
                fprintf('%s No expts\n',d(j).name);
                DATA.badnames{nbad} = d(j).name;
                DATA.badtimes(nbad) = d(j).datenum;
                nbad = nbad+1;
                badidx = 1:nbad-1;
            elseif iscell(exps)
                Expts{nexp} = exps{1};
                Expts{nexp}.gui.classified = 0;
                Expts{nexp}.gui.counted = 0;
                Expts{nexp}.gui.clustertype = 0;
                Expts{nexp}.gui.firsttrial = 1+length(Trialids);
                if isfield(trls,'bstimes')
                id = find(trls.bstimes > Expts{nexp}.Header.Start & trls.bstimes < Expts{nexp}.Header.End);
                Expts{nexp}.bstimes = trls.bstimes(id);
                end
                SpkId = [SpkId; trls.Spkid];
                newt = [Expts{nexp}.Trials.Trial];
                Trialids = [Trialids newt];
                for k = 1:length(Expts{nexp}.Trials)
                    news(k) = Expts{nexp}.Trials(k).Start(1);
                end
                TrialTimes = [TrialTimes news];
                Expts{nexp}.gui.ntrials = length(newt);
                if nexp > 1 && isfield(Spikes,'values') && size(Spikes.values,2) == size(All.Spikes.values,2)
                    Spikes.values = cat(1,Spikes.values, All.Spikes.values);
                    if isfield(All.Spikes,'dVdt')
                    Spikes.dVdt = cat(1,Spikes.dVdt, All.Spikes.dVdt);
                    end
                    Spikes.codes = [Spikes.codes; All.Spikes.codes];
                    Spikes.times = [Spikes.times; All.Spikes.times];
                else
                    Spikes = All.Spikes;
                end

                nexp = nexp+1;
                DATA.defaults.starttrial = 1+ exps{1}.Trials(end).Trial;
                
                         if isfield(trls,'Probes');
                             for j = 1:length(trls.Probes)
                                 probenames{trls.Probes(j).probe} = trls.Probes(j).var;
                                 Probes(trls.Probes(j).probe) = trls.Probes(j);
                             end
                         end
            else
                msgbox(sprintf('%s Incomplete Expt',d(j).name));         
            end
        end
        end
    end

    if exist('probenames','var') %can be empty for online relist
    np = 0;
    probes = [];
    for j = 1:length(probenames)
        if ~isempty(probenames{j})
            np = np+1;
            DATA.probes(np) = Probes(j);
            probes = [probes j];
        end
    end
                DATA.probelist = probes;
            DATA.probevars = {probenames{probes}};
            DATA.probenames = cellstr(int2str(probes'));
            if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
                DATA.probe = DATA.probelist(1);
            end
%            DATA.probes = Probes;
            if length(DATA.probelist) > 1 
                DATA = AddMultiProbeGUI(DATA);
                DATA.probesource = FindProbeSources(DATA);
            else
                DATA.probesource = 1;
            end
    elseif strncmp(DATA.filetype,'Grid',4)
        np = sscanf(DATA.filetype,'GridData %d');
        if isempty(np)
            np = 96; % default
        end
        DATA.probelist = 1:np;
        DATA.probesource = 100;
        DATA.probenames = cellstr(int2str(DATA.probelist'));
        DATA.state.applylastcluster = 1;
        DATA.state.classifyallprobes = 1;

        for j = 1:np
            DATA.probes(np).probe = j;
        end
            if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
                DATA.probe = DATA.probelist(1);
            end

        if np > 1
                DATA = AddMultiProbeGUI(DATA);
        end
    end
    it = findobj(DATA.toplevel,'Tag','ProbeId');
    set(it,'string',DATA.probenames,'value',DATA.probe);


    if nexp == 1
        questdlg(sprintf('No expts in %s',name),'test','OK','OK');
    else
        DATA = ListExpts(DATA,Expts);
        DATA.Expts = Expts;
        DATA.AllData = All;
        DATA.AllData.Spikes = Spikes;
        DATA.AllData.SpikeIdx = SpkId;
        DATA.AllData.Trialids = Trialids;
        DATA.AllData.TrialTimes = TrialTimes;
        DATA.AllData.pcs = [];
        if DATA.state.autosetlist && DATA.state.online == 1
            for j = lastn:nexp-1
              DATA = SetExptSpikes(DATA, j, 0);
            end
        end
    end
    DATA.name = name;
    DATA.state.online = 1;
    

    
function PlotSpike(DATA, ispk, probe)

    set(0,'CurrentFigure',DATA.svfig);
    if nargin == 2
    j = DATA.AllData.Spikes.codes(ispk,2)+1;
    if DATA.plot.dvdt == 2
          set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',DATA.AllData.Spikes.values(ispk,2:end));
    elseif DATA.plot.dvdt
          set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.dVdt,2)]);
    else
        set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.values(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.values(ispk,:),2)]);
    end
    title(sprintf('%d: Cl %d at %.3f',ispk,j-1,DATA.AllData.Spikes.times(ispk)/10000)); 
    else
        j = DATA.AllSpikes{DATA.probe}.codes(ispk,2)+1;
        set(DATA.svh(j), 'Ydata', DATA.AllSpikes{probe(1)}.values(ispk,:)-DATA.plot.SpikeMaxV/2,'Xdata',[1:size(DATA.AllSpikes{probe(1)}.values,2)]);
        id = find(ispk == DATA.sids{1});
        if length(id)
            jspk = DATA.sids{2}(id);
        set(DATA.svh(j+DATA.nclusters+1), 'Ydata', DATA.AllSpikes{probe(2)}.values(jspk,:)+DATA.plot.SpikeMaxV/2,'Xdata',[1:size(DATA.AllSpikes{probe(2)}.values,2)]);
        end
        title(sprintf('%d: Cl %d at %.3f',ispk,j-1,DATA.AllSpikes{DATA.probe}.times(ispk)/10000));
    end
    drawnow;

function PlayNextSpike(a,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(a);
if DATA.ISIpair
    PlotISIPair(DATA,DATA.ISIpair);
    DATA.ISIpair = DATA.ISIpair+1;
    set(DATA.toplevel,'UserData',DATA);
    return;
end
set(0,'CurrentFigure',DATA.svfig);
if isfield(DATA,'AllSpikes')
    PlotSpike(DATA,DATA.currentspike,[DATA.probe DATA.xprobes]);
else
PlotSpike(DATA,DATA.currentspike);
if isfigure(DATA.xyfig)
    set(0,'CurrentFigure',DATA.xyfig);
    PlotSpikeXY(DATA,DATA.currentspike,DATA.spkcolor{DATA.AllData.Spikes.codes(DATA.currentspike,2)+1});
end
end
DATA.currentspike = DATA.currentspike+1;

if length(DATA.probes) > 1 && DATA.plot.timebyspikeprobe > 0
    p = DATA.probe;
    probes = find(DATA.plot.useprobe);
    set(0,'CurrentFigure',DATA.timefig);
    hold off;
    tw = 100;
    set(gca,'xlim',[0 tw*2]);
    nsmp = size(DATA.AllSpikes{p}.values,2);
    hscale = DATA.AllSpikes{p}.interval * 10000;
    ts = DATA.AllSpikes{p}.times(DATA.currentspike);
    for j = 1:length(probes)
        [sids{j}, alltimes{j} allcodes{j}] = FindSpikes(DATA, [ts-tw ts+tw], probes(j),[]);
        voff(j) = DATA.plot.SpikeVsep*(j-1);
        x = [];
        y = [];
        id = find(alltimes{j} > ts-tw & alltimes{j} < ts+tw);
        for t = 1:length(id)
            x = (alltimes{j}(id(t))-ts+tw) +[1:nsmp]*hscale;
            y = DATA.AllSpikes{probes(j)}.values(sids{j}(id(t)),:)+voff(j);
            plot(x,y,'color',DATA.spkcolor{allcodes{j}(id(t))+1});
            hold on;
        end
    end
    set(gca,'xlim',[0 tw*2],'ylim',[-DATA.plot.SpikeMaxV max(voff)+DATA.plot.SpikeMaxV]);
    drawnow;
end


set(DATA.toplevel,'UserData',DATA);

function PlotISIPair(DATA, pair)

GetFigure(DATA.tag.spikev);
hold off;
spk = find(DATA.AllData.Spikes.times >= floor(DATA.isis(pair)));
spk = spk(1);
tsmp = [1:size(DATA.AllData.Spikes.values,2)] * DATA.AllData.Spikes.interval * 10000;
plot(tsmp,DATA.AllData.Spikes.values(spk-1,:));
[a,ta] = max(DATA.AllData.Spikes.values(spk-1,:));
[b,tb] = max(DATA.AllData.Spikes.values(spk,:));
    hold on;
isi = diff(DATA.AllData.Spikes.times([spk-1 spk]));
times = isi+tsmp+20000*DATA.AllData.Spikes.interval;
isi = isi + (tb-ta) * DATA.AllData.Spikes.interval * 10000;
plot(times,DATA.AllData.Spikes.values(spk,:));
title(sprintf('Time %.0f ISI %.1f',DATA.AllData.Spikes.times(spk-1),isi));

function CutTrial(a,b)
%DATA = combine('getstate');
DATA = GetDataFromFig(a);

t = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial;
if length(DATA.probelist) > 1 & isfield(DATA,'CellList')
    cells = CellsSelected(DATA);
     for j = 1:length(cells)
         DATA.CellQuality(cells(j),t) = 3;
     end
else
DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial = -abs(t);
end
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'string',sprintf('%d|',[DATA.Expts{DATA.currentexpt}.Trials.Trial]),'value',1);
set(DATA.toplevel,'UserData',DATA);
PlayOneTrial(DATA,a,1);
        
function PlayLastTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA, a,-1);

function PlayNextTrial(a, b)
DATA = GetDataFromFig(a);
PlayOneTrial(DATA,a,1);

function SelectTrial(a, b)
DATA = GetDataFromFig(a);
c = get(findobj(a,'Tag','ChooseTrial'),'value');
PlayOneTrial(DATA,c,0);

function PlayOneTrial(DATA, a, b, varargin)
setgui = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nosetgui',6)
        setgui = 0;
    end
    j =j+1;
end

%DATA = combine('getstate');
expid = DATA.currentexpt;

if b < 0  & DATA.currenttrial > 1 %step back
    DATA.currenttrial = DATA.currenttrial-1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 & DATA.currenttrial > 1 
        DATA.currenttrial = DATA.currenttrial-1;
    end
elseif b > 0 & DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) %step back
    DATA.currenttrial = DATA.currenttrial+1;
    while DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial < 0 &  ...
            DATA.currenttrial < length(DATA.Expts{DATA.currentexpt}.Trials) 
        DATA.currenttrial = DATA.currenttrial+1;
    end
elseif b == 0 % step to a trial on list
    DATA.currenttrial = find([DATA.Expts{DATA.currentexpt}.Trials.Trial] == a);    
    DATA.currenttrial = a;    
end
if DATA.currenttrial <= 0
    DATA.currenttrial = 1;
end
Trial = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial);
if Trial.Trial < 0 && b == 0 && setgui %%manually select -Trial = include again
    DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).Trial = abs(Trial.Trial);
    it = findobj(DATA.svfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
end
itrial = find(DATA.AllData.Trialids == abs(Trial.Trial));
set(0,'CurrentFigure',DATA.svfig);
hold off;
%DATA = PlotTrialSpikes(DATA, itrial, mycolors, DATA.clusters);
if DATA.state.recut
    nc = size(DATA.cluster,1)+1;
else
    nc = DATA.s2clusters;
end
Trial.ed = GetEval(DATA.Expts{DATA.currentexpt},'ed',DATA.currenttrial);
probes =  [DATA.probe DATA.xprobes];
if (length(probes) > 1 & DATA.plot.syncoverlay == 0) || DATA.plot.showwave
    Aargs = {'timemode'};
    if  ~isfield(DATA,'timefig')
        DATA.timefig = GetFigure('SpikeTime');
    end
    set(0,'CurrentFigure',DATA.timefig);
else
    Aargs = {};
end


if DATA.plot.voltxy
    set(gca,'xlim',[-5 5],'ylim',[-5 5]);
end
if DATA.syncsign < 2
    PlotTrialSyncSpikes(DATA, [Trial.Start(1) Trial.End(end)], [DATA.probe DATA.xprobes], mycolors,'Trial',Trial);
else
    for j = 2:length(probes)
        DATA = APlotTrialSpikes(DATA, [Trial.Start(1)-DATA.state.preperiod Trial.End(end)+DATA.state.postperiod], mycolors, nc, 1,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(DATA.svhn), Aargs{:});
    end
    DATA = APlotTrialSpikes(DATA, [Trial.Start(1)-DATA.state.preperiod Trial.End(end)+DATA.state.postperiod], mycolors, nc, 1,'Trial',Trial,'probe',[probes(1) 1 length(probes)],Aargs{:});
    set(0,'CurrentFigure',DATA.svfig);
    if isfield(DATA,'CellList') && size(DATA.CellList,2) >= Trial.Trial
    for j = 1:length(probes)
        id = find(DATA.CellList(:,Trial.Trial) == probes(j));
        if id
            text(10,DATA.plot.currentoffset(j),['cell' num2str(id(1))]);
        end
    end
    end
end
if length(probes) > 1 & isfield(DATA,'timefig') %also show spikes over time
    set(0,'CurrentFigure',DATA.timefig);
    vh = DATA.tvh;
    if DATA.plot.timebyspikeprobe > 0
        p = DATA.probe;
        tw = 100;
        set(gca,'xlim',[0 tw*2]);
        [sid, stimes, codes] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], p,[]);
        sid = sid(codes ==1);
        stimes = stimes(codes ==1);
        dt = [-110:10:110];
        dts = [];
        for j = 1:length(probes)
            [sids{j}, alltimes{j} allcodes{j}] = FindSpikes(DATA, [Trial.Start(1) Trial.End(end)], probes(j),[]);
            voff(j) = DATA.plot.SpikeVsep*(j-1);
            diffs{j} = [];
              for k = 1:length(stimes)
                  dts(k,:) = hist(stimes(k) - alltimes{j},dt);
              end
              corrs(j,:) = sum(dts);
        end
        nsmp = size(DATA.AllSpikes{probes(1)}.values,2);
        hscale = DATA.AllSpikes{probes(1)}.interval * 10000;
        hold off;
        for k = 1:length(sid)
            ts = stimes(k);
        for j = 1:length(probes)
            x = [];
            y = [];
            id = find(alltimes{j} > ts-tw & alltimes{j} < ts+tw);
            for t = 1:length(id)
                x = (alltimes{j}(id(t))-ts+tw) +[1:nsmp]*hscale;
                y = DATA.AllSpikes{probes(j)}.values(sids{j}(id(t)),:)+voff(j);
                plot(x,y,'color',DATA.spkcolor{allcodes{j}(id(t))+1});
                diffs{j} = [diffs{j} alltimes{j}(id(t))-ts];
                hold on;
            end
        end
%        hold off;
        set(gca,'xlim',[0 tw*2],'ylim',[-DATA.plot.SpikeMaxV max(voff)+DATA.plot.SpikeMaxV]);
        drawnow;
        pause(0.05);
        end
        DATA.currentspike = sid(1);
    else
    for j = 1:length(probes)
        DATA = APlotTrialSpikes(DATA, [Trial.Start(1) Trial.End(end)], mycolors, nc, 0,'Trial',Trial,'probe',[probes(j) j length(probes)],'lineoff',(j-1)*(1+DATA.svhn),'timemode');
    end
    end
end
if DATA.state.uselfp
   GetFigure('TrialLFP');
    PlotLFPRaw(DATA.state,Trial,DATA.Expts{expid}.Header.LFPsamplerate);
end
CheckSpoolButton(DATA);
set(DATA.toplevel,'UserData',DATA);
set(0,'CurrentFigure',DATA.svfig);
tid = DATA.Expts{DATA.currentexpt}.Trials(DATA.currenttrial).id;
if setgui
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'value',DATA.currenttrial);
end
    
function PlotLFPRaw(state, Trial, crrate);
  
    times = [1:size(Trial.LFP,1)] .* crrate;
    times = times - (Trial.Start(1)-Trial.lfptime)/10000;
    plot(times,Trial.LFP);
    
    
function DATA = PlotTrialSpikes(DATA, itrial, colors, clusters)
    spka = DATA.AllData.SpikeIdx(itrial,1);
    spkb = DATA.AllData.SpikeIdx(itrial,2);
    Spks = DATA.AllData.Spikes;
    if spka
        for spk = spka:spkb;
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            hold on;
            energy  = sum(diff(adc).^2);
            DATA.Spikes.energy(spk)= energy;
            svar(spk) = var(adc);
            DATA.Spikes.vw(spk) = svar(spk)/energy;
        end
        title(sprintf('Trial %d (id%d)',DATA.AllData.Trialids(itrial),...
            DATA.AllData.Trialids(itrial)));

        drawnow;


        set(0,'CurrentFigure',DATA.xyfig);
        for j = 1:length(clusters)
            sp = intersect(clusters{j},[spka:spkb]);
            plot(DATA.Spikes.energy(sp),DATA.Spikes.vw(sp),'.','color',colors{j});
        end
    end

function [x,DATA] = GetSpikeVals(DATA, ispk, values, dVdt, type, recalc, pcs)

%Do NOT change the order of these definiitons. Cluster params are saved
%as integers....
SPKENERGY=1;
SPKVARE=2;
SPKMAXRATE = 3;
SPKPREMINRATE = 4;
SPKPEAKTIME = 5;
SPKPREMIN = 6;
SPKPEAK = 7;
SPKVAR = 8;
SPKSYMMETRY = 9;
SPKCENTROID = 10;

SPKMINRATE = 11;
SPKMAXRATEA = 12;
SPKMINRATEA = 13;
SPKMEANRATEA = 14;
SPKMAXRATEB = 15;
SPKMINRATEB = 16;
SPKMEANRATEB = 17;
SPKMEANA = 18;
SPKMEANB = 19;
SPKMAXA = 20;
SPKMAXB = 21;
SPKMIN =22;
SPKMINA = 23;
SPKMINB = 24;
SPKHEIGHT = 25;
SPKENERGYA=26;
SPKENERGYB=27;
SPKPREMAXRATE = 28;
TEMPLATEA = 29;
TEMPLATEB = 30;
TEMPLATEC = 31;
SPKSYMNEG = 32;
SPKPCA1= 33;
SPKPCA2= 34;
SPKVARTESTA= 35;
SPKVARTESTB= 36;
SPKPCA3= 37;
SPKPCA4= 38;
SPKPCA12= 39;
SPKMEAN =40;
SPKMAXACCEL = 41;
SPKMINACCEL = 42;
SPKMEANACCEL = 43;
SPKVAREA = 44;
SPKVAREB = 45;
SPKRMSV = 46; %RMS velocity
SPKADC1=47;
SPKADC2=48;
SPKADC3=49;
SPKADC4=50;
ENERGY1 = 51;
ENERGY2 = 52;
ENERGY3 = 53;
ENERGY4 = 54;
if isnan(values) %% return names. Names must match order of actual variable (above). Order can be anything
    if DATA.subprobes > 1
        addsub = 1;
    else
        addsub = 0;
    end
    x = {'Energy' 'Var/Energy' 'MaxRate' 'PreMinRate' 'PeakTime' 'PreMin' 'Peak' 'Var' 'Symmetry' 'Centroid' 'Minrate' 'Maxrate(A)' ... 
    'Minrate(A)' 'Meanrate(A)' 'Maxrate(B)' 'Minrate(B)' 'Meanrate(B)' 'Mean(A)'...
    'Mean(B)' 'Max(A)' 'Max(B)' 'Min' 'Min(A)' 'Min(B)' 'Height' 'EnergyA' 'EnergyB' 'PreMaxRate' ...
    'Template 1' 'Template 2' 'Template 3' '-Symmetry' 'PCA1' 'PCA2' 'test1' 'test2' 'PCA3' 'PCA4' 'PCA1-PCA2' 'Mean' 'AccelMax' 'AccelMin' 'AccelMean'...
    'Var/sqrt(energy)' 'sqrt(Var/Energy)' 'sqrt(Energy)' 'ADC1' 'ADC2' 'ADC3' 'ADC4'};
   DATA = [SPKENERGY SPKVARE SPKMAXRATE SPKVAR SPKPEAK SPKMIN SPKPREMINRATE SPKPEAKTIME SPKPREMIN ...
       SPKSYMMETRY SPKCENTROID SPKMINRATE ...
       SPKMAXRATEA SPKMINRATEA SPKMEANRATEA SPKMAXRATEB SPKMINRATEB ...
    SPKMEANRATEB SPKMEAN SPKMEANA SPKMEANB SPKMAXA SPKMAXB SPKMIN SPKMINA SPKMINB ...
    SPKHEIGHT SPKENERGYA SPKENERGYB SPKPREMAXRATE TEMPLATEA TEMPLATEB ...
    TEMPLATEC SPKSYMNEG SPKPCA1 SPKPCA2 SPKVARTESTA SPKVARTESTB SPKPCA3 SPKPCA4 SPKPCA12  SPKMEAN SPKMAXACCEL SPKMINACCEL SPKMEANACCEL ...
    SPKVAREA SPKVAREB SPKRMSV SPKADC1 SPKADC2 SPKADC3];
if addsub
DATA = [DATA ENERGY1 ENERGY2  ENERGY3 ENERGY4 ];
x = {x{:} 'Energy 1' 'Energy 2' 'Energy 3' 'Energy 4'};
end
return;
end


p = DATA.probe;
if p > size(DATA.cluster,2)
    DATA.cluster{DATA.currentcluster,p} = {};
end

if DATA.currentcluster > size(DATA.cluster,1) || (~isfield(DATA.cluster{DATA.currentcluster,p},'Arange') & isfield(DATA.cluster{1,p},'Arange'))
    DATA.cluster{DATA.currentcluster,p}.Arange = DATA.cluster{1,p}.Arange;
    DATA.cluster{DATA.currentcluster,p}.Brange = DATA.cluster{1,p}.Brange;
    DATA.cluster{DATA.currentcluster,p}.Erange = DATA.cluster{1,p}.Erange;
end
arange = DATA.clusterArange;
brange = DATA.clusterBrange;
splen = size(values,2).*size(values,3);
erange = intersect(DATA.clusterErange,1:splen-1);    
    
if ismember(type,[SPKMAXACCEL SPKMINACCEL SPKMEANACCEL])
    acc = diff(values,2,2);
end
%dVdt and values are already just teh values for ispk
%ispk is only sent so that the correct entries in DATA.Spikes are 
%filled in
if type == SPKENERGY
    if recalc
        x  = sum(dVdt(:,erange).^2,2);
        DATA.Spikes.energy(ispk)= x;
    else
        x = DATA.Spikes.energy(ispk);
    end
elseif type == SPKRMSV
        x = sqrt(DATA.Spikes.energy(ispk));
elseif type == ENERGY1
    erange = 1:(splen/4)-1;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY2
    erange = (splen/4):(splen/2)-2;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY3
    erange = (splen/2)-1:3*splen/4-3;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY4
    erange = (splen-splen/4)-3:splen-4;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == SPKENERGYA
    [mins, imins] = min(values'); 
    imins(find(imins<6)) = 6;
    imins(find(imins>26)) = 26;
    for j = 1:size(dVdt,1)
        x(j)  = sum(dVdt(j,imins(j)-5:imins(j)+5).^2,2);
    end
        DATA.Spikes.energy(ispk)= x;
elseif type == SPKVARE
    if recalc
        x = var(values(:,erange)');
        x = x./DATA.Spikes.energy(ispk);
        DATA.Spikes.vw(ispk) = x;
    else
        x = DATA.Spikes.vw(ispk);
    end
elseif type == SPKVAREA
    x = var(values(:,erange)');
    x = x./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKVAREB
    x = var(values(:,erange)');
    x = sqrt(x)./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKMAXRATE
    x = max(dVdt(:,:)');
elseif type == SPKPREMINRATE %min rate value preceding max rate value;
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(dVdt(j,1:t(j)));
    end
elseif type == SPKPREMAXRATE %MAX rate preceding peak V
    [x,t] = max(values(:,:)');
    t(find(t > size(dVdt,2))) = size(dVdt,2);
    for j = 1:length(t)
        x(j) = max(dVdt(j,1:t(j)));
    end
elseif type == SPKMAXRATEA
    x = max(dVdt(:,arange)');
elseif type == SPKMAXRATEB
    x = max(dVdt(:,brange)');
elseif type == SPKMINRATE
    x = min(dVdt(:,:)');
elseif type == SPKMINRATEA
    x = min(dVdt(:,arange)');
elseif type == SPKMINRATEB
    x = min(dVdt(:,brange)');
elseif type == SPKMEANA
    x = mean(values(:,arange)');
elseif type == SPKMEAN
    x = mean(values');
elseif type == SPKMEANB
    x = mean(values(:,brange)');
elseif type == SPKMEANRATEA
    x = mean(dVdt(:,arange)');
elseif type == SPKMEANRATEB
    x = mean(dVdt(:,brange)');
elseif type == SPKMAXA
    x = max(values(:,arange)');
elseif type == SPKPCA1
    x = pcs(:,1);
elseif type == SPKPCA2
    x = pcs(:,2);
elseif type == SPKPCA3
    x = pcs(:,3);
elseif type == SPKPCA4
    x = pcs(:,4);
elseif type == SPKPCA12
    x = pcs(:,1)./pcs(:,2);
elseif type == SPKMAXB
    x = max(values(:,brange)');
elseif type == SPKMIN
    x = min(values(:,:)');
    if DATA.plot.nodc
        x = x - mean(values');
    end
elseif type == SPKMINA
    x = min(values(:,arange)');
elseif type == SPKMINB
    x = min(values(:,brange)');
elseif type == SPKPREMIN
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(values(j,1:t(j)));
    end
elseif type == SPKPEAKTIME
    [x,t] = max(dVdt(:,:)');
    zc = diff(sign(dVdt(:,:)')); 
    len = size(zc,1);
    for j = 1:length(t)
        tm = find(zc(:,j) < 0 & [1:len]' >= t(j));
        if isempty(tm)
            x(j) = 0;
        else
            x(j) = tm(1); % first zero crossing after peakrate
        end
    end
elseif type == SPKPEAK
    [x,t] = max(values(:,:)');
elseif type == SPKHEIGHT
    [x,t] = max(values(:,:)');
    [y,t] = min(values(:,:)');
    x = x-y;
elseif type == SPKVAR
    x = var(values(:,:)');
elseif type == SPKSYMMETRY || type == SPKSYMNEG
  len =  size(dVdt,2);
  cn  = round(sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2));
  for j = 1:length(ispk)
  z(j) = mean(values(j,max([cn(j)-5 1]):cn(j)),2);
  y(j) = mean(values(j,cn(j):min([cn(j)+5 len])),2);
  end
  if type == SPKSYMNEG
  x = atan2(-y,-z);
  else
  x = atan2(y,z);
  end
elseif type == SPKMAXACCEL
    x = max(acc');
elseif type == SPKMINACCEL
    x = min(acc');
elseif type == SPKMEANACCEL
    x = mean(abs(acc'));
elseif type == SPKCENTROID
  len =  size(dVdt,2);
  x  = sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2);
elseif ismember(type, [TEMPLATEA TEMPLATEB TEMPLATEC])
    j = 1+type-TEMPLATEA;
    x = values * DATA.Templates(j,:)';
elseif ismember(type, [SPKADC1 SPKADC2 SPKADC3 SPKADC4])
    j = 1+type-SPKADC1;
    x = values(:,DATA.adcid(j));
end

x = double(reshape(x,1,length(x)));

function [x,DATA] = OldGetSpikeVals(DATA, ispk, type, recalc)

SPKENERGY=1;
SPKVARE=2;
SPKMAXRATE = 3;
SPKPREMINRATE = 4;
SPKPEAKTIME = 5;
SPKPREMIN = 6;
SPKPEAK = 7;
SPKVAR = 8;
SPKSYMMETRY = 9;
SPKCENTROID = 10;
SPKCENTROID = 10;

SPKMINRATE = 11;
SPKMAXRATEA = 12;
SPKMINRATEA = 13;
SPKMEANRATEA = 14;
SPKMAXRATEB = 15;
SPKMINRATEB = 16;
SPKMEANRATEB = 17;
SPKMEANA = 18;
SPKMEANB = 19;
SPKMAXA = 20;
SPKMAXB = 21;
SPKMIN =22;
SPKMINA = 23;
SPKMINB = 24;
SPKHEIGHT = 25;
SPKPREMAXRATE = 26;
SPKVARTESTA= 27;
SPKVARTESTB= 28;


if isnan(ispk) %% return names
    x = {'Energy' 'Var/Energy' 'MaxRate' 'PreMinRate' 'PeakTime' 'PreMin' 'Peak' 'Var' 'Symmetry' 'Centroid' 'Minrate' 'Maxrate(A)' ... 
    'Minrate(A)' 'Meanrate(A)' 'Maxrate(B)' 'Minrate(B)' 'Meanrate(B)' 'Mean(A)'...
    'Mean(B)' 'Max(A)' 'Min(A)' 'Min' 'Max(B)' 'Min(B)' 'Height' 'test1' 'test2'};
  return;
end


p = DATA.probe;
if p > size(DATA.cluster,2)
    DATA.cluster{DATA.currentcluster,p} = {};
end

if DATA.currentcluster > size(DATA.cluster,1) || (~isfield(DATA.cluster{DATA.currentcluster,p},'Arange') & isfield(DATA.cluster{1,p},'Arange'))
    DATA.cluster{DATA.currentcluster,p}.Arange = DATA.cluster{1,p}.Arange;
    DATA.cluster{DATA.currentcluster,p}.Brange = DATA.cluster{1,p}.Brange;
    DATA.cluster{DATA.currentcluster,p}.Erange = DATA.cluster{1,p}.Erange;
end
arange = DATA.clusterArange;
brange = DATA.clusterBrange;
splen = size(DATA.AllData.Spikes.values,2);
erange = intersect(DATA.clusterErange,1:splen-1);    
    

if type == SPKENERGY
    if recalc
        x  = sum(DATA.AllData.Spikes.dVdt(ispk,erange).^2,2);
        DATA.Spikes.energy(ispk)= x;
    else
        x = DATA.Spikes.energy(ispk);
    end

elseif type == SPKVARE
    if recalc
        x = var(DATA.AllData.Spikes.values(ispk,erange)');
        x = x./DATA.Spikes.energy(ispk);
        DATA.Spikes.vw(ispk) = x;
    else
        x = DATA.Spikes.vw(ispk);
    end
elseif type == SPKMAXRATE
    x = max(DATA.AllData.Spikes.dVdt(ispk,:)');
elseif type == SPKPREMINRATE
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    for j = 1:length(t)
        x(j) = min(DATA.AllData.Spikes.dVdt(ispk(j),1:t(j)));
    end
elseif type == SPKPREMAXRATE
    [x,t] = max(DATA.AllData.Spikes.values(ispk,:)');
    for j = 1:length(t)
        x(j) = max(DATA.AllData.Spikes.dVdt(ispk(j),1:t(j)));
    end
elseif type == SPKMAXRATEA
    x = max(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMAXRATEB
    x = max(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMINRATE
    x = min(DATA.AllData.Spikes.dVdt(ispk,:)');
elseif type == SPKMINRATEA
    x = min(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMINRATEB
    x = min(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMEANA
    x = mean(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMEANB
    x = mean(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKMEANRATEA
    x = mean(DATA.AllData.Spikes.dVdt(ispk,arange)');
elseif type == SPKMEANRATEB
    x = mean(DATA.AllData.Spikes.dVdt(ispk,brange)');
elseif type == SPKMAXA
    x = max(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMAXB
    x = max(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKMIN
    x = min(DATA.AllData.Spikes.values(ispk,:)');
elseif type == SPKMINA
    x = min(DATA.AllData.Spikes.values(ispk,arange)');
elseif type == SPKMINB
    x = min(DATA.AllData.Spikes.values(ispk,brange)');
elseif type == SPKPREMIN
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    for j = 1:length(t)
        x(j) = min(DATA.AllData.Spikes.values(ispk(j),1:t(j)));
    end
elseif type == SPKPEAKTIME
    [x,t] = max(DATA.AllData.Spikes.dVdt(ispk,:)');
    zc = diff(sign(DATA.AllData.Spikes.dVdt(ispk,:)')); 
    len = size(zc,1);
    for j = 1:length(t)
        tm = find(zc(:,j) < 0 & [1:len]' >= t(j));
        if isempty(tm)
            x(j) = 0;
        else
            x(j) = tm(1); % first zero crossing after peakrate
        end
    end
elseif type == SPKPEAK
    [x,t] = max(DATA.AllData.Spikes.values(ispk,:)');
elseif type == SPKHEIGHT
    [x,t] = max(DATA.AllData.Spikes.values(ispk,:)');
    [y,t] = min(DATA.AllData.Spikes.values(ispk,:)');
    x = x-y;
elseif type == SPKVAR
    x = var(DATA.AllData.Spikes.values(ispk,:)');
elseif type == SPKSYMMETRY
  len =  size(DATA.AllData.Spikes.dVdt,2);
  cn  = round(sum((DATA.AllData.Spikes.dVdt(ispk,:).^2)*[1:len]',2)./sum(DATA.AllData.Spikes.dVdt(ispk,:).^2,2));
  for j = 1:length(ispk)
  z(j) = mean(DATA.AllData.Spikes.values(ispk(j),max([cn(j)-5 1]):cn(j)),2);
  y(j) = mean(DATA.AllData.Spikes.values(ispk(j),cn(j):min([cn(j)+5 len])),2);
  end
  x = atan2(y,z);
elseif type == SPKCENTROID
  len =  size(DATA.AllData.Spikes.dVdt,2);
  x  = sum((DATA.AllData.Spikes.dVdt(ispk,:).^2)*[1:len]',2)./sum(DATA.AllData.Spikes.dVdt(ispk,:).^2,2);
end


function DATA = PlotSpikes(DATA, ispk)

    classify = 1;
    ctype = 2;
    [cx, DATA] = GetSpikeVals(DATA,ispk, DATA.AllSpikes.values(ispk,:),DATA.AllSpikes.dVdt(ispk,:), DATA.plot.clusterX, classify, []);
    DATA.Spikes.cx(ispk) = cx;
    %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
    nc =  length(DATA.cluster)+1;
    DrawSpikeWaves(DATA, ispk, nc, 2);
    nc =  length(DATA.cluster)+1;
    DrawSpikeWaves(DATA, ispk, nc, 2);
    [cy, DATA] = GetSpikeVals(DATA,ispk, DATA.AllSpikes.values(ispk,:),DATA.AllSpikes.dVdt(ispk,:),DATA.plot.clusterY, classify, []);
    DATA.Spikes.cy(ispk) = cy;
    DATA.currentspike = ispk(1);
    DATA = SetSpkCodes(DATA,ispk,DATA.probe,0);
    nc =  length(DATA.cluster)+1;
    DrawSpikeWaves(DATA, ispk, nc, 2);
    
    set(0,'CurrentFigure',DATA.xyfig);
    for j = 0:nc
        sp = find(DATA.AllData.Spikes.codes(ispk, ctype) == j);
        plot(cx(sp),cy(sp),...
            '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
        hold on; %% need this when called from PlotOneTrial
    end

    function DrawSpikeWaves(DATA, ispk, nclusters, ctype)
if  DATA.svfig == 0
    return;
end
        
        for j = 1:nclusters+1
        vs{j} = [];
        xs{j} = [];
        end
    for j = 1:length(DATA.svh)
        vs{j} = [];
        xs{j} = [];
    end
        needdv = 0;
    
    if isfield(DATA,'AllSpikes')
        splen = size(DATA.AllSpikes{DATA.probe}.values,2).*size(DATA.AllSpikes{DATA.probe}.values,3);
        adc = DATA.AllSpikes{DATA.probe}.values(ispk,:);
        codes = DATA.AllSpikes{DATA.probe}.codes(ispk,:);
        if needdv
        dvdt = DATA.AllSpikes{DATA.probe}.dVdt(ispk,:);
        end
    else
    splen = size(DATA.AllData.Spikes.values,2).*size(DATA.AllData.Spikes.values,3);
    adc = DATA.AllData.Spikes.values(ispk,:);
    codes = DATA.AllData.Spikes.codes(ispk,:);
    dvdt = DATA.AllData.Spikes.dVdt(ispk,:);
    end
    for spk = 1:length(ispk);
        j = codes(spk, ctype)+1;
        if DATA.plot.dvdt
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} [1:splen-1] NaN];
        elseif DATA.plot.voltxy
            vs{j} = [vs{j} adc(spk,33:64) NaN];
            xs{j} = [xs{j} adc(spk,65:96) NaN];
        else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
        end
    end
    set(0,'CurrentFigure',DATA.svfig);
    nc = min([nclusters+1 length(DATA.svh)]);
    for j = 1:max([nc length(DATA.svh)])
        if ~isempty(xs{j})
            if ishandle(DATA.svh(j))
                set(DATA.svh(j),'Xdata' , xs{j}, 'Ydata', vs{j});
            else
                DATA.svh(j) = line('Xdata' , xs{j}, 'Ydata', vs{j});
            end
        elseif ishandle(DATA.svh(j)) %no spikes, wipe clean
            set(DATA.svh(j),'Xdata' , 0, 'Ydata', 0);
        end
    end
    drawnow;

function [ispk, sspk,cx, cy] = PlotTrialSyncSpikes(DATA, times, probes, colors, varargin)

SPKMIN=22;
SPKPEAK=7;
    
    nv = 32;
    dt = 2;
    j=1;
    ctype = 2;
    step = DATA.plot.SpikeMaxV;
    if length(probes) > 3
        step = 3 * DATA.plot.SpikeMaxV/4;
    end
    if DATA.plot.SpikeVsep > 0
        step = DATA.plot.SpikeVsep;
    end
    if DATA.plot.syncoverlay
        timemode = 0;
    else
        timemode = 1;
    end
    while j <= length(varargin)
    if strncmpi(varargin{j},'lineoff',6)
        j = j+1;
        lineoff = varargin{j};
    elseif strncmpi(varargin{j},'Probe',4)
        j = j+1;
        probe = varargin{j};
    elseif strncmpi(varargin{j},'Trial',4)
        j = j+1;
        Trial = varargin{j};
    end
    j = j+1;
    end

    if isfield(DATA,'sids') && isfield(DATA,'AllSpikes')
        ai = find(DATA.AllSpikes{probes(1)}.times(DATA.sids{1}) > times(1) ...
            & DATA.AllSpikes{probes(1)}.times(DATA.sids{1}) < times(2));
        for j = 1:length(probes)
            sids{j} = DATA.sids{j}(ai);
            ispks{j} = find(DATA.AllSpikes{probes(j)}.times > times(1) ...
                & DATA.AllSpikes{probes(j)}.times < times(2));
            itimes{j} = DATA.AllSpikes{probes(j)}.times(sids{j});
        end
        ispk = ispks{1};
    else
    for j = 1:length(probes)
        if DATA.spikelist == -1
            spikelist = [0 1 2 3 4];
        else
            spikelist = DATA.spikelist;
        end
        ispks{j} = find(DATA.AllSpikes{probes(j)}.times > times(1) &...
            DATA.AllSpikes{probes(j)}.times < times(2) & ...
            ismember(DATA.AllSpikes{probes(j)}.codes(:,2),spikelist));
        itimes{j} = DATA.AllSpikes{probes(j)}.times(ispks{j});
    end
    iids = [];
    dts{1} = zeros(size(itimes{1}));
    for p = 2:length(ispks)
        ids = [];
        for j = 1:length(itimes{1})
            id = find(abs(itimes{p}-itimes{1}(j)) < dt);
            if length(id)
                ids = [ids id(1)];
                dts{p}(length(ids)) = (itimes{p}(id(1))-itimes{1}(j)) * nv/10;
                iids = [iids j];
            end
        end
        sids{p} = ispks{p}(ids);
    end
    sids{1} = ispks{1}(iids);
    end
    if ~isempty(sids{1})
        if timemode
            set(0,'CurrentFigure',DATA.timefig);
            vh = DATA.tvh;
        else
            set(0,'CurrentFigure',DATA.svfig);
            vh = DATA.svh;
        end
        nc = (length(probes)-1)/2;
        n = DATA.svhn;
    for j = 1:length(probes)
        for k = 1:n
        vs{(j-1)*n+k} = [];
        xs{(j-1)*n+k} = [];
        end
        dts{j} = itimes{j} - itimes{1};
       for k = 1:length(sids{j});
           spk = sids{j}(k);
           nl = (j-1)*n + DATA.AllSpikes{probes(j)}.codes(spk,2)+1;
           if DATA.plot.syncoverlay
           vs{nl} = [vs{nl} DATA.AllSpikes{probes(j)}.values(spk,:)+(j-nc-1)*step NaN];
            xs{nl} = [xs{nl} [1:nv]-dts{j}(k) NaN];
           else
           vs{nl} = [vs{nl} DATA.AllSpikes{probes(j)}.values(spk,:) NaN];
            xs{nl} = [xs{nl} (k-1)*nv + [1:nv]-dts{j}(k) NaN];
            end
       end

                
    for k = 1:n
        nl = (j-1)*n + k;
    if nl > length(DATA.svh)
        vh(nl) = line('Xdata' , xs{nl}, 'Ydata', vs{nl});
    else
        if ishandle(DATA.svh(j))
            set(vh(nl),'Xdata' , xs{nl}, 'Ydata', vs{nl});
        end
    end
    end
    end
    if length(xs{1})
        set(gca,'xlim',[1 max(xs{1})]);
    end
    end
    sspk = cat(2,sids{:});
   nspk = length(sids{1});

   title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f %d/%d spks',abs(Trial.Trial),...
       Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,nspk,length(ispk)));

   cx = [];
   cy = [];
   if DATA.plot.synccluster == 1 %min/min
        [cx, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),SPKMIN, 0,[]);
        [cy, DATA] = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),SPKMIN, 0,[]);
   elseif DATA.plot.synccluster == 4 %PCA1 spk1 vs 2
        cx = DATA.pca(3,ai);
        cy = DATA.pca(4,ai);
   elseif DATA.plot.synccluster == 5 %PCA2 1 vs 2
        cx = DATA.pca(5,ai);
        cy = DATA.pca(6,ai);
   elseif DATA.plot.synccluster == 6 %PCA 1 vs PCA2
        cx = DATA.pca(1,ai);
        cy = DATA.pca(2,ai);
   elseif DATA.plot.synccluster == 7 %PCA vs xcorr
        cx = DATA.pca(1,ai);
        cy = DATA.AllSpikes{probes(1)}.values(sids{1},:) .* DATA.AllSpikes{probes(2)}.values(sids{2},:);
        cy = sum(cy');
   elseif DATA.plot.synccluster == 8 %Current Cluster X, for each probe
        cx  = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterX, 1,[]);
        cy  = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),DATA.plot.clusterX, 1,[]);
   elseif DATA.plot.synccluster == 9 %Current Cluster Y, for each probe
        cx  = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterY, 1,[]);
        cy  = GetSpikeVals(DATA,sids{2}, DATA.AllSpikes{probes(2)}.values(sids{2},:), DATA.AllSpikes{probes(2)}.dVdt(sids{2},:),DATA.plot.clusterY, 1,[]);
   elseif DATA.plot.synccluster == 10 %Current Cluster Plot for sum
       V = (DATA.AllSpikes{probes(1)}.values(sids{1},:)+DATA.AllSpikes{probes(2)}.values(sids{2},:))/2;
       dV= ( DATA.AllSpikes{probes(1)}.dVdt(sids{1},:)+ DATA.AllSpikes{probes(2)}.dVdt(sids{2},:))/2;
        cx  = GetSpikeVals(DATA,sids{1}, V, dV, DATA.plot.clusterX, 1);
        cy  = GetSpikeVals(DATA,sids{1}, V, dV, DATA.plot.clusterY, 1);
   else
        [cx, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterX, 0,[]);
        [cy, DATA] = GetSpikeVals(DATA,sids{1}, DATA.AllSpikes{probes(1)}.values(sids{1},:), DATA.AllSpikes{probes(1)}.dVdt(sids{1},:),DATA.plot.clusterY, 0,[]);
   end
   if length(cx) & length(cy)
   set(0,'CurrentFigure',DATA.xyfig);
   for j = 0:length(DATA.cluster{probes(1)})+1
       sp = find(DATA.AllSpikes{probes(1)}.codes(sids{1}, ctype) == j);
       plot(cx(sp),cy(sp),...
           '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
       hold on; %% need this when called from PlotOneTrial
   end
   end
   drawnow;
    
function [DATA, ispk] = APlotTrialSpikes(DATA, times, colors, nclusters, classify, varargin)

j = 1;
Trial = [];
probe = DATA.probe;
lineoff = 0;
step = DATA.plot.SpikeMaxV;
step = DATA.vstep;
timemode = 0;
nprobes = 1;
if DATA.plot.showartifacts
    maxcl = 9;
else
maxcl = 8;
end
voff = NaN;
ip = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'lineoff',6)
        j = j+1;
        lineoff = varargin{j};
    elseif strncmpi(varargin{j},'Probes',6)
        j = j+1;
        probe = varargin{j};
    elseif strncmpi(varargin{j},'Probe',4)
        j = j+1;
        probe = varargin{j}(1);
        if length(varargin{j}) > 1
            ip = varargin{j}(2); % # in list
        end
        if length(varargin{j}) > 2
            nprobes = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'timemode',4)
        timemode = 1;
    elseif strncmpi(varargin{j},'Trial',4)
        j = j+1;
        Trial = varargin{j};
    elseif strncmpi(varargin{j},'voff',4)
        j = j+1;
        voff = varargin{j};
    end
    j = j+1;
end
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{probe(1)};
    if isfield(Spks,'pcs')
        PCs = Spks.pcs;
    else
    PCs = [];
    end
elseif DATA.probe == 100
    Spks = DATA.AllData.UstimV;
else
    Spks = DATA.AllData.Spikes;
    if isempty(DATA.AllData.pcs)
        PCs = DATA.AllData.Spikes.codes;
    else
    PCs = DATA.AllData.pcs;
    end
end

if ~isfield(Spks,'values')
    ispk = [];
    return;
end
splen = size(Spks.values,2).*size(Spks.values,3);
% try calculating energy for all spikes in Expt in one step.
ispk = find(Spks.times > times(1) &...
        Spks.times < times(2) & Spks.codes(:,2) < maxcl);
if DATA.TriggerSign == -1
    id = find(Spks.values(ispk,11) < 0);
    ispk = ispk(id);
end
if isfield(DATA,'AllClusters')
    return;
end
if ispk
    if length(PCs) < max(ispk)
        PCs = Spks.codes;
    end
    if isfield(DATA,'sids')
        syncspk = find(ismember(ispk,DATA.sids{ip}));
        syncspklst = ismember(ispk,DATA.sids{ip});
        if length(probe) == 2
            Spks.values = (DATA.AllSpikes{probe(1)}.values(DATA.sids{1},:)+DATA.AllSpikes{probe(2)}.values(DATA.sids{2},:))/2;
            Spks.times =  DATA.AllSpikes{probe(1)}.times(DATA.sids{1});
            ispk = find(Spks.times > times(1) &...
                Spks.times < times(2));
        end
        if ismember(DATA.syncsign,[3 4 5])
            Spks.codes(ispk,3) = 0;
            Spks.codes(ispk(syncspk),3) = 1;
            ctype = 3;
        end
    end

    if DATA.probe < 100
        [cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs(ispk,:));
        DATA.Spikes.cx(ispk) = cx;
        %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
        [cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs(ispk,:));
        DATA.Spikes.cy(ispk) = cy;
        dvdt = Spks.dVdt(ispk,:);
        if length(ispk)
            DATA.currentspike = ispk(1);
        end
        adc = Spks.values(ispk,:);
        % recut == 2 means that clusters are not set here, but clusters have
        % been defined (previous list), so use those properites.
        if DATA.state.recut == 2 && DATA.probe == probe
            DATA = SetSpkCodes(DATA,ispk, DATA.probe,0);
            if isfield(DATA,'AllSpikes') %copy any changes into spks
                Spks.codes = DATA.AllSpikes{DATA.probe}.codes;
            else
                Spks.codes = DATA.AllData.Spikes.codes;
            end
        end
    else
        Spks.codes = DATA.AllData.UstimV.codes;
        adc = Spks.values(ispk,:);
        cx = zeros(size(ispk));
        cy = zeros(size(ispk));
    end


    for j = [1:nclusters+1]
        vs{j} = [];
        xs{j} = [];
    end
    for spk = 1:length(ispk);
        if DATA.syncsign < 2 & isfield(DATA,'sids')
            j = syncspklst(spk)+1;
        else
            j = Spks.codes(ispk(spk), ctype)+1;
        end
        if DATA.plot.dvdt ==2
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} adc(spk,1:splen-1) NaN];
        elseif DATA.plot.dvdt
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} [1:splen-1] NaN];
        elseif timemode
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen]+(Spks.times(ispk(spk))-times(1)).*timemode NaN];
        elseif DATA.plot.nodc
            vs{j} = [vs{j} adc(spk,:)-mean(adc(spk,:)) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
        elseif DATA.plot.voltxy == 3
            vs{j} = [vs{j} adc(spk,65:96) NaN];
            xs{j} = [xs{j} adc(spk,97:128) NaN];
        elseif DATA.plot.voltxy == 2
            vs{j} = [vs{j} adc(spk,33:64) NaN];
            xs{j} = [xs{j} adc(spk,65:96) NaN];
        elseif DATA.plot.voltxy == 1
            vs{j} = [vs{j} adc(spk,1:32) NaN];
            xs{j} = [xs{j} adc(spk,33:64) NaN];
        elseif DATA.plot.voltxy == 4
            vs{j} = [vs{j} adc(spk,65:96)-adc(spk,33:64) NaN];
            xs{j} = [xs{j} [1:splen/4] NaN];
        else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
        end
    end
    if length(DATA.plot.voffsets) >= ip
        for j = [1:nclusters+1]
            vs{j} = vs{j}+DATA.plot.voffsets(ip);
        end
        DATA.plot.currentoffset(ip) = DATA.plot.voffsets(ip);
    else
        DATA.plot.currentoffset(ip) = (ip-(nprobes+1)/2).*step;
        for j = [1:nclusters+1]
            vs{j} = vs{j}+(ip-(nprobes+1)/2).*step;
        end
    end

    if timemode && isfield(DATA,'timefig')
        set(0,'CurrentFigure',DATA.timefig);
        text(0,(ip-(nprobes+1)/2).*step,sprintf('%d',probe));
        vh = DATA.tvh;
    else
        set(0,'CurrentFigure',DATA.svfig);
        vh = DATA.svh;
    end
    nc = min([nclusters+1 length(DATA.svh)]);


    for j = 1:nc
        k = j+lineoff;
        if k > length(vh)
            vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
        elseif ~isempty(xs{j})
            if ishandle(vh(k))
                set(vh(k),'Xdata' , xs{j}, 'Ydata', vs{j});
            else
                vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
            end
        elseif vh(k) && ishandle(vh(k)) %no spikes, wipe clean
            set(vh(k),'Xdata' , 0, 'Ydata', 0);
        end
    end
    if ~isempty(Trial)
        xc = ExtraLabels(Trial);
        nspk = sum(Spks.codes(ispk,2) == DATA.currentcluster);
        title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f%s %d/%d spks',abs(Trial.Trial),...
            Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,xc,nspk,length(ispk)));
    end

    if probe == DATA.probe
        drawnow;
        set(0,'CurrentFigure',DATA.xyfig);
        for j = 0:nclusters
            sp = find(Spks.codes(ispk, ctype) == j);
            plot(cx(sp),cy(sp),...
                '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
            hold on; %% need this when called from PlotOneTrial
        end
    end
else
        title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f',abs(Trial.Trial),...
            Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed));
    
end
    DATA.minplottime = 0.00;
    if DATA.minplottime > 0.001
        while toc < DATA.minplottime
        end
    end
            

function xc = ExtraLabels(Trial)
    xc = '';
    if isfield(Trial,'uStim') && Trial.uStim
        xc = [xc 'Ustm '];
    end

function spikelist = WhichClusters(top,varargin)

spikelist = [];

if get(findobj(top,'Tag','UseCluster0'),'value')
    spikelist = [spikelist 0];
end
if get(findobj(top,'Tag','UseCluster1'),'value')
    spikelist = [spikelist 1];
end
if get(findobj(top,'Tag','UseCluster2'),'value')
    spikelist = [spikelist 2];
end
if get(findobj(top,'Tag','UseCluster3'),'value')
    spikelist = [spikelist 3];
end
if get(findobj(top,'Tag','UseCluster4'),'value')
    spikelist = [spikelist 4];
end
if isempty(spikelist)
    spikelist = -1;
end


function [meanrate, proportion] = SumExptSpikes(DATA, expid, cluster)
    cspks = DATA.Expts{expid}.gui.spks;
    ctype = 2;
    j = 1;
    for trial = [DATA.Expts{expid}.Trials]
        if isfield(DATA,'AllClusters')
        ispk = find(DATA.AllClusters(DATA.probe).times(cspks) > trial.Start(1) & DATA.AllClusters(DATA.probe).times(cspks) < trial.End(end)+500);
        ispks(j) = sum(DATA.AllClusters(DATA.probe).codes(cspks(ispk)) == cluster);
        elseif isfield(DATA,'AllSpikes')
        ispk = find(DATA.AllSpikes{DATA.probe}.times(cspks) > trial.Start(1) & DATA.AllSpikes{DATA.probe}.times(cspks) < trial.End(end)+500);
        ispks(j) = sum(DATA.AllSpikes{DATA.probe}.codes(cspks(ispk),ctype) == cluster);
        else
        ispk = find(DATA.AllData.Spikes.times(cspks) > trial.Start(1) & DATA.AllData.Spikes.times(cspks) < trial.End(end)+500);
    ispks(j) = sum(DATA.AllData.Spikes.codes(cspks(ispk),ctype) == cluster);
        end
    aspks(j) = length(ispk);
    dur(j) = trial.End(end)-trial.Start(1);
    j = j+1;
    end
    meanrate = 10000.* mean(ispks)./mean(dur);
    proportion = sum(ispks)./sum(aspks);
    
function [DATA, counts] = CountSpikes(DATA, expid, varargin)

    if isfield(DATA,'AllClusters')
%can be using AllData.Spikes, but have no spikes
    elseif ~isempty(DATA.AllData.Spikes) && length(DATA.AllData.Spikes.codes) == 0
        counts = 0;
        return;
    end
replot = 0;
reclassify = DATA.state.recount;
j = 1;
while j <= nargin -2
    if strncmp(varargin{j},'reclassify',3)
        reclassify = 1;
    elseif strncmp(varargin{j},'replot',3)
        replot = 1;
    end
    j = j+1;
end

nt = 1;
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
if  DATA.state.online == 2
return;
end
if DATA.Expts{expid}.gui.classified == 0  || reclassify
   DATA = SetExptSpikes(DATA, expid, 0);
end

if isfield(DATA,'AllSpikes')
Spks = DATA.AllSpikes{DATA.probe};
elseif isfield(DATA,'AllClusters')
Spks = DATA.AllClusters(DATA.probe);
else
Spks = DATA.AllData.Spikes;
end

if isempty(Spks) | ~isfield(Spks,'codes')
    for nt = 1:length(DATA.Expts{expid}.Trials)
        DATA.Expts{expid}.Trials(nt).count = 0;
        DATA.Expts{expid}.Trials(nt).Spikes = [];
    end
        
    counts = 0;
    return;
end
spikelist = WhichClusters(DATA.toplevel);
alli = [];
%restrict finds to the range of this expt - speeds things up
espk = ExptSpikeListAll(DATA,expid,Spks.times);
nt = 1;
for trial = [DATA.Expts{expid}.Trials]
    ispk = find(Spks.times(espk) > trial.Start(1)-DATA.state.preperiod & Spks.times(espk) < trial.End(end)+DATA.state.postperiod);
    ispk = espk(ispk);
    alli = [alli; ispk];
    ispks = find(ismember(Spks.codes(ispk,ctype),spikelist));
    DATA.Expts{expid}.Trials(nt).Spikes = round(Spks.times(ispk(ispks)) - trial.Start(1));
    spks = DATA.Expts{expid}.Trials(nt).Spikes;
    DATA.Expts{expid}.Trials(nt).count = sum(spks > 500);
    ispks = find(~ismember(Spks.codes(ispk,ctype),spikelist));
    DATA.Expts{expid}.Trials(nt).OSpikes = round(Spks.times(ispk(ispks)) - trial.Start(1));
    if DATA.state.recut
        DATA.Expts{expid}.Trials(nt).Ocodes = Spks.codes(ispk(ispks),2);
    else
        DATA.Expts{expid}.Trials(nt).Ocodes = Spks.codes(ispk(ispks),1);
    end
    nt = nt+1;
end
if replot
    GetFigure(DATA.tag.dataplot);
    args = PlotArgs(DATA, DATA.Expts{expid});
    PlotExpt(DATA.Expts{expid},args{:});
    t = get(get(gca,'title'),'String');
    title([t sprintf(' P%d',DATA.probe) sprintf(' Cl %d',spikelist)]);
end
DATA.Expts{expid}.gui.counted = sum(spikelist);
if ~isfield(DATA.Expts{expid}.gui,'setispk') || length(alli) ~= DATA.Expts{expid}.gui.setispk
DATA.Expts{expid}.gui.spks = alli;
DATA.Expts{expid}.gui.setispk = 2;
end
DATA.spikelist = spikelist;

nc = max(Spks.codes(alli,2));
counts(1) = sum(Spks.codes(alli,2) == 0);
if nc
for j = 1:nc
    counts(j+1) = sum(Spks.codes(alli,2) == j);
end
end
DATA.Expts{expid}.gui.spkcounts = counts;
DATA.Expts{expid}.gui.nclusters = nc;

    
function espk = ExptSpikeListAll(DATA, eid, stimes)
    t(1) = DATA.Expts{eid}.Trials(1).Start(1)-DATA.state.preperiod;
    t(2) = DATA.Expts{eid}.Trials(end).End(end)+DATA.state.postperiod;
    espk = find(stimes > t(1) & stimes < t(2));


function nc = CountClusters(Cluster)
    nc = 0;
if isempty(Cluster)
    nc = 0;
else
    for k = 1:size(Cluster,2) 
        ncs(k) = 0;
    for j = 1:min([size(Cluster,1) 7])  %Clusters > 7 are artifact
        if isfield(Cluster{j,k},'x')
            ncs(k) = j;
        end
    end
    end
    nc = max(ncs);
end

function outname = Expt2FileName(DATA,Expt, cluster)

    [expname, extype, suff] = Expt2Name(Expt);
    it = strmatch(extype,DATA.expstrs,'exact');
    if length(it)
        expname = strrep(expname,extype,DATA.expnames{it(1)});
    end
    if DATA.state.includeprobename
        cs = ['.p' num2str(DATA.probe) 'c'];
    else
        cs = '.c';
    end
    outname = [strrep(DATA.datafilename,'.mat',cs) num2str(cluster) '.' expname suff '.mat'];
    
    
function outname = CombinedName(DATA,eid,cluster, varargin)

probe = DATA.probe;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'probe',5)
        j = j+1;
        probe = varargin{j};
    end
    j = j+1;
end
    
suff = '';
it = strmatch(DATA.exptypelist{eid(1)},DATA.expstrs,'exact');
if ~isempty(it)
    expname = DATA.expnames{it};
else
    expname = DATA.exptypelist{eid(1)};
end
stimname = strrep(DATA.explist{eid(1)},['.' DATA.exptypelist{eid(1)}],'');
if strfind(stimname,'CRC')
    stimname = strrep(stimname,'CRC','');
    suff = 'CRC';
elseif strfind(stimname,'RC')
    stimname = strrep(stimname,'RC','');
    suff = 'RC';
end
if isempty(strfind(expname,stimname))
    expname = [stimname '.' expname];
end
if DATA.state.includeprobename
    cs = ['.p' num2str(probe) 'c'];
else
    cs = '.c';
end
if DATA.state.online
outname = [DATA.datafilename '/' expname suff cs num2str(cluster) '.mat'];
else
outname = [strrep(DATA.datafilename,'.mat',cs) num2str(cluster) '.' expname suff '.mat'];
end


function savecl = SetExptClusters(caller,b, varargin)
%DATA = combine('getstate');

savecl = 1;
DATA = GetDataFromFig(caller);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nosave',5)
        if length(DATA.probes) > 2
            savecl = 0;
        end
    end
    j = j+1;
end

nc = CountClusters(DATA.cluster);
DATA.nclusters = nc;
eid = DATA.currentexpt;
DATA.Expts{eid}.Cluster{1,DATA.probe}.touched = 2;
DATA.Expts{DATA.currentexpt}.gui.ncluster = nc;
DATA.cluster{DATA.currentcluster,DATA.probe}.autocut = 0;
if isfield(DATA.cluster{DATA.currentcluster,DATA.probe},'deleted') & ...
    DATA.cluster{DATA.currentcluster,DATA.probe}.deleted == 1
    DATA.cluster{DATA.currentcluster,DATA.probe}.quality = 0;
else
    DATA.cluster{DATA.currentcluster,DATA.probe}.quality = DATA.state.currentclusterquality;
end

touched = zeros(size(DATA.probelist));
if isfield(DATA,'AllClusters') || isfield(DATA,'AllSpikes') % only save clusters that have been set
for k = 1:nc
    for j = 1:min([size(DATA.cluster,2) length(DATA.AllSpikes)])
        if  isfield(DATA.cluster{k,j}, 'touched') && DATA.cluster{k,j}.touched > 0 && isfield(DATA.AllSpikes{j},'cx')
        DATA.Expts{DATA.currentexpt}.Cluster{k,j} = DATA.cluster{k,j};
        if sum(ismember(DATA.cluster{k,j}.params,[29 30 31])) %templates used
            DATA.Expts{DATA.currentexpt}.Cluster{k,j}.Templates = DATA.Templates;
        end
          touched(j) = 1;
        end
    end
end
    fprintf('Probes %s\n',num2str(find(touched >0)));
else
        DATA.Expts{DATA.currentexpt}.Cluster = DATA.cluster;
        for k = 1:nc
        if isfield(DATA.cluster{k,1},'params') & sum(ismember(DATA.cluster{k,1}.params,[29 30 31])) %templates used
            DATA.Expts{DATA.currentexpt}.Cluster{k,1}.Templates = DATA.Templates;
        end
        end
%
%  If the set button is hit, it means that user is happy with cluster
%  currently shown in the XY plot. Need to track this so that ONLY these
%  clusters are recorded in OnlineClusters. Set touched to 2, then onlhy
% Add these to AllClusters below
    DATA.Expts{eid}.Cluster{1,DATA.probe}.touched = 2;
end
DATA.Expts{DATA.currentexpt}.gui.clustertype = 1;
DATA.Expts{DATA.currentexpt}.gui.classified = 1;
spkvarnames= DATA.spkvarnames;
p = DATA.probe;
for j = 1:length(DATA.Expts)
    if isfield(DATA.Expts{j},'Cluster') 
        enc =  CountClusters(DATA.Expts{j}.Cluster);
        if DATA.state.online == 1
            for k = 1:size(DATA.Expts{j}.Cluster,2)
%             for k = p;
                if ~isfield(DATA.Expts{j}.Cluster{1,k},'touched')
                    DATA.Expts{j}.Cluster{1,k}.touched = 0;
                end
%online temp cuts "just top see" get save here. User needs to delete diligenlyt
%tried only saving current probe, but that is no good, esp for Grid Data
                if DATA.Expts{j}.Cluster{1,k}.touched > 0
                    AllClusters{j}.Cluster(1:enc,k) = DATA.Expts{j}.Cluster(1:enc,k);
                    touched(p) = 1;
                end
            end
        else
        AllClusters{j}.Cluster = DATA.Expts{j}.Cluster(1:enc,:);
        end
        AllClusters{j}.ids = [DATA.Expts{j}.Trials(1).id DATA.Expts{j}.Trials(end).id];
        clustertypes(j) = DATA.Expts{j}.gui.clustertype;
    end
    excludelist{j}  = find([DATA.Expts{j}.Trials.Trial] < 0);
    if isfield(DATA.Expts{j}.gui,'spks')
        expispk{j} = DATA.Expts{j}.gui.spks;
    else
        expispk{j} = [];
    end
    ispklen(j) = length(expispk{j});
end

if savecl
    a = findobj(DATA.xyfig,'String','Set+Next');
    set(a,'backgroundcolor','w')
    DATA.savedclusters = 1;
else
    DATA.savedclusters = 0;
end
Templates = DATA.Templates;
if isfield(DATA,'TemplateInfo')
TemplateInfo = DATA.TemplateInfo;
else
TemplateInfo = [];
end

if DATA.currentcluster > 7 %defining artifact
    id = PlotArtifacts(DATA);
    if length(id)
        if isfield(DATA,'AllSpikes')
        elseif isfield(DATA.AllData,'Spikes')
            DATA.AllData.Spikes.codes(id,2) = 8;
        end
    end
end
if isfield(DATA,'AllSpikes')  %need to save all probes with a new cluster
    for j = find(touched > 0)
        if isfield(DATA.AllSpikes{j},'times')
        f = DATA.AllSpikes{j}.firstspki;
        clid(f:f+size(DATA.AllSpikes{j}.times,1)-1) = DATA.AllSpikes{j}.codes(:,2);
        save(ClusterFile(DATA,'probe',j),'AllClusters','clustertypes','excludelist','clid','spkvarnames','Templates','TemplateInfo');
        end
    end
    if strncmp(DATA.filetype,'Grid',4)
        for j = 1:length(DATA.AllSpikes)
            if isempty(DATA.AllSpikes{j})
                codes{j} = [];
                times{j} = [];
            else
                codes{j} = DATA.AllSpikes{j}.codes(:,2);
                times{j} = DATA.AllSpikes{j}.times;
            end
        end
        save(ClusterFile(DATA,'codes'),'codes','times');
    end
elseif isfield(DATA,'AllClusters')  %need to save all probes with a new cluster
    for j = find(touched > 0)
        clid = DATA.AllClusters(j).codes(:,2);
        save(ClusterFile(DATA),'AllClusters','clustertypes','excludelist','clid','spkvarnames','Templates','TemplateInfo');
    end
elseif isfield(DATA.AllData.Spikes, 'codes')
    if savecl
        clid = DATA.AllData.Spikes.codes(:,2);
        cfile = ClusterFile(DATA);
        save(cfile,'AllClusters','clustertypes','excludelist','clid','spkvarnames','ispklen','Templates','TemplateInfo');
        ifile = regexprep(cfile,'\.p([0-9]*)cl\.','.p$1ispk.');
        if ~exist(ifile)
            save(ifile,'expispk','ispklen')
        end
    end
else
save(ClusterFile(DATA),'AllClusters','clustertypes','excludelist','clid','spkvarnames','ispklen','Templates','TemplateInfo');
end

if length(DATA.probelist) > 1 && savecl > 1
save(ClusterFile(DATA,'allprobes'),'AllClusters','spkvarnames','Templates','TemplateInfo');
end

if DATA.logfid
    fprintf(DATA.logfid,'Cluster set for P%d Expt %d %s by %s\n',DATA.probe,DATA.currentexpt,datestr(now),DATA.user);
end
% now re-do list of su-expts to reflect cut clusters
eid = get(DATA.clst,'value');
if DATA.currentcluster > 7 %defining artifact
    DATA.currentcluster = 1;
    figure(DATA.xyfig);
    hold off;
    DATA = DrawXYPlot(DATA, DATA.Expts{DATA.currentexpt}.gui.spks);
    cid = findobj('Tag','Clusterid');
    set(cid,'value',1);
end
DATA = ListSubExpts(DATA,eid,'relist');
set(DATA.toplevel,'UserData',DATA);
cid = findobj('Tag','ClusterIsSet');
set(cid,'value',1);
if DATA.state.autoplotcells
    GetFigure(DATA.tag.celllist);
    PlotCellList(DATA);
end


function DelClusterButton(caller,b)
it = findobj(get(caller,'parent'),'Tag','Clusterid');
c = get(it,'value');
DeleteCluster(c, caller);

function ClrSpkWin(caller,b)
%DATA = combine('getstate');
if isstruct(caller)
    DATA = caller;
else
    DATA = GetDataFromFig(caller);
end
GetFigure(DATA.tag.clusterxy);
ym = get(gca,'ylim');
xm = get(gca,'xlim');
hold off;
plot(0,0,'+');
set(gca,'ylim',ym);
set(gca,'xlim',xm);
hold on;
if ~isstruct(caller) % from a mouse button
DATA = DrawClusters(DATA,DATA.cluster, 0);
set(DATA.toplevel,'UserData',DATA);
end


function cfile = ClusterFile(DATA,varargin)
getonline = 0;
getauto = 0;
allprobes = 0;
spkcodefile = 0;
probe = DATA.probe;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'auto',4) %Read Clusters from Online mat file
        getauto = 1;
    elseif strncmpi(varargin{j},'allprobes',4) %Read Clusters from Online mat file
        allprobes = 1;
    elseif strncmpi(varargin{j},'getonline',4) %Read Clusters from Online mat file
        getonline = 1;
    elseif strncmpi(varargin{j},'codes',4) %Save codes/times
        spkcodefile = 1;
    elseif strncmpi(varargin{j},'probe',4) %Read Clusters from Online mat file
        j = j+1;
        probe = varargin{j};
    end
    j = j+1;
end
if DATA.state.online == 0 %%not online data
    if getonline
        [a,b] = splitpath(DATA.datafilename);
        if length(b)
        cfile = [b '/OnlineClusters.mat'];
        else
        cfile = 'OnlineClusters.mat';
        end
    elseif getauto
        cfile = strrep(DATA.datafilename,'.mat','.autocl.mat');
    elseif length(DATA.probelist) > 1 && allprobes
        cfile = strrep(DATA.datafilename,'.mat','.allcl.mat');
    elseif spkcodefile 
        cfile = strrep(DATA.datafilename,'.mat','codes.mat');
    elseif length(DATA.probelist) > 1 
        cfile = strrep(DATA.datafilename,'.mat',['.p' num2str(probe) 'cl.mat']);
    else
        cfile = strrep(DATA.datafilename,'.mat','.cl.mat');
    end
else
    if spkcodefile 
        cfile = sprintf('%s/SpikeCodes%d.mat',DATA.datafilename,DATA.Expts{DATA.currentexpt}.Trials(end).id);
    else
        cfile = [DATA.datafilename '/Clusters.mat'];
    end
end

function cfile = CombinerLst(DATA)
if isfield(DATA,'datafilename') %%not online data
    cfile = strrep(DATA.datafilename,'.mat','.combine.mat');
else
    cfile = [DATA.datafilename '/combine.mat'];
end

function DATA = CalcClusterVars(DATA, ispk, varargin)
SPKENERGY=1;
SPKVARE = 2;
noforce  = 1;
j = 1;
probe = DATA.probe;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j},'values')
        j = j+1;
        Spikes = varargin{j};
    elseif strncmpi(varargin{j},'force',5)
        noforce = 0;
    elseif strncmpi(varargin{j},'probe',5)
        j =j+1;
        probe = varargin{j};
    end
    j = j+1;
end

if isfield(DATA,'AllClusters') && noforce % All calculated and stored, but No ADCs, so
  return;
end

if length(ispk) == 2
    ispk = ispk(1):ispk(2);
end
if ispk
    DATA = CheckForPCA(DATA,ispk, 1);
    if isfield(DATA,'AllSpikes')
        Spikes = DATA.AllSpikes{probe};
        if isfield(DATA.AllSpikes{probe},'pcs')
        PCs = DATA.AllSpikes{probe}.pcs;
        else
        PCs = DATA.AllSpikes{probe}.codes;
        end
        adc = Spikes.values(ispk,:);
        if ~isfield(Spikes,'dVdt')
            Spikes.dVdt(ispk,:) = diff(adc,[],2);
        end
    else
        Spikes = DATA.AllData.Spikes;
        adc = Spikes.values(ispk,:);
        if isempty(DATA.AllData.pcs) || length(DATA.AllData.pcs) < max(ispk) %might have pcs computed for earlier spikes
            PCs = DATA.AllData.Spikes.codes;
        else
        PCs = DATA.AllData.pcs;
        end
    end
    adc = Spikes.values(ispk,:);
    if isfield(Spikes,'dVdt')
        energy  = sum(Spikes.dVdt(ispk,:)'.^2);
    else
    energy  = sum(diff(adc').^2);
    end
    svar = var(adc');
    DATA.Spikes.energy(ispk)= energy;
    if DATA.plot.clusterX == SPKENERGY && noforce
        DATA.Spikes.cx(ispk)= energy;
    else
        DATA.Spikes.cx(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterX, 1,PCs(ispk,:));
    end
    DATA.Spikes.vw(ispk) = svar./energy;
    if DATA.plot.clusterY == SPKVARE & noforce
        DATA.Spikes.cy(ispk)= DATA.Spikes.vw(ispk);
    else
        DATA.Spikes.cy(ispk)= GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),DATA.plot.clusterY, 1,PCs(ispk,:));
    end
    if isfield(DATA,'AllSpikes')
        lastspk = min([length(DATA.AllSpikes{probe}.times) length(DATA.Spikes.cx)]);
        DATA.AllSpikes{probe}.cx = DATA.Spikes.cx(1:lastspk);
        DATA.AllSpikes{probe}.cy = DATA.Spikes.cy(1:lastspk);
    end
end

function [x,y, DATA] = GetClusterSpace(DATA, Expt)


    x = DATA.plot.clusterX;
    y = DATA.plot.clusterY;
    p = DATA.probe;
    if isfield(Expt,'Cluster') & ~isempty(Expt.Cluster) & size(Expt.Cluster,2) >= p
        j = 1;
        while j < size(Expt.Cluster,1) && isempty(Expt.Cluster{j,p})
            j = j+1;
        end
        if j <= size(Expt.Cluster,1) & p <= size(Expt.Cluster,2) & ~isempty(Expt.Cluster{j,p}) & isfield(Expt.Cluster{j,p},'params')
            x = Expt.Cluster{j,p}.params(1);
            y = Expt.Cluster{j,p}.params(2);
            DATA.clusterArange = Expt.Cluster{j,p}.Arange;
            DATA.clusterBrange = Expt.Cluster{j,p}.Brange;
            if isfield(Expt.Cluster{j,p},'Erange')
                DATA.clusterErange = Expt.Cluster{j,p}.Erange;
            end
        end
    else
        x = DATA.plot.clusterX;
        y = DATA.plot.clusterY;
    end
       
function set = ClusterIsSet(Expt, probe)

if isfield(Expt,'Cluster') && size(Expt.Cluster,2) >= probe && ...
        ~isempty(Expt.Cluster{1,probe}) && isfield(Expt.Cluster{1,probe},'x')
    set = 1;
    if isfield(Expt.Cluster{1,probe},'autocut') & Expt.Cluster{1,probe}.autocut > 0
        set = 2;
    end
else
    set = 0;
end
   
function [ispk, spktimes, codes] = FindSpikes(DATA, times, probe, range)

if DATA.plot.showartifacts
    maxcl = 9;
else
    maxcl = 8;  %>8 is artifacts. Set this to 9 to include artifacts
end
if isfield(DATA,'AllSpikes')
    if ~isfield(DATA.AllSpikes{probe},'codes') | isempty(DATA.AllSpikes{probe}) | ~isfield(DATA.AllSpikes{probe},'times')
        ispk = [];
        return;
    end
    ispk = find(DATA.AllSpikes{probe}.times > times(1) &...
    DATA.AllSpikes{probe}.times < times(2) & ...
    DATA.AllSpikes{probe}.codes(:,2) < maxcl);
    if nargout > 1
    spktimes = DATA.AllSpikes{probe}.times(ispk);
    codes = DATA.AllSpikes{probe}.codes(ispk,2);
    end
elseif isfield(DATA,'AllClusters')
    ispk = find(DATA.AllClusters(probe).times > times(1) &...
    DATA.AllClusters(probe).times < times(2));
elseif isempty(DATA.AllData.Spikes) || ~isfield(DATA.AllData.Spikes,'times')
    ispk = [];
else
    if isempty(range)
        ispk = find(DATA.AllData.Spikes.times > times(1) &...
            DATA.AllData.Spikes.times < times(2) ...
            & DATA.AllData.Spikes.codes(:,2) < maxcl);
    else
        ispk = find(DATA.AllData.Spikes.times(range) > times(1) &...
            DATA.AllData.Spikes.times(range) < times(2)...
            & DATA.AllData.Spikes.codes(range,2) < maxcl);
        ispk = range(ispk);
    end
end
        

function [DATA, ispk, dprime, details] = SetExptSpikes(DATA, expid, show, varargin)

    details = [];
    dprime = 0;
    ispk = [];
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'spks',4)
            j = j+1;
            ispk = varargin{j};
        elseif strncmpi(varargin{j},'useexpall',9)
            if isfield(DATA.Expts{expid},'Cluster')
                DATA.cluster = DATA.Expts{expid}.Cluster;
            else
                DATA.cluster = [];
            end
            if isfield(DATA.Expts{expid}.gui,'spks')
                ispk = DATA.Expts{expid}.gui.spks;
            end
        elseif strncmpi(varargin{j},'useexp',5)
%            if isfield(DATA.Expts{expid},'Cluster') &
%            size(DATA.Expts{expid}.Cluster,2) >= DATA.probe & size(DATA.Expts{expid}.Cluster,1) > 0 & ~isfield(DATA.Expts{expid}.Cluster{1,DATA.probe},'touched')
% Sep 2010 why the check for touched? Messes up showing existng expt
% clusters
if isfield(DATA.Expts{expid},'Cluster') & size(DATA.Expts{expid}.Cluster,2) >= DATA.probe & size(DATA.Expts{expid}.Cluster,1) > 0 
                DATA.cluster = DATA.Expts{expid}.Cluster;
            end
        end
        j = j+1;
    end
    if isempty(ispk) & ~DATA.state.psychonly
        if isfield(DATA,'AllSpikes')
            espk = [];
        elseif isfield(DATA,'AllClusters')
            espk = [];
        else
            espk = ExptSpikeListAll(DATA, expid, DATA.AllData.Spikes.times);
        end
    for j = 1:length(DATA.Expts{expid}.Trials)
        times(1) = DATA.Expts{expid}.Trials(j).Start(1)-DATA.state.preperiod;
        times(2) = DATA.Expts{expid}.Trials(j).End(end)+DATA.state.postperiod;
        tspk = FindSpikes(DATA, times, DATA.probe, espk);
        ispk = [ispk; tspk];
    end
    end
DATA.Expts{expid}.gui.spkrange = [min(ispk) max(ispk)];
DATA.Expts{expid}.gui.spks = ispk;
if (ischar(show) & strcmp(show,'setrange')) || DATA.state.nospikes;
    return;
end

if ~isfield(DATA,'spklist') | isempty(DATA.spklist)
    DATA.spklist = ispk;
end

DATA = CalcClusterVars(DATA,  ispk);
% if a cluster is set for this expt, use it.
% otherwise use the current one
cl = DATA.currentcluster;
isset = ClusterIsSet(DATA.Expts{expid}, DATA.probe); 
if isset == 1
    if ~isempty(DATA.Expts{expid}.gui.spkrange) && iscluster(DATA.Expts{expid}.Cluster,cl,DATA.probe) && ...
        DATA.Expts{expid}.Cluster{cl,DATA.probe}.lastspk < DATA.Expts{expid}.gui.spkrange(1)
        DATA.Expts{expid}.Cluster{cl,DATA.probe}.lastspk = DATA.Expts{expid}.gui.spkrange(2);
    end
    for k = 1:size(DATA.Expts{expid}.Cluster,2);
        for m = 1:size(DATA.Expts{expid}.Cluster,1);
            if ~isempty(DATA.Expts{expid}.Cluster{m,k})
            DATA.cluster{m,k} = DATA.Expts{expid}.Cluster{m,k};
            end
        end
    end
    DATA.Expts{expid}.gui.classified = 1;
elseif size(DATA.cluster,2) < DATA.probe
    DATA.cluster{1,DATA.probe} = NewClusterStruct(DATA);
elseif iscluster(DATA.cluster,1,DATA.probe) 
% actually is a cluster defined, but not yet applied to this expt
% make sure spk range is correct
    DATA.Expts{expid}.gui.classified = 2;
    if length(ispk)
    DATA.cluster{1,DATA.probe}.firstspk = ispk(1);
    DATA.cluster{1,DATA.probe}.lastspk = ispk(end);
    else
    DATA.cluster{1,DATA.probe}.firstspk = NaN;
    DATA.cluster{1,DATA.probe}.lastspk = NaN;
    end
end
if DATA.state.recut
    DATA.currentexpt = expid;
    DATA = CheckForPCA(DATA, ispk, 0);
    [DATA, dprime, details] = SetSpkCodes(DATA,ispk,DATA.probe,show);
end


function C = NewClusterStruct(DATA)
        C.firstspk = NaN;
        C.lastspk = NaN;

function DATA = DrawClusters(DATA, cluster, setfig)
if setfig
    set(0,'CurrentFigure',DATA.xyfig);
end
p = DATA.probe;
for j = 1:min([size(cluster,1) 7])
    if p > size(cluster,2) 
        C = [];
    else
        C = cluster{j,p};
    end
    while ~isempty(C)
    if isfield(C,'params')
        color = DATA.spkcolor{j+1};
        if isfield(C,'touched') 
            if C.touched == 0
            color = color./2;
            elseif C.touched == -1 %online cluster
               color = (1- color)/2;
            end
        end
        if C.params(1) == DATA.plot.clusterX & ...
                C.params(2) == DATA.plot.clusterY
            h = DrawCluster(C,color);
            if(h)
                DATA.cluster{j,p}.h = h;
                if isfield(C,'autocut') & C.autocut
                    set(h,'LineStyle','--');
                end
            end
            hold on;
        end
    end
    if isfield(C,'Cluster')
        C = C.Cluster;
    else
        C ={};
    end
    end
end

function CheckSpoolButton(DATA)
    if ~isfigure(DATA.xyfig)
        return;
    end
  spoolbutton = findobj(DATA.xyfig,'Tag','SpoolSpikes');
  if  isempty(spoolbutton)
      return;
  end
if DATA.spooling == 2 || DATA.firsttrial > 1
    set(spoolbutton,'backgroundcolor','r');
else
    set(spoolbutton,'backgroundcolor','w');
end

function cluster = CopyClusters(cluster, ec)
%
% Copy clusters from an expt (ec) to main cluster. Only copy if set.
% that way newly defined clusters on a probe are not wiped out.

for j = 1:size(ec,1)
for k = 1:size(ec,2)
    if isfield(ec{j,k},'params')
        cluster{j,k} = ec{j,k};
    end
end
end


    function h =  DrawCluster(cluster, color)

h =0;
if isempty(cluster) | ~isfield(cluster,'x')
    return;
end
if isfield(cluster,'h') & ishandle(cluster.h)
    delete(cluster.h);
end
tmp.r = [cluster.x(2) cluster.y(2)];
tmp.c = [cluster.x(1) cluster.y(1)];
tmp.xrange = cluster.x(3);
tmp.yrange = cluster.y(3);
tmp.angle = -cluster.angle;
tmp.color = color;
tmp = myellipse(tmp);
h = tmp.lasth;


function SetClusterZ(a,b)
DATA = GetDataFromFig(a);
b = get(findobj(a,'Tag','SetClusterZ'),'value');
DATA.plot.clusterZ = b;
Plot3DClusters(DATA,0);
set(DATA.toplevel,'UserData',DATA);


function DATA = SpkVMenu(a,b, type)
DATA = GetDataFromFig(a);

if type == 1
    a = max(DATA.AllData.Spikes.values(DATA.spklist,:),[],2);
    vm = min([prctile(a,99.5) .* 1.1 5]);
    DATA.plot.SpikeMaxV = vm;
    a = min(DATA.AllData.Spikes.values(DATA.spklist,:),[],2);
    vm = max([prctile(a,0.5) .* 1.1 -5]);
    DATA.plot.SpikeMinV = vm;
    if isempty(b) %called from gui, so fix the ranges
        DATA.plot.autoVrange = 0;
    end
    set(DATA.toplevel,'UserData',DATA);
elseif type ==2
    Spks = GetSpikeStruct(DATA);
    id = find(Spks.codes(:,2) < 8);
    a = FitGaussMeans([DATA.Spikes.cx(DATA.spklist)' DATA.Spikes.cy(DATA.spklist)'],2, 'verbose');
    hold on;
    ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
    plot(a.obj.mu(1,1),a.obj.mu(1,2),'+');
    plot(a.obj.mu(2,1),a.obj.mu(2,2),'+');
    
    if isfield(Spks,'pcs')
        a = FitGaussMeans(Spks.pcs(DATA.spklist,:),2);
        fprintf('Distance for %d PCs %.2f\n',size(Spks.pcs,2),a.mahal);
    end
elseif type ==3  %recalc PCA
    DATA = CheckForPCA(DATA,DATA.spklist, 1);
    set(DATA.toplevel,'UserData',DATA);
end

 function res = FitGaussMeans(X,N, varargin)
    verbose = 0;
    dprime = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'verbose',4)
            verbose = 1;
        end
        j = j+1;
    end
           
    G = gmdistribution.fit(X,N,'Options',statset('MaxIter',1000));
    res.obj = G;
    distance = sqrt(mahal(G,G.mu));
    distance = sqrt(1./(1./distance(1,2)+1./distance(2,1)));
    sigmas = diag(mean(G.Sigma,3)); % sum variances across componets, in each dimenions

    if N == 2
    nsd = diff(G.mu)./sqrt(sigmas)';
    dprime = sqrt(sum(nsd.^2));
    res.dprime = dprime;
    end
    res.mahal = distance;
    if verbose
        fprintf('Distance %.2f (%.2f)\n',distance,dprime);
    end
    

    
function SetTrialRange(a,b, type)
DATA = GetDataFromFig(a);

if type == 2
    DATA.firsttrial = DATA.currenttrial;
elseif type ==3
    DATA.firsttrial = 0;
    DATA.lasttrial = 0;
elseif type ==1
    DATA.lasttrial = DATA.currenttrial;
end

set(DATA.toplevel,'UserData',DATA);
CheckSpoolButton(DATA);
RePlotCellList(DATA);



function SetXYCluster(a,b, type, plot, val,varargin)
DATA = GetDataFromFig(a);
recalcxy = 0;
if type == 1
DATA.plot.clusterX = val;
recalcxy = 1;
elseif type == 2
DATA.plot.clusterY = val;
recalcxy =1;
elseif type == 3
DATA.plot.clusterZ = val;
elseif type == 4
DATA.plot.clusterX = val(1);
if length(val) > 1
    DATA.plot.clusterY = val(2);
else
DATA.plot.clusterY = varargin{1};
end
elseif type == 5
DATA.plot.clusterX = val(1);
if length(val) > 2
    DATA.plot.clusterY = val(2);
    DATA.plot.clusterZ = val(3);
else
DATA.plot.clusterY = varargin{1};
DATA.plot.clusterZ = varargin{2};
end
elseif type == 6
    DATA.clusterErange = [1:31] + (plot-1) * 31;
end
ispk = DATA.spklist;
if (plot == 3 | type == 3 ) & DATA.plot.clusterZ > 0
    Plot3DClusters(DATA,recalcxy);
else
        DATA = CalcClusterVars(DATA, DATA.Expts{DATA.currentexpt}.gui.spks,'force');
        hold off;
    DATA = DrawXYPlot(DATA, ispk); %returned DATA has chnges to clusterX/Yrange
end
set(DATA.toplevel,'UserData',DATA);
if sum(ismember([DATA.plot.clusterY DATA.plot.clusterX DATA.plot.clusterZ],[33 34 37 38]))
    PlotEig(DATA);
end

    
SetGui(DATA);


function PlotEig(DATA)
    if isfield(DATA.AllData,'EigVec')
   GetFigure('EigenVectors');
   plot(DATA.AllData.EigVec(:,end-5:end));
    end
   
function DATA = Plot3DClusters(a,b);
recalc = 0;
    if isfield(a,'AllData');
        DATA = a;
        recalc = b;
    else
    DATA = GetDataFromFig(a);
    end


[xyzfig, isnew] = GetFigure('Cluster3D');
if isnew
    if ~isfield(DATA.plot,'clusterZ')
        DATA.plot.clusterZ = 5;
    end
    z.parentfigtag = DATA.tag.top;
    z.parentfig = DATA.toplevel;
    set(xyzfig,'UserData',z);
    bp = [10 10 100 DATA.plot.ch];
%   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'Position', bp,...
%      'Tag','SetClusterZ','Callback',@SetClusterZ,'value',DATA.plot.clusterZ);
    xyfig = xyzfig;
  hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 1, 3, k});
    end
    hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 2,3, k});
    end
    hm = uimenu(xyfig,'Label','Z','Tag','3DplotMenu');
     uimenu(hm,'Label','None','Callback',{@SetXYCluster, 3, 0});
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 3, 3, k});
    end

  
  DATA.figs.xyz = xyzfig;
    end
Spks = DATA.AllData.Spikes;
ispk = DATA.spklist;
DATA = CheckForPCA(DATA, ispk, 0);
if isempty(DATA.AllData.pcs)
    PCs = [];
else
PCs = DATA.AllData.pcs(ispk,:);
end
classify = 1;
[cz, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterZ, classify,PCs);
DATA.Spikes.cz(ispk) = cz;
if recalc
[cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs);
DATA.Spikes.cx(ispk) = cx;
[cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs);
DATA.Spikes.cy(ispk) = cy;
end
cl = unique(Spks.codes(ispk,2));
hold off;
for j = 1:length(cl)
    id = find(Spks.codes(ispk,2) == cl(j));
plot3(DATA.Spikes.cx(ispk(id)),DATA.Spikes.cy(ispk(id)),cz(id),'.',...
    'color',DATA.spkcolor{cl(j)+1},'markersize',DATA.ptsize);
     hold on;
end
xlabel(['X:' DATA.spkvarnames{DATA.plot.clusterX}]);
ylabel(['Y:' DATA.spkvarnames{DATA.plot.clusterY}]);
zlabel(['Z:' DATA.spkvarnames{DATA.plot.clusterZ}]);

function AddClusterIdButton(a,b)
    DATA = GetDataFromFig(a);
    it = findobj(DATA.xyfig,'Tag','OptimizeDprime');
    fc = findobj(DATA.xyfig,'Tag','ForceClusterid');
    if isempty(fc)
        set(a,'label','Stop cluster ids forcing');
    bp = get(it,'position');
    bp(2) = bp(2)+bp(4);
   
    uicontrol(DATA.xyfig,'style','pop','string','1|2|3|4|5|6|7','Position',bp,'Tag','ForceClusterid',...
        'Callback',@Update);
    else
        set(a,'label','Force cluster ids');
        delete(fc);
    end

    

function DATA = PlaySpikes(DATA, expid)

    
mode = 2;
cw = DATA.plot.cw;
ch = DATA.plot.ch;
rh = ch+10;
[xyfig, isnew] = GetFigure(DATA.tag.clusterxy);
DATA.xyfig = xyfig;
if isnew
        if ~isfield(DATA,'figpos') | isempty(DATA.figpos{2})
        bp = get(DATA.toplevel,'Position');
        fp = get(xyfig,'Position');
       DATA.figpos{2} =  [bp(1)+bp(3) bp(2)+bp(4)-fp(4) fp(3) fp(4)];
       if fp(4) == 0 %can happen with double clicks
           fp(4) = 420;
       end
       if fp(3) == 0
           fp(3) = 420;
       end
       if sum(DATA.figpos{2}([2 4]))+100 > DATA.gui.scrsz(4)
           DATA.figpos{2}(2) = DATA.gui.scrsz(4) - fp(4) - 100;
       end
        set(xyfig,'Position',DATA.figpos{2});
        end
    bp = [5 5 40 20];
    cp = bp;
    uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7|Artifact','Position',bp,'Tag','Clusterid',...
        'Callback',@Update);
    cp(3) = cw * 5;
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'checkbox', 'Callback', @DensityPlot, ...
'String', 'Dens','Tag','Density', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@NextList, 'setfirst'}, ...
'String', 'Set+Next', 'Position', cp,'tag','Set+Next');
    cp(2) = cp(2)+rh;

    
    if length(DATA.probelist) > 1
        cp(3) = cw*1.5;
        uicontrol(xyfig,'style','pushbutton','string','+','Position',cp, 'Callback',@AddOneCellToList);

        cp(1) = cp(1) + cp(3);
        cp(3) = cw*3;
        uicontrol(xyfig,'style','pop','string','1|2|3|4|5|6|7|8|9|10','Position',cp,'Tag','AddOneCellToList',...
            'Callback',@AddOneCellToList);
        cp(2) = cp(2)+rh;
        cp(1) = 5;
        cp(3) = 40;
    end



    uicontrol(gcf,'Style', 'pushbutton', 'Callback', @NextList, ...
'String', 'Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@NextList, 'clearfirst'}, ...
'String', 'Clr+Next', 'Position', cp);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','spool','Position',cp,'Tag','SpoolSpikes',...
        'Callback', @SpoolSpikes);
    cp(2) = cp(2)+rh;
    uicontrol(gcf,'style','pushbutton','string','Optim','Position',cp,'Tag','OptimizeDprime',...
        'Callback', @OptimizeDprimeHit);
    cp(2) = cp(2)+rh;   
    
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pop','string','Not Set|Nothing|MU-|MU+|Poor|OK|Good|V Good|Excellent|FromDprime','Position',bp,'Tag','ClusterQuality',...
        'Callback',@Update);

    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Set','Position',bp,'Callback', @SetExptClusters);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw;
    uicontrol(xyfig,'style','CheckBox','Position',bp,'Tag','ClusterIsSet');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*3;
    uicontrol(xyfig,'style','pushbutton','string','Del','Position',bp,'Callback', @DelClusterButton);
    bp(1) = bp(1) + bp(3) + 10;
    uicontrol(xyfig,'style','pushbutton','string','Clr','Position',bp,'Callback', @ClrSpkWin);
    bp(1) = bp(1) + bp(3) + 10;
    bp(3)=cw*5;
    uicontrol(xyfig,'style','text','string','Max: X','Position',bp);
    bp(1) = bp(1) + bp(3);
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterXrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterXmax');

    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Y','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterYrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterYmax');
    bp(1) = bp(1) + bp(3) + 10;
    bp(3) = cw*1;
    uicontrol(xyfig,'style','text','string','Z','Position',bp);
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*4;
    uicontrol(xyfig,'style','edit','string',sprintf('%.2f',DATA.plot.clusterZrange(2)),'Position',bp,...
    'Callback', @RescaleClusterPlot,'Tag','ClusterZmax');
    bp(1) = bp(1) + bp(3);
    bp(3) = cw*5;
    uicontrol(xyfig,'style','CheckBox','string','auto','Position',bp,'Tag','AutoScale',...
        'value',(DATA.plot.autoscale > 0),'Callback',@Update);
    bp(1) = bp(1) + bp(3) + 10;
    hm = uimenu(xyfig,'Label','Plot3D','Tag','3DplotMenu');
     uimenu(hm,'Label','None','Callback',{@SetXYCluster, 3, 2, 0});
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 3, 2,k});
    end
    hm = uimenu(xyfig,'Label','X','Tag','XClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 1, 2,k});
    end
    hm = uimenu(xyfig,'Label','Y','Tag','YClusterMenu');
    for j = 1:length(DATA.spkvarorder)
        k = DATA.spkvarorder(j);
        uimenu(hm,'Label',DATA.spkvarnames{k},'Callback',{@SetXYCluster, 2, 2, k});
    end
    hm = uimenu(xyfig,'Label','Ops','Tag','XYAddMenu');
     uimenu(hm,'Label','Force Cluster Ids','Callback',{@AddClusterIdButton});
     a = strmatch('Energy',DATA.spkvarnames,'exact');
     b = strmatch('Var/Energy',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','Energy-Var/E','Callback',{@SetXYCluster, 4, 2, a(1), b(1)});
     b = strmatch('sqrt(Var/Energy)',DATA.spkvarnames,'exact');
     a = strmatch('sqrt(Energy)',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','sqrt(Energy)-sqrt(Var/E)','Callback',{@SetXYCluster, 4, 2, a(1), b(1)});
     b = strmatch('Var/Energy',DATA.spkvarnames,'exact');
     a = strmatch('sqrt(Energy)',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','sqrt(Energy)-Var/E','Callback',{@SetXYCluster, 4, 2, a(1), b(1)});
     b = strmatch('PCA1',DATA.spkvarnames,'exact');
     a = strmatch('PCA2',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','PCA1 - 2','Callback',{@SetXYCluster, 4, 2, a(1), b(1)});
     b = strmatch('PCA2',DATA.spkvarnames,'exact');
     a = strmatch('PCA3',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','PCA2 - 3','Callback',{@SetXYCluster, 4, 2, a(1), b(1)});
     if DATA.subprobes > 1
     b = strmatch('Energy 1',DATA.spkvarnames,'exact');
     a = strmatch('Energy 2',DATA.spkvarnames,'exact');
     c = strmatch('Energy 3',DATA.spkvarnames,'exact');
     uimenu(hm,'Label','Energy1-2-3','Callback',{@SetXYCluster, 5, 3, a(1), b(1), c(1)});
     uimenu(hm,'Label','Pts trode 1','Callback',{@SetXYCluster, 6,1});
     uimenu(hm,'Label','Pts trode 2','Callback',{@SetXYCluster, 6,2});
     uimenu(hm,'Label','Pts trode 3','Callback',{@SetXYCluster, 6,3});
     uimenu(hm,'Label','Pts trode 4','Callback',{@SetXYCluster, 6,4});
     end     
     uimenu(hm,'Label','recalc PCA','Callback',{@SpkVMenu, 3});
     uimenu(hm,'Label','Gmix Distance','Callback',{@SpkVMenu, 2});
    tmpdat.parentfigtag = DATA.tag.top;
    set(xyfig,'UserData',tmpdat);
end
hold off;
CheckSpoolButton(DATA);
%if isempty(findobj('Tag',DATA.tag.spikev
 %   sfig = figure('Renderer','painters','Tag','DATA.tag.spikev');
[sfig, isnew] = GetFigure(DATA.tag.spikev);
if isnew 
   if ~isfield(DATA,'figpos') | length(DATA.figpos) < 3 | isempty(DATA.figpos{3})
        bp = get(DATA.toplevel,'Position');
        fp = get(sfig,'Position');
       DATA.figpos{3} =  [bp(1) bp(2)-fp(4)-80 fp(3) fp(4)];
        set(sfig,'Position',DATA.figpos{3});
    end

    x = 10;
    c = 10;
    bp = [10 10 40 20];
    uicontrol(sfig,'style','pushbutton','string','>>','Position',bp,'Tag','NextTrial',...
        'Callback', @PlayNextTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','<<','Position',bp,'Tag','LastTrial',...
        'Callback', @PlayLastTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','X','Position',bp,'Tag','CutTrial',...
        'Callback', @CutTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pop','string','1|2|3','Position',bp,'Tag','ChooseTrial',...
        'Callback', @SelectTrial);
    bp(1) = bp(1)+bp(3);
    uicontrol(sfig,'style','pushbutton','string','+','Position',bp,'Tag','NextSpike',...
        'Callback', @PlayNextSpike);
    bp(1) = bp(1)+bp(3);
    bp(3)=c*6;
    uicontrol(sfig,'style','pushbutton','string','spool','Position',bp,'Tag','SpoolSpikes',...
        'Callback', @SpoolSpikes);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
    'String', 'dVdt', 'Tag', 'dVdt', 'Position', bp,'value',DATA.plot.dvdt,...
    'Callback',@Update);
    bp(1) = bp(1) + bp(3);
    
    bp(3)=c*9;
    uicontrol(gcf,'Style', 'checkbox',...
    'String', 'Stop', 'Tag', 'StopSpool', 'Position', bp,'value',0,...
    'Callback',@Update);
    tmpdat.parentfigtag = DATA.tag.top;
    set(sfig,'UserData',tmpdat);
    set(sfig, 'WindowScrollWheelFcn',@ScrollTrial);

      hm = uimenu(sfig,'Label','Set','Tag','SpikeVMenu');
      uimenu(hm,'Label','All Trials','Callback',{@SetTrialRange, 3});
      uimenu(hm,'Label','End Range','Callback',{@SetTrialRange, 1});
      uimenu(hm,'Label','Range->','Callback',{@SetTrialRange, 2});
      uimenu(hm,'Label','Vrange','Callback',{@SpkVMenu, 1});


    DATA.hline = 0;
end
    it = findobj(sfig,'Tag','ChooseTrial');
    set(it,'string',sprintf('%d|',[DATA.Expts{expid}.Trials.Trial]),'value',1);
    DATA.svfig = sfig;
if DATA.spooling ~= 2 %2 == spool from current trial
    DATA.currenttrial = 0;
end
DATA.ISIpair = 0;
DATA.currentexpt = expid;
if DATA.test.fastplot
    set(DATA.svfig,'renderer','painters');
    ax = findobj(DATA.svfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    set(ax,'Xlim',[1 46],'Ylim',[-5 5]);
end
ax = findobj(DATA.svfig,'type','ax');



if DATA.state.uselfp
    DATA.state.lfig = GetFigure('LFP');
    hold off;
    DATA.Expts{expid} = LoadSpike2LFP(DATA.Expts{expid});
end
colors = mycolors;
if DATA.probe == 100
    Spks = DATA.AllData.UstimV;
elseif isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
    Spks = DATA.AllData.Spikes;
end
tstart = now;
if isfield(Spks,'values')
nsamples = size(Spks.values,2).* size(Spks.values,3);
else
    nsamples = 0;
end
if size(Spks.codes,2) > 1
nclusters = 1+max(Spks.codes(:,2));
else
nclusters = 1;
end
set(0,'CurrentFigure',xyfig);

if isfield(DATA.plot,'useprobe')
    xprobes = find(DATA.plot.useprobe);
    xprobes = setdiff(xprobes,DATA.probe);
else
    xprobes = [];
end
DATA.xprobes = xprobes;
probes = [DATA.probe xprobes];
nprobes = length(probes);



reclassify = 0;
cid = findobj('Tag','ClusterIsSet');
if DATA.spooling
    DATA.state.recut = 1;
elseif ClusterIsSet(DATA.Expts{expid},DATA.probe) & DATA.state.recut
    DATA.cluster = CopyClusters(DATA.cluster,DATA.Expts{expid}.Cluster);
    DATA.state.recut = 1;
    if ismember(DATA.Expts{expid}.gui.clustertype, [0 2 3]) %% The online cut or autocut
        set(cid,'value',0);
    else
        set(cid,'value',1);
    end
    if ~isfield(DATA.Expts{expid}.gui,'classified') | DATA.Expts{expid}.gui.classified ~= 1
     DATA = SetExptSpikes(DATA,expid,0);
    end
    if length(xprobes) & isfield(DATA,'AllSpikes')
        for j = 1:length(probes)
            if ~isfield(DATA.AllSpikes{probes(j)},'spklist')
                times(1) = DATA.Expts{expid}.Trials(1).Start(1)-500;
                times(2) = DATA.Expts{expid}.Trials(end).End(end)+500;
                DATA.AllSpikes{probes(j)}.spklist = FindSpikes(DATA, times, probes(j),[]);
            end
            DATA = SetSpkCodes(DATA, DATA.AllSpikes{probes(j)}.spklist, probes(j),0);
        end
    end
    
elseif ~isfield(DATA,'cluster')
    DATA.cluster = {};
elseif DATA.state.recut %No cluster yet defined
   DATA.state.recut = 2;
% when inheriting a cluster, don't inherit any spike ranges
  nclusters = 1+max(Spks.codes(:,2));
  p = DATA.probe;
  for j = 1:nclusters
      if j <= size(DATA.cluster,1) & size(DATA.cluster,2) >= p & ~isempty(DATA.cluster{j,p})
      if isfield(DATA.Expts{expid}.gui,'spkrange') & ~isempty(DATA.Expts{expid}.gui.spkrange)
      DATA.cluster{j,p}.firstspk = DATA.Expts{expid}.gui.spkrange(1);
      DATA.cluster{j,p}.lastspk = DATA.Expts{expid}.gui.spkrange(2);
      else
      DATA.cluster{j,p}.firstspk = 1;
      DATA.cluster{j,p}.lastspk = length(Spks.times);
      end
      end
  end
    set(cid,'value',0);
end
nclusters = double(nclusters);
if DATA.plot.synccluster == 0
    ClrSpkWin(DATA);
    DATA = DrawClusters(DATA,DATA.cluster,0);
    if isfield(DATA.Expts{expid},'OnlineCluster') & DATA.state.showonlineclusters
        DATA = DrawClusters(DATA,DATA.Expts{expid}.OnlineCluster,0);
    end
else
    hold off;
    plot(0,0,'+');
    hold on;
    DATA = DrawClusters(DATA,DATA.cluster,0);
end
if DATA.plot.autoscale == 0
  set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange); 
end 

set(0,'CurrentFigure',sfig);
hold off;

    if DATA.spikelist == -1
       spikelist = [0 1 2 3 4];
    else
       spikelist = DATA.spikelist;
    end



if isfield(DATA,'AllSpikes')
for k = 1:length(probes)
        nclusters = max([nclusters max(DATA.AllSpikes{probes(k)}.codes(:,2))]);
end
end
n = CountClusters(DATA.cluster);
nclusters = max([n nclusters]);
x = 32; %should calculate real size..
for k = 0:length(xprobes)
    ids{k+1} = FindSpikes(DATA, DATA.Expts{expid}.Header.trange,probes(k+1),[]);
    for j = 1:nclusters+1
        DATA.svh(j+k*(nclusters+1)) = plot([1:nsamples],[1:nsamples] * 0,'color',DATA.spkcolor{j});
        hold on;
    end
    if isfield(DATA,'AllSpikes')
        %nclusters needs to be max # used across all probes
        if length(DATA.plot.voffsets) > k & DATA.plot.prettyfigs
        h = text(25,DATA.plot.voffsets(k+1)+1,sprintf('Probe %d',probes(k+1)));
        set(h,'FontWeight','bold');
        else
        text(x,(k-(nprobes-1)/2) * DATA.plot.SpikeVsep,num2str(probes(k+1)));
        end
        ncut = sum(DATA.AllSpikes{probes(k+1)}.codes(:,2) > 0);
        if ncut
            id = find(ismember(DATA.AllSpikes{probes(k+1)}.codes(ids{k+1},2),spikelist));
%            ids{k+1} = ids{k+1}(id);
        end
        if ~isfield(DATA.AllSpikes{probes(k+1)},'dVdt')
            DATA.AllSpikes{probes(k+1)} = CleanSpikes(DATA.AllSpikes{probes(k+1)});
        end
    end
end
DATA.svhn = nclusters+1;

expspks = ids{1}; 
if ~isfield(DATA.Expts{expid}.gui,'spks') %haven't loaded this expt yet
    times = [DATA.Expts{expid}.Trials(1).Start(1) DATA.Expts{expid}.Trials(end).End(end)];
    if isfield(DATA,'AllSpikes')
    nspk = sum(DATA.AllSpikes{DATA.probe}.times > times(1) & DATA.AllSpikes{DATA.probe}.times < times(2)); 
    else
    nspk = sum(DATA.AllData.Spikes.times > times(1) & DATA.AllData.Spikes.times < times(2)); 
    end
    if nspk < 10000
        DATA.ptsize = 6;
    elseif nspk < 2000
        DATA.ptsize = 10;
    else
        DATA.ptsize = 4;
    end
    expspks = ids{1}; 
else
%    expspks = DATA.Expts{expid}.gui.spks; % no good if > 1 probe loaded
end

if DATA.plot.autoVrange
    DATA.spklist = expspks;
    DATA = SpkVMenu(DATA,0,1);
end
DATA.ptsize = CheckPtSize(DATA, expspks);

%FindSpikes might list some spikes just past the list recorded in
%Expts.gui.spkrange.
if isfield(DATA,'Spikes') & isfield(DATA.Spikes,'cx')
if max(ids{1}) > length(DATA.Spikes.cx)
    ispk = length(DATA.Spikes.cx):max(ids{1});
    DATA = CalcClusterVars(DATA,  ispk);
end
end

if length(ids{1}) & isfield(DATA.AllData.Spikes,'codes') & size(DATA.AllData.Spikes.codes,1) >= ids{1}(end)
id = find(DATA.AllData.Spikes.codes(ids{1},2) > nclusters);
if length(id)
    DATA.AllData.Spikes.codes(ids{1}(id),2) = 0;
end
end

DATA = CheckForPCA(DATA, expspks, 0);
        
        
DATA.vstep = DATA.plot.SpikeMaxV * 0.8;
    if DATA.plot.SpikeVsep > 0
        DATA.vstep = DATA.plot.SpikeVsep;
    end
if DATA.plot.dvdt == 2
    set(gca,'Xlim',[-5 5],'Ylim',[-5 5]);
elseif length(xprobes)
    x = (length(xprobes)-1) .* DATA.plot.SpikeMaxV/2;
    if length(xprobes) > 3
        x = (length(xprobes)-1) .* DATA.plot.SpikeVsep;
    end
    x = (length(xprobes)) .* DATA.plot.SpikeVsep;
    if length(xprobes) 
    set(gca,'Xlim',[1 nsamples+1],'Ylim',[-(DATA.plot.SpikeMaxV+x/2) DATA.plot.SpikeMaxV+x/2]);
    end
elseif nsamples > 1        
    set(gca,'Xlim',[1 nsamples],'Ylim',[-DATA.plot.SpikeMaxV DATA.plot.SpikeMaxV]);
end
DATA.nclusters = CountClusters(DATA.cluster);
nt = 1;
allspks = [];
firstspk = 0;
hline = 0;
if ~isfield(DATA,hline)
    DATA.hline = 0;
end

DATA.playingspk = 1;
DATA.nclusters = CountClusters(DATA.cluster);
if DATA.nclusters > nclusters
    nclusters = DATA.nclusters;
end
DATA.svhn = nclusters+1;

set(DATA.toplevel,'UserData',DATA);
if isempty(DATA.currentcluster) || length(DATA.cluster) < DATA.currentcluster 
    DATA.currentcluster = 1;
end

Aargs = {}; timemode = 0;
trialdur = DATA.Expts{expid}.Trials(1).End(end) - DATA.Expts{expid}.Trials(1).Start(1);
timemode = DATA.plot.showwave;

if isfield(DATA,'AllSpikes') & length(probes) > 1 & DATA.plot.syncoverlay == 0
    Aargs = {Aargs{:} 'timemode'};
    timemode = 1;
end



if timemode
        if  ~isfield(DATA,'timefig')
        DATA.timefig = GetFigure('SpikeTime');
    end
    set(0,'CurrentFigure',DATA.timefig);
    set(gca,'Xlim',[0 trialdur/10]);
    for k = 0:length(xprobes)
        for j = 1:nclusters+1
            DATA.tvh(j+k*(nclusters+1)) = plot([1:46],[1:46] * 0,'color',DATA.spkcolor{j});
            hold on;
        end
    end
end


nprobes = length(probes);
step = DATA.plot.SpikeMaxV;
DATA.syncprobes = [];
if length(probes) > 1 && isfield(DATA,'AllSpikes')
    if DATA.plot.synccluster > 0
        DATA.Spikes.cx(1:end) = 0;
        DATA.Spikes.cy(1:end) = 0;
        DATA.syncprobes = probes(1:2);
    else
    end
%syncsign = -1 show only sychronous spikes with negetive trigger
%syncsign = -1 show only sychronous spikes with negetive trigger
%syncsign = 11 show only sychronous spikes with positive trigger

    
% even if plotting non sync spikes, want to know which are synced
%    if DATA.syncsign < 2
        dt = 2;


        %find synchronous spikes

        for j = 2:length(probes)
            if isfield(DATA,'sids') && sum(ismember(DATA.sids{1},DATA.AllSpikes{probes(1)}.times(ids{1}))) > 100
                aid = DATA.sids{1};
                bid = DATA.sids{2};
            else
            [aid,bid] = FindSync(DATA.AllSpikes{probes(1)}.times(ids{1}),...
                DATA.AllSpikes{probes(j)}.times(ids{j}),dt);
            aid = ids{1}(aid);
            bid = ids{j}(bid);
            end
            dcs(1,1:length(aid)) = mean(DATA.AllSpikes{probes(1)}.values(aid,:),2);
            dcs(j,1:length(bid)) = mean(DATA.AllSpikes{probes(j)}.values(bid,:),2);
            if ismember(DATA.syncsign,[-1 3])
                id = find(DATA.AllSpikes{probes(1)}.values(aid,9) < 0 & ...
                    DATA.AllSpikes{probes(j)}.values(bid,9) < 0);
                aid = aid(id);
                bid = bid(id);
            elseif ismember(DATA.syncsign,[1 4])
                id = find(DATA.AllSpikes{probes(1)}.values(aid,9) > 0 & ...
                    DATA.AllSpikes{probes(j)}.values(bid,9) > 0);
                aid = aid(id);
                bid = bid(id);
            end
            DATA.sids{1} = aid;
            DATA.sids{j} = bid;
            aids{j} = aid;
            if ismember(DATA.plot.synccluster,[3 4 5 6]) % need PCA
                [E, V, DATA.pca] = SpikePCA(DATA.AllSpikes,[probes(1) probes(j)],{aid bid} );
            end
        end
        if length(probes) == 3 %intersection of all three
            [DATA.sids{1}, bi, ci] = intersect(aids{2},aids{3});
            DATA.sids{2} = DATA.sids{2}(bi);
            DATA.sids{3} = DATA.sids{3}(ci);        
        end

        if DATA.plot.prettyfigs ~= 1
        o = (length(probes)-1)/2;
        for j = 1:length(probes)
            if ~isempty(DATA.sids{j})
            mspk(j,:) = mean(DATA.AllSpikes{probes(j)}.values(DATA.sids{j},:),1);
            plot(mspk(j,:)+((j-1-o) * DATA.vstep),':');
            end
        end
        end

end
start = max([DATA.currenttrial 1]);
last = length(DATA.Expts{expid}.Trials);
% negative Trial numbers indicate manual exclusion
uset = start-1+find([DATA.Expts{expid}.Trials(start:last).Trial] > 0);
nsync = 0;nspks = 0;
tc = 1;
syncspikes = [];

while tc <= length(uset)  %not for loop, so can change tc inside loop
    trial = DATA.Expts{expid}.Trials(uset(tc)).Trial;
    stop = get(findobj(DATA.svfig,'Tag','StopSpool'),'value');
    if stop
        set(findobj(DATA.svfig,'Tag','StopSpool'),'value',0);
        if 0
            DATA.playingspk = 0;
            return;
        else
            tc = length(uset);
            trial = DATA.Expts{expid}.Trials(uset(nt)).Trial;
        end
    end
    if DATA.state.uselfp
        GetFigure('TrialLFP');
    PlotLFPRaw(DATA.state,DATA.Expts{expid}.Trials(uset(nt)),DATA.Expts{expid}.Header.LFPsamplerate);
    end
    if timemode
    set(0,'CurrentFigure',DATA.timefig);
    else
    set(0,'CurrentFigure',sfig);
    end
    hold off;
    itrial = find(DATA.AllData.Trialids == trial);
    if mode == 1
        DATA = PlotTrialSpikes(DATA,itrial,colors, clusters);
    elseif mode == 2
        times(1) = DATA.Expts{expid}.Trials(uset(nt)).Start(1)-DATA.state.preperiod;
        times(2) = DATA.Expts{expid}.Trials(uset(nt)).End(end)+DATA.state.preperiod;
        Trial = DATA.Expts{expid}.Trials(uset(nt));
        if isfield(Trial,'uStim') && Trial.uStim > 0
            ustim = 1;
        else
            ustim = 0;
        end
        Trial.ed = GetEval(DATA.Expts{expid},'ed',DATA.currenttrial);
        
        if DATA.syncsign < 2 & length(probes) > 1
            [spks, sspks, cx, cy] = PlotTrialSyncSpikes(DATA, times, [DATA.probe xprobes], colors,'Trial',Trial);
            if DATA.plot.synccluster > 0
                syncspikes = cat(1,syncspikes,sspks);;
            DATA.Spikes.cx(sspks(:,1)) = cx;
            DATA.Spikes.cy(sspks(:,1)) = cy;
            DATA.AllSpikes{DATA.probe}.cx(sspks(:,1)) = cx;
            DATA.AllSpikes{DATA.probe}.cy(sspks(:,1)) = cy;
            DATA.AllSpikes{xprobes(1)}.cx(sspks(:,2)) = cx;
            DATA.AllSpikes{xprobes(1)}.cy(sspks(:,2)) = cy;
            end
            nspks = nspks+length(spks);
            nsync = nsync+size(sspks,2);
            if DATA.plot.synccluster == 10
                    [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probes', [probes(1:2)]);
            end

        else

            if length(probes) > 1
                if DATA.plot.syncoverlay
                for j= 1:length(xprobes)
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1));
                end
                [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)]);
                end
                if timemode
                for j= 1:length(xprobes)
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe',[xprobes(j) 1+j length(probes)],'lineoff',j*(nclusters+1),'timemode');
                end
                [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],'timemode');
                end
            else %single electrode
                if ustim == 0 || (probes(1) > 99 && ustim == 1)
                    [DATA, spks] = APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:});
                else
                    spks = [];
                end
                if timemode
                    APlotTrialSpikes(DATA,times,colors, nclusters,1,'Trial',Trial,'Probe', [probes(1) 1 length(probes)],Aargs{:},'timemode');
                end
            end
        end
        allspks = [allspks spks'];
        nt = nt+1;
    end
    hold on;
    tc = tc+1;
end
if DATA.plot.synccluster > 0
    DATA.syncspikes = syncspikes;
end
if DATA.plot.showsync
    fprintf('%d/%d (%.3f) spikes Synchronized\n',nsync,nspks,nsync./nspks);
end
            

if DATA.probe == 100;
    DATA.Expts{expid}.MeanPulse = mean(Spks.values(allspks,:));
    plot(DATA.Expts{expid}.MeanPulse);
end
if DATA.state.uselfp  % look for calibration spikes

 GetFigure('LFP');
 hold off;
 CalcLFPPulse(DATA.Expts{expid},DATA.AllData,'plot');
 GetFigure('SpikeV');
end
if isempty(allspks)
    DATA.spkrange(1) = 1;
    DATA.spkrange(2) = 1;
else
    DATA.spkrange(1) = min(allspks);
    DATA.spkrange(2) = max(allspks);
    DATA.spklist = allspks;
end
DATA.Expts{expid}.gui.spkrange = DATA.spkrange;
DATA.Expts{expid}.gui.spks = allspks;
DATA.Expts{expid}.gui.s2clusters = 1+max(Spks.codes(allspks,1));
DATA.s2clusters = DATA.Expts{expid}.gui.s2clusters;
DATA = FinishXYPlot(DATA);
set(xyfig, 'KeyPressFcn',@KeyPressed);
set(xyfig, 'KeyReleaseFcn',@KeyReleased);
set(xyfig, 'WindowButtonDownFcn',@ButtonPressed);
set(xyfig, 'WindowButtonMotionFcn',@ButtonDragged);
set(xyfig, 'WindowButtonUpFcn',@ButtonReleased);
set(xyfig, 'WindowScrollWheelFcn',@ScrollWheel);
axis('manual');
ClearMouse;
DATA.densityplot = 0;
xrange = get(gca,'Xlim');
yrange = get(gca,'Ylim');


if DATA.test.fastplot
    set(DATA.xyfig,'renderer','painters');
    ax = findobj(DATA.xyfig,'type','ax');
    set(ax,'DrawMode','fast','XTickMode','manual','YTickMode','manual');
    %set(,'DrawMode','fast');
end
if DATA.state.fixrange
    if length(allspks) > 1000
        prc = 99;
    else
        prc = 95;
    end
    xm = prctile(DATA.Spikes.energy(allspks),prc);
    ym = prctile(DATA.Spikes.vw(allspks),prc);
    if yrange(2) > 2 * ym
        set(gca,'Ylim',[0 ym * 1.5]);
    end
    if xrange(2) > 2 * xm
        set(gca,'Xlim',[0 xm * 1.5]);
    end
end
mousept.xrange = diff(xrange);
mousept.yrange = diff(yrange);
fprintf('Spool Took %.2f\n',(now-tstart) * 24 * 60 * 60);
if isfield(DATA,'AllClusters')
    set(0,'currentfig',DATA.xyfig);
    plot(DATA.Spikes.cx(allspks),DATA.Spikes.cy(allspks),'.','markersize',DATA.ptsize);
    DATA.allexp = DATA.currentexpt;
    DATA = PlotAllProbeXY(DATA,0);
end
if isfield(DATA,'AllSpikes') && nprobes == 2 && DATA.plot.xcorr
    GetFigure('Xcorr');
    xc = CalcXcorr(DATA,DATA.currentexpt,probes(1),probes(2));
end

if DATA.plot.voltxy == 4
    PlotTrodeXcorr(DATA,0);
end

DATA.playingspk = 0;
SetGui(DATA);


    function PlotPCAs(DATA, eid)
    for e = 1:length(eid)
        [DATA, allspks] = SetExptSpikes(DATA,eid(e),'setrange');
    end
    
    GetFigure('PCA');
    pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;];
    pw=0.33;
    ph = 0.5;
    cid = unique(DATA.AllData.Spikes.codes(allspks,2));
    for j = 1:size(pairs,1);
        x = mod(j-1,3) * 0.33;
        y = floor((j-1)/3) * 0.5;
        subplot('Position' ,[x y pw ph]);
        hold off;
        for k = 1:length(cid)
            id = find(DATA.AllData.Spikes.codes(allspks,2) == cid(k));
            plot(DATA.AllData.pcs(allspks(id),pairs(j,1)),DATA.AllData.pcs(allspks(id),pairs(j,2)),'.','markersize',1,'color',DATA.spkcolor{k});
            hold on;
        end
        set(gca,'Xtick',[],'YTick',[]);
    end

function OUT = CalcTrodeCorrs(DATA,eid)
    cell = 0;
    mu=0;
    clabel = {};
    spikes = [];
    GetFigure('MeanSpike');
    colors = mycolors;
    hold off;
    GetFigure('MeanMU');
    subplot(2,1,1);
    hold off;
    
    for p = 1:length(DATA.probelist)
        DATA = SetProbe(DATA,DATA.probelist(p));
        eid = eid(find(eid <= length(DATA.Expts)));
        for e = 1:length(eid)
        [DATA, allspks] = SetExptSpikes(DATA,eid(e),'setrange');
        for j = 1:length(allspks);
            xc = corrcoef(reshape(DATA.AllData.Spikes.values(allspks(j),:),32,4));
            xcs(j,:) = [xc(1,2) xc(2,3) xc(3,4) xc(1,3) xc(2,4) xc(1,4)];
        end
        tx(p,e,:) = prctile(xcs,50);
        ds(p,e,:) = [1 1 1 2 2 3];
        cid = unique(DATA.AllData.Spikes.codes(allspks,2));
        for j = 2:length(cid)
            cell = cell+1;
            id = find(DATA.AllData.Spikes.codes(allspks,2) ==cid(j));
            mspk = mean(DATA.AllData.Spikes.values(allspks(id),:));
            spikes(cell).v = mspk;
            spikes(cell).probe = p;
            spikes(cell).cluster = cid(j);
            GetFigure('MeanSpike');
            plot(mspk,'color',colors{cell});
            hold on;
            mspk = reshape(mspk,32,4);
            [a,c] = max(var(mspk));
            amp = mspk' * mspk(:,c);
            amp = amp./amp(c);
            cd = abs([1:4] - c);
            cds(cell,:) = cd(cd > 0);
            ccs(cell,:) = amp(cd > 0);
            clabel{cell} = sprintf('P%dC%d',p,cid(j));
        end
        ts = [9 41 73 105];
            v = prctile(DATA.AllData.Spikes.values(allspks,ts),5);
        for j = 1:length(ts)
            mu = mu+1;
            id = find(DATA.AllData.Spikes.values(allspks,ts(j)) < v(j));
            mspk = mean(DATA.AllData.Spikes.values(allspks(id),:));
            mus(mu).v = mspk;
            mus(mu).probe = p;
            mus(mu).subprobe = j;
            GetFigure('MeanMU');
            subplot(2,1,1);
            plot(mspk,'color',colors{j});
            hold on;
            mspk = reshape(mspk,32,4);
            amp = mspk' * mspk(:,j);
            amp = amp./amp(j);
            cd = abs([1:4] - j);
            mus(mu).amp = amp(cd > 0);
            mus(mu).d = cd(cd > 0);
            muds(mu,:) = cd(cd > 0);
            mucs(mu,:) = amp(cd > 0);
        end
        end
    end
    GetFigure('MeanSpike');
    legend(clabel);
    GetFigure('MeanMU');
    subplot(2,1,2);
    hold off;
    for j = 1:4
    id = find([mus.subprobe] == j);
    plot(mean(cat(1,mus(id).v)),'linewidth',2,'color',colors{j});
    hold on;
    end
    GetFigure('TrodeCorr');
    hold off;
    plot(ds(:),tx(:),'o');
    hold on;
    if exist('cds','var')
    plot(cds(:),ccs(:),'ro');
    else 
        cds = [];
        ccs = [];
    end
    dval = unique(ds);
    plot(muds(:),mucs(:),'go');
    for j = 1:length(dval)
        id = find(cds == dval(j));
        cmean(j) = mean(ccs(id));
        id = find(ds == dval(j));
        dmean(j) = mean(tx(id));
        id = find(muds == dval(j));
        mumean(j) = mean(mucs(id));
    end
    plot(dval,dmean,'-');
    plot(dval,cmean,'r-');
    plot(dval,mumean,'g-');
    
    DATA.TrodeMean.spikes = spikes;
    DATA.TrodeMean.mu = mus;
    OUT = DATA;

function PlotTrodeXcorr(a,b)
    DATA = GetDataFromFig(a);
    allspks = DATA.spklist;
        for j = 1:length(allspks);
        xc = corrcoef(reshape(DATA.AllData.Spikes.values(allspks(j),:),32,4));
        xcs(j,:) = [xc(1,2) xc(2,3) xc(3,4)];
    end
    GetFigure('SpkieCorr');
    hold off;
    plot(xcs(:,2),xcs(:,3),'.','markersize',DATA.ptsize);
    hold on; plot(median(xcs(:,2)),median(xcs(:,3)),'r+');

    
function SetSetBox(DATA, ei)
        
cid = findobj('Tag','ClusterIsSet');
if ismember(DATA.Expts{ei}.gui.clustertype, [0 2 3]) %% The online cut or autocut
    set(cid,'value',0);
else
    set(cid,'value',1);
end

function DATA = GetDataFromFig(a)
% for functions that are called with DATA directly, not from a gui, need to
%retutn first argument
    if isstruct(a) & isfield(a,'toplevel') 
        DATA = a;
        return;
    end
    DATA = get(a,'UserData');

    b = a;
    while isempty(DATA) & ~isempty(b)
        b = get(b,'parent');
        DATA = get(b,'UserData');
    end
    if isfield(DATA,'parentfigtag')
        DATA = get(findobj('Tag',DATA.parentfigtag),'UserData');
    end
    
function RescaleClusterPlot(a,b)
    DATA = GetDataFromFig(a);
it = findobj(DATA.xyfig,'Tag','ClusterXmax');
if it
    DATA.plot.clusterXrange(2) = str2num(get(it,'string'));
end

it = findobj(DATA.xyfig,'Tag','ClusterYmax');
if it
    DATA.plot.clusterYrange(2) = str2num(get(it,'string'));
end
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
if it
    DATA.plot.clusterZrange(2) = str2num(get(it,'string'));
end
ax = findobj(DATA.xyfig,'type','ax');

%
%if the user has typed in a value, surely they want to manually scales
DATA.plot.autoscale = 0;
set(findobj(DATA.xyfig,'Tag','AutoScale'),'value',DATA.plot.autoscale);
if DATA.plot.autoscale == 0
    set(ax,'Xlim',DATA.plot.clusterXrange);
    set(ax,'Ylim',DATA.plot.clusterYrange);
end

    
    
if DATA.densityplot
    caxis([0 DATA.plot.clusterZrange(2)]);
end
set(DATA.toplevel,'UserData',DATA);

function SpoolSpikes(a,b)
DATA = GetDataFromFig(a);
DATA.spooling = 1;
if isfigure(a) & Bstrmatch(get(a,'Tag'),'SpoolSpikes') & DATA.currenttrial > 1
    DATA.spooling = 2; %spool from current spike to end
end
DATA = PlaySpikes(DATA,DATA.currentexpt);
set(DATA.toplevel,'UserData',DATA);

function NextList(a,b, varargin)
%DATA = combine('getstate');

j = 1;
saved = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'clearfirst',6)
        DeleteCluster(0,a,'nodraw');
        SetExptClusters(a,b,'nosave');
    elseif strncmpi(varargin{j},'setfirst',6)
        saved = SetExptClusters(a,b,'nosave');
        if saved == 0
        set(a,'backgroundcolor','r');  %as its not saved
        end
    end
    j = j+1;
end

DATA.savedclusters = saved;
DATA = GetDataFromFig(a);
playspk = get(findobj(DATA.toplevel,'Tag','ShowSpikes'),'value');

if DATA.playingspk
    set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
    return;
end
strs = get(DATA.elst,'string');
val = get(DATA.elst,'value');
if val < length(strs)
    set(DATA.elst,'value',val+1)
    if isfield(DATA,'AllClusters') & playspk == 0
        combine('setexpt',DATA);
        if 0 % need to use setexp to clear cluster.touched, etc
        id = get(DATA.elst,'value');
        DATA.currentexpt = DATA.expid(id(1));
        DATA.allexp = DATA.currentexpt;
        PlotAllProbeXY(DATA);
        end
    else
        combine('setexpt',DATA);
    end
end


function ClearMouse()
global mousept;

mousept.start = [];
mousept.lasth = [];
mousept.angle = 0;
mousept.r = [1 1];
mousept.down = 0;
mousept.mode = 0;


 function  PlotXYDensity(energy,vw)

     if length(vw) < 10
         return;
     end
if length(vw) > 100000
    lprc = 0.01;
    hprc = 99.99;
    sx=3;
    sy=3;
elseif length(vw) > 10000
    lprc = 0.1;
    hprc = 99.9;
    sx=3;
    sy=3;
elseif length(vw) > 1000
    lprc = 1;
    hprc = 99;
    sx=5;
    sy=5;
else
    lprc = 5;
    hprc = 95;
    sx=8;
    sy=8;
end    

erange = [prctile(energy,lprc) prctile(energy,hprc)];
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
mode = 2;
tic;
if mode ==1 % add real gaussian to grid for each
    z = zeros(size(x));
    sx = (diff(erange)/100)^2;
    sy = (diff(vrange)/100)^2;
    for j=1:length(energy)
        z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
    end
elseif mode ==2 %build fine 2-D histogram, then smooth
    [gx,gy] = meshgrid(-10:10,-10:10);

    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(size(x));
    vi = 1+floor(nbins * (vw-vrange(1))/diff(vrange));
    ei = 1+floor(nbins * (energy-erange(1))/diff(erange));
    idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
    for j =idx
        z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
    end
    z = conv2(z,G,'same');
end
toc
pcolor(x,y,z);
shading('interp')


function [x,y,z] = CalcDensity(DATA, expspks, mode)

if ismember(DATA.syncsign,[-1 1])
energy = DATA.Spikes.cx(DATA.sids{1});
vw = DATA.Spikes.cy(DATA.sids{1});
elseif isfield(DATA,'AllClusters')
    energy  = DATA.AllClusters(DATA.probe).cx(expspks);
    vw  = DATA.AllClusters(DATA.probe).cy(expspks);
else
energy = DATA.Spikes.cx(expspks);
vw = DATA.Spikes.cy(expspks);
end
DATA.plot.DensitySigma = [3 3];

if length(vw) > 10000
    lprc = 0.01;
    hprc = 99.99;
elseif length(vw) > 10000
    lprc = 0.1;
    hprc = 99.9;
    
elseif length(vw) > 1000
    lprc = 1;
    hprc = 99;
    DATA.plot.DensitySigma = [5 5];
else
    lprc = 5;
    hprc = 95;
    DATA.plot.DensitySigma = [10 10];
end    

if DATA.xyfig == gcf && mode < 3
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
else
    if prctile(energy,hprc) > prctile(energy,hprc-1) * 3
        erange = [prctile(energy,lprc) prctile(energy,hprc-1)];
    else
        erange = [prctile(energy,lprc) prctile(energy,hprc)];
    end
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
end
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
tic;
if mode ==1 % add real gaussian to grid for each
    z = zeros(size(x));
    sx = (diff(erange)/100)^2;
    sy = (diff(vrange)/100)^2;
    for j=1:length(energy)
        z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
    end
elseif mode ==2 || mode == 3 %build fine 2-D histogram, then smooth
    [gx,gy] = meshgrid(-10:10,-10:10);
    sx= DATA.plot.DensitySigma(1);
    sy= DATA.plot.DensitySigma(2);
    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(size(x));
% ignore spikes where both are set to 0
    idx = find(vw ~=0 | energy ~=0);
    
    vi = 1+floor(nbins * (vw(idx)-vrange(1))/diff(vrange));
    ei = 1+floor(nbins * (energy(idx)-erange(1))/diff(erange));
    idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
    for j =idx
        z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
    end
    z = conv2(z,G,'same');
end
         
function DensityPlot(a,b)
global mousept;

%DATA = combine('getstate');
DATA = GetDataFromFig(a);

if isfield(DATA,'spklist') & length(DATA.spklist) > 10
    expspks = DATA.spklist;
else
    expspks = DATA.spkrange(1):DATA.spkrange(2);
end

    
    figure(DATA.xyfig);
if DATA.densityplot
    DATA.densityplot = 0;
    hold off;
    DATA = DrawXYPlot(DATA,expspks);
%    ClearMouse;
    set(DATA.toplevel,'UserData',DATA);
    SetGui(DATA);
    return;
end

[x,y,z] = CalcDensity(DATA, expspks, 2);

erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
toc
pcolor(x,y,z);
shading('interp')
hold on;
if 0 %code for trackin peak ridge. May help track clusters/auto cut
[a,b] = max(z);
py = smooth(y(b,1),3,'gauss');
plot(x(1,:), py);
base = 0.1;
for j = length(py):-1:1
    dip(j) = a(j)./(base+max(a(j:end)).^2);
end
plot(dip);
end
DATA = DrawClusters(DATA, DATA.cluster, 0);
DATA.densityplot = 1;
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
x = caxis;
set(it,'string',sprintf('%.2f',x(2)));
SetGui(DATA);
set(DATA.toplevel,'UserData',DATA);


function [expname, exptypename, suff] = OldExpt2Name(Expt)
% take an expt and return a name identifying the
% exptype, of the form rds.dxXce  
SpkDefs;
suff = [];
        stimname = stimnames{Expt.Stimvals.st+1};
    if strcmp(Expt.Stimvals.et,'or') && Expt.Stimvals.ei > 100
        Expt.Stimvals.et = 'dir';
    end
    
    if strmatch(Expt.Stimvals.e2, 'e0')
        if strmatch(Expt.Stimvals.e3, 'e0')
        exptypename = Expt.Stimvals.et;
        else
            exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e3];
        end       
%        expname = [Expt.Stimvals.et];
    else
        exptypename = [Expt.Stimvals.et 'X' Expt.Stimvals.e2];
        if ~strcmp(Expt.Stimvals.e3,'e0')
            exptypename = [exptypename 'X' Expt.Stimvals.e3];
        end
    end
    if isfield(Expt.Header,'Options') & strfind(Expt.Header.Options,'+cr')
        exptypename = ['C' exptypename];
    end
    if Expt.Header.psych
        exptypename = [exptypename 'P'];
    end
    if isfield(Expt.Stimvals,'rb') & Expt.Stimvals.rb ~= 0
        exptypename = [exptypename 'RB'];
    end
    if Expt.Header.rc
        suff = 'RC';
    end
    if strmatch(exptypename,{'dxXId' 'orXId' 'dirXId'},'exact')
        jx = GetEval(Expt,'jx');
        if jx > 0
         exptypename = [exptypename 'D'];
        end
    end
    if strncmp(exptypename,'dxXce',5)
        if Expt.Stimvals.n2 > 2
         exptypename = strrep(exptypename,'dxXce','dxXces');
        end
    end
    expname = [stimname '.' exptypename];

function DATA = ListExpts(DATA, Expts);

SpkDefs;
explist = {};
na = 1;
nb = 1;

explist{1} = 'All';
explabel{1} = 'All';

exptypelist = [];
na= 2;
for j = 1:length(Expts)
% expname is set in APlaySpkFile. Expt2Name n use unless it is made a
% separate function to be used by both.
%BUT exptype name IOS calculated fo real here, so they need to match!!
    [expname, exptypename, suff] = Expt2Name(Expts{j});
    expname = [expname suff];
% if this expt is not in the list, add it
    if isempty(strmatch(Expts{j}.Header.expname,explist,'exact'))
        explist{na} = expname;
        if DATA.state.online == 2
            explist{na} = Expts{j}.Header.expname;
%            Expts{j}.Header.expname = expname;
        else
            explist{na} = Expts{j}.Header.expname;
        end
        exptypelist{na} = exptypename;
        DATA.explist = explist;
        DATA.exptypelist = exptypelist;
        cls = '';
        for k = 1:3
            outname = CombinedName(DATA, na,k);
            if exist(outname,'file')
                cls = [cls ' c' num2str(k)];
            end
        end
        if na > length(DATA.allcombineids)
            DATA.allcombineids{na} = [];
        end
        explabel{na} = [explist{na} cls];
        DATA.explist = explist;
        na = na+1;
    else
%        explist{na} = 'unknown';
    end    
end
if isfield(DATA,'clst')
p  = get(DATA.clst,'Listboxtop');
if p > length(explabel)
    set(DATA.clst,'ListboxTop',1);
end
p  = get(DATA.clst,'value');
if p > length(explabel)
    set(DATA.clst,'value',1);
end
set(DATA.clst,'string',explabel);
end
DATA.explist = explist;
DATA.exptypelist = exptypelist;

function [id, ok] = CheckListForName(list,name)

fid = fopen(list,'r');
if fid >0 
    ok = 1;
    names = textscan(fid,'%s');
    files = splitpath(names{1});
    id = strmatch(splitpath(name),files);
    fclose(fid);
else
    ok = 0;
    id = [];
end


function CheckLists(DATA)
       
 pairs = {'ORBW' '/bgc/bgc/anal/orbw/lemORBW.lst' ; ...
     'OTRC' '/bgc/bgc/anal/orbw/otrc.lst'};
        
        [a,dp] = splitpath(DATA.datafilename);
        d = dir(dp);
        nf = 0;
        for j = 1:length(d)
            if strfind(d(j).name,'.mat') 
              for k = 1:size(pairs,1)
                if strfind(d(j).name,pairs{k,1})
                    if isempty(CheckListForName(pairs{k,2},d(j).name))
                    nf = nf+1;
                    dat.files{nf} = [dp '/' d(j).name];
                    dat.lists{nf} = pairs{k,2};
                    else
                        fprintf('%s Already in %s\n',d(j).name,pairs{k,2});
                    end
                    dat.nf = nf;
                end
              end
            end
        end
if nf
    SPACE = 10;
    cw=8;
    ch=10;
    figure('Position',[10 10 cw*60, (ch+SPACE) * (nf+2)]);
    fp = get(gcf,'Position');
    bp = [10 10 fp(4) 20];
    for j = 1:nf
        bp = [10 bp(2)+bp(4)+SPACE fp(3) ch * 2];
        uicontrol(gcf,'Style', 'checkbox',...
            'String', [dat.files{j} dat.lists{j}], 'Tag', 'SuffList', 'Position', bp,'UserData',j);
    end
    bp = [10 10 cw * 5 20];
        uicontrol(gcf,'Style', 'pushbutton',...
            'String','go', 'callback',@cplists);
     set(gcf,'UserData',dat);
    
end


function cplists(caller,b)

    lfig = get(caller,'parent');
    dat = get(lfig,'UserData');
    it = findobj(lfig,'Tag','SuffList');
    for j = 1:length(it)
        go = get(it(j),'value');
        if go
            fprintf('Adding %s to %s\n',dat.files{j},dat.lists{j});
            AddToList(dat.lists{j},dat.files{j});
        end
    end
    close(lfig);
 
function DATA = ListSubExpts(DATA, id, varargin)

setv = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'relist',4)
        setv = get(DATA.elst,'value');
    end
    j = j+1;
end
suffs = DATA.suffs;
na = 1;
subexplist = {};
expid = [];
nrp = zeros(size(DATA.explist));
for j = 1:length(DATA.Expts);
    eid = strmatch(DATA.Expts{j}.Header.expname,{DATA.explist{id}},'exact');
    if( id == 1 | ~isempty(eid)) & isfield(DATA.Expts{j},'Trials')
        tid = regexp(DATA.Expts{j}.Header.Name,'Expt[0-9]*.mat');
        if ~isempty(tid) %online file
            label = sprintf(' (%s:%d-%d)',DATA.Expts{j}.Header.Name(tid:end-4),...
                DATA.Expts{j}.Trials(1).Trial,...
                DATA.Expts{j}.Trials(end).Trial');
        elseif DATA.show.times
            label = sprintf(' %.1f-%.1f (id %d-%d)',DATA.Expts{j}.Trials(1).Start(1)./10000,...
                DATA.Expts{j}.Trials(end).Start(1)./10000,DATA.Expts{j}.Trials(1).id,DATA.Expts{j}.Trials(end).id);
        else
            label = sprintf(' %d-%d (id %d-%d)',DATA.Expts{j}.Trials(1).Trial,...
                DATA.Expts{j}.Trials(end).Trial,DATA.Expts{j}.Trials(1).id,DATA.Expts{j}.Trials(end).id);
        end
        if DATA.Expts{j}.Header.psych
            label = [label ' P'];
        end
        if DATA.Expts{j}.Header.rc
            label = [label ' RC'];
        end
        label = [label ShowString(DATA,DATA.Expts{j})];
        if id == 1
            expi = strmatch(DATA.Expts{j}.Header.expname, DATA.explist,'exact');
            nrp(expi) = nrp(expi)+1;
            subexplist{na} = [DATA.Expts{j}.Header.expname num2str(na) label];
        else
            subexplist{na} = [DATA.explist{id(eid)} suffs(mod(na-1,length(suffs))+1) label];
            subexplist{na} = [DATA.explist{id(eid)} num2str(j) label];
        end
        DATA.explabels{j} = subexplist{na};
        if DATA.Expts{j}.gui.clustertype == 0
            subexplist{na} = [subexplist{na} '*'];
        elseif DATA.Expts{j}.gui.clustertype == 2
            subexplist{na} = [subexplist{na} '(O)'];
        elseif DATA.Expts{j}.gui.ncluster == 0
            subexplist{na} = [subexplist{na} '(Z)'];            
        end
        expid(na) = j;
        na = na+1;
    elseif strmatch('unknown',{DATA.explist{id}})
        subexplist{na} = [DATA.explist{id} suffs(na)];
        expid(na) = j;
        na = na+1;        
    end
end
p  = get(DATA.elst,'Listboxtop');
if p > length(subexplist)
    set(DATA.elst,'ListboxTop',1);
end
    set(DATA.elst,'string',subexplist,'value',setv);
    DATA.expid = expid;
    DATA.subexplist = subexplist;
    if ~isfield(DATA,'currentexpt')
        DATA.currentexpt = 1;
    end

    
function timerfna(tim, varargin)
    DATA = get(findobj('Tag',get(tim,'Tag')),'UserData');
%DATA = get(findobj('Tag','Combiner'),'UserData');
if DATA.state.autolist
tic;
    combine('relist');
%toc
end

function timerfn(tim, varargin)
    DATA = get(findobj('Tag',get(tim,'Tag')),'UserData');
    d = dir(DATA.datafilename);
    if d.datenum > DATA.lastread || d.bytes > DATA.lastsize;
 %       fprintf('Reading from %d\n',DATA.linesread);
        [Expts, DATA.linesread] = ReadOnlineTxt(DATA.datafilename, DATA);
%        fprintf('Read %d expts\n',length(Expts));
         Expts = CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
%         fprintf('List\n');
         if length(Expts) > length(DATA.Expts)
             DATA.Expts = Expts;
             DATA = ListExpts(DATA,Expts);
 %            fprintf('listexp\n');
             DATA = combine('listexps',DATA,'Tag',DATA.tag.top);
%             fprintf('Set\n');
             set(DATA.elst,'value',length(get(DATA.elst,'string')))
             fprintf('Replot (%d,%d) at %s\n',length(Expts),DATA.linesread,datestr(now));
         else
             DATA.Expts = Expts;
             fprintf('Replot (%d,%d,%d) at %s\n',length(Expts),length(Expts{end}.Trials),DATA.linesread,datestr(now));
         end
        DATA.lastread = d.datenum;
        DATA.lastsize = d.bytes;
        combine('setexp',DATA,'Tag',DATA.tag.top);
%        fprintf('Expt Set\n');
    end

    
function LoadAllProbes(a,b, varargin)
    DATA = GetDataFromFig(a);
    expid = DATA.currentexpt;
    j = 1; 
    while j <= length(varargin)
        j = j+1;
    end
    times(1) = DATA.Expts{expid}.Trials(1).Start(1)-10000;
    times(2) = DATA.Expts{expid}.Trials(end).End(end) + 10000;
    DATA.allexp = expid;
    Spikes = DATA.AllData.Spikes;
    tic;
    DATA = GetAllProbeFig(DATA);
    if ~isfield(DATA,'ptsize')
        DATA.ptsize = 1;
    end
    for j = 1:length(DATA.probelist)
        if DATA.state.online == 0
            DATA.AllData.Spikes = GetProbeFiles(DATA, DATA.probelist(j),DATA.subprobe,'trange',times/10000);
        else
            filename = ['C:' DATA.Expts{DATA.currentexpt}.Header.Name];
            if DATA.probelist(j) > 16 || DATA.probesource(j) == 1
                filename = strrep(filename,'/Expt','A/Expt');
            end
            [DATA.AllData.Spikes]= GetProbeSpikes(DATA.AllData, filename , DATA.probevars{j},DATA.subprobe);
        end
        if isempty(DATA.AllData.Spikes) | isempty(DATA.Expts{expid}.gui.spks)
        DATA.clustervals(j).x = [];
        DATA.clustervals(j).y = [];
        DATA.clustervals(j).spkrange = [0 0 ];
        else
        DATA = SetExptSpikes(DATA, DATA.currentexpt, 0);
        ispk = DATA.Expts{expid}.gui.spks;
        DATA.clustervals(j).x = DATA.Spikes.cx(ispk);
        DATA.clustervals(j).y = DATA.Spikes.cy(ispk);
        DATA.clustervals(j).spkrange = [ispk(1) ispk(end)];
        subplot(6,4,j);
        if DATA.densityplot
        PlotXYDensity(DATA.Spikes.cx(ispk),DATA.Spikes.cy(ispk));slave
      
        else
        plot(DATA.Spikes.cx(ispk),DATA.Spikes.cy(ispk),'.','markersize',DATA.ptsize);
        end
        title(sprintf('%d: %d',j,length(ispk)));
        end
    end
    subplot(6,4,2);
    text(1,1.5,sprintf('Expt %d',DATA.currentexpt),'units','norm');
    toc
    DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
    set(DATA.toplevel,'UserData',DATA);

function DATA = GetAllProbeFig(DATA)
    if DATA.plot.prettyfigs
        GetFigure('AllProbeXYFig');
    elseif isfigure(DATA.figs.allprobes)
        GetFigure(DATA.figs.allprobes);
    else
        [DATA.figs.allprobes, isnew] = GetFigure(DATA.tag.allprobes);
        dat.parentfigtag = DATA.tag.top;
        set(DATA.figs.allprobes,'UserData',dat);
        if isnew
            bp = [10 10 DATA.plot.cw*7 20];
                uicontrol(DATA.figs.allprobes,'Style', 'checkbox',...
        'String', 'Density', 'Tag', 'AllDensity', 'Position', bp,'value',0,...
        'Callback',@Update);

         abp = bp;
         bp(2) = bp(2)+DATA.plot.ch * 2.2;
         uicontrol(DATA.figs.allprobes,'Style', 'pushbutton','String','Set','Position',bp,'Callback', @SetExptClusters);
         bp(1) = bp(1)+bp(3);
         bp(3) = DATA.plot.cw * 6;
         uicontrol(DATA.figs.allprobes,'Style', 'pushbutton', 'Callback', @NextList, ...
             'String', 'Next', 'Position',bp);
         bp = abp;
       for j = 1:length(DATA.probelist)
           bp(1) = bp(1)+bp(3);
           if j < 10
               bp(3) = DATA.plot.cw * 3;
           else
               bp(3) = DATA.plot.cw * 4;
           end
         uicontrol(DATA.figs.allprobes,'Style', 'checkbox',...
        'String', sprintf('%d',j), 'Tag', sprintf('MarkProbe%d',j), 'Position', bp,'value',0,...
        'Callback',@Update);
       end
       bp(1) = bp(1)+bp(3);
       bp(3) = DATA.plot.cw * 5;
       
       uicontrol(DATA.figs.allprobes,'Style', 'pushbutton',...
        'String', 'Mark', 'Position', bp,'value',0,...
        'Callback',@MarkProbes,'UserData',dat);

    bp(1) = bp(1)+bp(3);
       bp(3) = DATA.plot.cw * 5;
       uicontrol(DATA.figs.allprobes,'Style', 'pushbutton',...
        'String', 'Hide SpikeV', 'Position', bp,'value',0,...
        'Callback',@HideSpikes,'UserData',dat);
        end
    end
    
    
function PlotAllExptsAndProbes(a,b, varargin)
    if isfield(a, 'state') && isfield(a,'AllData')
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end
    tstart = now;
    
    DATA.state.autoplotnewprobe = 0;

    for j = 1:length(DATA.probelist)
        DATA = SetProbe(DATA, DATA.probelist(j));
        PlotAllExpts(DATA, 0, 'figlabel', ['AllExpts Probe ' num2str(j)]);
    end
    mytoc(tstart);
    DATA.state.autoplotnewprobe = 1;
    
function PlotSubspace(a,b, varargin)
    if isfield(a, 'state') && isfield(a,'AllData')
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end
    oldrc = DATA.plot.condenseRC;
    DATA.plot.condenseRC = 0;
    PlotCombined(DATA, DATA.Expt);
    DATA.plot.condenseRC = oldrc;
    
function SpoolAllExpts(a,b, toplevel, type)
DATA = get(toplevel,'UserData');
if type == 1
    skip = 10;
elseif type == 2
    skip = 20;
elseif type == 3
    skip = 50;
end    
h = 0;
stop = findobj(DATA.svfig,'Tag','StopSpool');

for j = 1:length(DATA.Expts)
    DATA.currentexpt = j;
   trials = [1:skip:length(DATA.Expts{j}.Trials) length(DATA.Expts{j}.Trials)];
    for k = trials
        PlayOneTrial(DATA,k,0,'nosetgui');
        if h == 0
            h = text(10,0,'trial','color','r');
        else
            set(h,'string',sprintf('%d',DATA.Expts{j}.Trials(k).Trial));
        end
    end
    if get(stop,'value') > 0 
        break;
    end
end
delete(h);

function PlotAllExpts(a,b, varargin)
    if isfield(a, 'state') && isfield(a,'AllData')
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end
    figlabela = 'AllExpts';
    DATA.firsttrial = 0;
    DATA.lasttrial = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'figlabel',4)
            j = j+1;
            figlabela = varargin{j};
        end
        j = j+1;
    end
    nex = length(DATA.Expts);
    [nr, nc] = Nsubplots(nex);
    pw = 0.9./nc;
    ph = (1-0.02.*nr)./nr;
    hstep = 0.02;
    vstep = 0.02;
    set(DATA.toplevel,'Name','Busy......');
    drawnow;

    DATA.ptsize = 1;
    [f, isnew]  = GetFigure(figlabela);
    if isnew  
        hm = uimenu(f,'Label','Spool','Tag','XClusterMenu');
        uimenu(hm,'Label','1/10','Callback',{@SpoolAllExpts, DATA.toplevel, 1});
        uimenu(hm,'Label','1/20','Callback',{@SpoolAllExpts, DATA.toplevel, 2});
        uimenu(hm,'Label','1/50','Callback',{@SpoolAllExpts, DATA.toplevel, 3});
    end

    tstart = now;
    for xi = 1:nex
        y = vstep + floor((xi-1)/nc) .* (ph+vstep);
        x = hstep + (mod(xi-1, nc))./(nc * 1.02);
    if isfield(DATA.Expts{xi}.gui,'spks')  %force, in case reloaded clusters
        [DATA, ispk] = SetExptSpikes(DATA,xi, 0,'useexpall');
    elseif 1
        [DATA, ispk] = SetExptSpikes(DATA,xi, 0,'useexp');
    end
    DATA.nclusters = CountClusters(DATA.cluster);
    set(0,'CurrentFigure',f);
    subplot('Position' ,[x y pw ph]);
    hold off;
    DrawXYPlot(DATA, DATA.Expts{xi}.gui.spks);
    xl = get(gca,'Xlim');
    yl = get(gca,'ylim');
    if iscluster(DATA.cluster,1,DATA.probe) & isfield(DATA.cluster{1,DATA.probe},'dprime')
        dp = DATA.cluster{1,DATA.probe}.dprime;
        text(xl(1)+diff(xl)/5,yl(2)-diff(yl)/10,sprintf('%d:%.1f',xi,dp),'color','r');
    else
        text(xl(1)+diff(xl)/5,yl(2)-diff(yl)/10,sprintf('%d ed=%.1f',xi,GetEval(DATA.Expts{xi},'ed')),'color','r');
    end
 %   drawnow;
    xlabel('');
    ylabel('');
    title('');
    end
    mytoc(tstart);
NotBusy(DATA);
set(DATA.toplevel,'UserData',DATA);

function MarkProbes(a,b)
    DATA = GetDataFromFig(a);
    for j = 1:length(DATA.probelist)
        it = findobj(get(a,'parent'),'Tag',sprintf('MarkProbe%d',j));
        if length(it) == 1 && get(it,'value') > 0
            DATA.markexp(DATA.currentexpt,j) = 1;
            DATA.plot.useprobe(j) = 1;
        else
            DATA.plot.useprobe(j) = 0;
        end
    end
    imagesc(DATA.markexp);
    set(DATA.toplevel,'UserData',DATA);
    
function HideSpikes(a,b);
    DATA = GetDataFromFig(a);
    x = GetCheck('NoSpikes')
    DATA.state.nospikes = ~x;
    SetCheck('NoSpikes',DATA.state.nospikes);
    if x
        set(a,'String','Show Waveforms');
    else
        set(a,'String','Hide Waveforms');
    end
    set(DATA.toplevel,'UserData',DATA);
      
        
function LoadAllProbeParams(a,b, varargin)
    DATA = GetDataFromFig(a);
    expid = DATA.currentexpt;
    j = 1; 
    while j <= length(varargin)
        j = j+1;
    end
    DATA.allexp = expid;
    Spikes = DATA.AllData.Spikes;
    tic;
    DATA = GetAllProbeFig(DATA);
    if ~isfield(DATA,'ptsize')
        DATA.ptsize = 1;
    end
    for j = 1:length(DATA.probelist)
        set(DATA.toplevel,'name',sprintf('Loading %d',DATA.probelist(j)));
        drawnow;
        if DATA.state.online == 0
            DATA.AllData.Spikes = GetProbeFiles(DATA, DATA.probelist(j),DATA.subprobe);
        else
            filename = ['C:' DATA.Expts{DATA.currentexpt}.Header.Name];
            if DATA.probelist(j) > 16
                filename = strrep(filename,'/Expt','A/Expt');
            end
            [DATA.AllData.Spikes]= GetProbeSpikes(DATA.AllData, filename , DATA.probevars{j},DATA.subprobe);
        end
        nspk = length(DATA.AllData.Spikes.times);
        DATA = CalcClusterVars(DATA,  1:nspk,'force');
        DATA.AllClusters(j).cx = DATA.Spikes.cx(1:nspk);
        DATA.AllClusters(j).cy = DATA.Spikes.cy(1:nspk);
        DATA.AllClusters(j).times = DATA.AllData.Spikes.times;
        DATA.AllClusters(j).codes = DATA.AllData.Spikes.codes(:,2);
    end
    DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
    set(DATA.toplevel,'UserData',DATA);
    NotBusy(DATA);
    
    
function DATA = LoadAllProbeSpikes(a,b, varargin)
    if isstruct(a)
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end
    expid = DATA.currentexpt;
    times(1) = DATA.Expts{expid}.Trials(1).Start(1)-10000;
    times(2) = DATA.Expts{expid}.Trials(end).End(end) + 10000;
    DATA.allexp = expid;
    DATA.AllData.Spikes = [];
    probelist = DATA.probelist;
    j = 1;
    while j <= length(varargin)
%        'select' load just selected probes for selected expt range
        if strncmpi(varargin{j},'select',5)
            AddProbeList(DATA);
            if ~isfield(DATA.plot,'useprobe')
                return;
            end
            probelist = find(DATA.plot.useprobe);
            if isempty(probelist)
                return;
            end
            if sum(probelist == DATA.probe) == 0
                DATA.probe = probelist(1);
            end
%            DATA.probe = max(probelist);
            eid = DATA.exabsid;
            times(1) = DATA.Expts{eid(1)}.Trials(1).Start(1) - 10000;
            times(2) = DATA.Expts{eid(end)}.Trials(end).End(end) - 10000;
        elseif  strncmpi(varargin{j},'allselect',6) %setlect all
            DATA.plot.useprobe = ones(size(DATA.probelist));
            AddProbeList(DATA);
            probelist = find(DATA.plot.useprobe);
            eid = DATA.exabsid;
            times(1) = DATA.Expts{eid(1)}.Trials(1).Start(1) - 10000;
            times(2) = DATA.Expts{eid(end)}.Trials(end).End(end) - 10000;
        end
        j = j+1;
    end
    tic;
    DATA = GetAllProbeFig(DATA);
    if ~isfield(DATA,'ptsize')
        DATA.ptsize = 1;
    end
    if isfield(DATA,'AllClusters') %% Can't have both
        DATA = rmfield(DATA,'AllClusters');
    end
    for j = 1:length(probelist)
        set(DATA.toplevel,'Name',sprintf('Loading probe %d',probelist(j)));
        drawnow;
        if DATA.state.online == 0
            DATA.AllSpikes{probelist(j)} = GetProbeFiles(DATA, probelist(j),DATA.subprobe,'trange',times/10000,'nodv');
%            DATA.AllSpikes{probelist(j)} = GetProbeFiles(DATA, probelist(j),'trange',times/10000);
            DATA = LoadClusters(DATA, ClusterFile(DATA,'probe', probelist(j)),DATA.subprobe,'probe', probelist(j));
        else
            filename = ['C:' DATA.Expts{DATA.currentexpt}.Header.Name];
            if DATA.probelist(j) > 16
                filename = strrep(filename,'/Expt','A/Expt');
            end
            DATA.AllSpikes{probelist(j)} = GetProbeSpikes(DATA.AllData, filename , DATA.probevars{probelist(j)},DATA.subprobe);
            DATA.AllSpikes{probelist(j)}.firstspki = 1;
        end
    end
    toc
    if DATA.state.recut
        tic
        fprintf('Classifying Spikes');
        set(DATA.toplevel,'Name',sprintf('Classifying Spikes'));
        drawnow;
        DATA = ReClassifyAll(DATA);
        toc;
    end
    
%   DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
    set(DATA.toplevel,'UserData',DATA);
    AddProbeList(DATA);
    if isfield(DATA,'currentexpt') & ismember(DATA.probe, probelist)
        id = find(probelist == DATA.probe);
        if id < length(probelist)
            id = probelist(id+1);
        elseif id > 1
            id = probelist(id-1);
        else
            id = DATA.probe;
        end
        xc = CalcXcorr(DATA,DATA.currentexpt,DATA.probe,id);
    end

    function LoadSpikeTimes(a,b, varargin)

        DATA = GetDataFromFig(a);
        eid = DATA.currentexpt;
        ns = ones(size(DATA.probelist));
        for j = 1:length(DATA.probelist)
            DATA.AllSpikes{j}.times = [];
        end
        for e = 1:length(DATA.Expts)
            DATA.currentexpt = e;
            sfile = ClusterFile(DATA, e, 'codes');
            if exist(sfile,'file')
                load(sfile);
                for j = 1:length(DATA.probelist)
                    n = length(codes{j});
                    if n
                    DATA.AllSpikes{j}.codes(ns(j):ns(j)+n-1,2) = codes{j};
                    DATA.AllSpikes{j}.times(ns(j):ns(j)+n-1,1) = times{j};
                    ns(j) = ns(j)+n;
                    end
                end
            end
            DATA.Expts{e}.gui.classified = 1;
        end
        for j = 1:length(DATA.probelist)
            if isfield(DATA.AllSpikes{j},'codes')
            ns(j) = sum(DATA.AllSpikes{j}.codes(:,2) > 0);
            else
                ns(j) = NaN;
            end
        end
       DATA.currentexpt = eid;
       DATA.state.nospikes = 1;
       DATA.state.recount = 1;
       DATA.state.showspikes = 0;
       DATA.state.showspkxy = 0;
       SetGui(DATA);
        set(DATA.toplevel,'UserData',DATA);
        
    function xc = CalcXcorrDC(DATA, eids, sa, sb)
    dt = 2;
    colors = 'bmr';
    fprintf('Calculating DC Cross Correlation %d,%d\n',sa,sb);
    tic;
    for j = 1:length(eids)
        Expt =  DATA.Expts{eids(j)};
        trange = [Expt.Trials(1).Start(1) Expt.Trials(end).End(end)];
        if DATA.syncsign > 0
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
            DATA.AllSpikes{sa}.values(:,9) > 0);
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
            DATA.AllSpikes{sb}.values(:,9) > 0);
        elseif DATA.syncsign < 0
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
            DATA.AllSpikes{sa}.values(:,9) < 0);
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
            DATA.AllSpikes{sb}.values(:,9) < 0);
        else 
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2));
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2));
        end

        [ai,bi] = FindSync(DATA.AllSpikes{sa}.times(aid),...
                DATA.AllSpikes{sb}.times(bid),dt);
            dcs(1,:) = mean(DATA.AllSpikes{sa}.values(aid(ai),24:end),2);
            dcs(2,:) = mean(DATA.AllSpikes{sb}.values(bid(bi),24:end),2);
    end
    scatter(dcs(1,:),dcs(2,:),'.');
    

function DATA = ReClassifyAll(DATA, varargin)
% apply defind clusters to all expts/probes
% useful if saved list of classification codes is missing/ wrong
% need to run this when loading all spikes (because, err....)
mkmean = 0;
    j = 1;
     probes = DATA.probelist;
    while j <= length(varargin)
        if strncmpi(varargin{j},'mkmean',5)
            mkmean = 1;
        elseif strncmpi(varargin{j},'probes',5)
            j = j+1;
            probes = varargin{j};
        end
        j = j+1;
    end
if 0
    [DATA, ispk] = SetExptSpikes(DATA, DATA.currentexpt,'setrange');
    DATA = SetSpkCodes(DATA, ispk, DATA.probe, 0);
end
    oldprobe = DATA.probe;
    ts = now;
    if length(probes) == 1 && DATA.probe == probes(1)
        DATA = CalcClusterVars(DATA, 1:length(DATA.AllData.Spikes.times));
    end
    for k = 1:length(probes)
        if ~isfield(DATA,'AllSpikes')
            if DATA.probe ~= probes(k) || length(probes) > 1
            DATA = SetProbe(DATA, probes(k));
            DATA = CalcClusterVars(DATA, 1:length(DATA.AllData.Spikes.times));
            end
        elseif k <= length(DATA.AllSpikes) && ~isempty(DATA.AllSpikes{k})
            DATA.probe = k;
            DATA = CalcClusterVars(DATA, 1:length(DATA.AllSpikes{k}.times));
        end
        if DATA.probe == probes(k) % have spikes for this
            ts = [ts now];
            for j = 1:length(DATA.Expts)
                if isfield(DATA.Expts{j},'Cluster')
                    DATA.cluster = DATA.Expts{j}.Cluster;
                    if iscluster(DATA.cluster, 1, k) == 1
                        [DATA, ispk] = SetExptSpikes(DATA, j,'setrange');
                        tic;
                        DATA = SetSpkCodes(DATA, ispk, DATA.probe, 0);
                    end
                end
                ts = [ts now];
                ctimes(k,j) = ts(end)-ts(end-1);
            end
            if mkmean
                DATA = CalcMeanSpike(DATA,1:length(DATA.Expts));
            end
        end
    end
    DATA.probe = oldprobe;
        
    
function xc = CalcXcorrV(DATA, eids, sa, sb)
    dt = 2;
    pts = 24:32;
    colors = 'bmr';
    fprintf('Calculating V Cross Correlation %d,%d\n',sa,sb);
    tic;
    for j = 1:length(eids)
        Expt =  DATA.Expts{eids(j)};
        trange = [Expt.Trials(1).Start(1) Expt.Trials(end).End(end)];
        if DATA.syncsign > 0
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
            DATA.AllSpikes{sa}.values(:,9) > 0);
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
            DATA.AllSpikes{sb}.values(:,9) > 0);
        elseif DATA.syncsign < 0
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
            DATA.AllSpikes{sa}.values(:,9) < 0);
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
            DATA.AllSpikes{sb}.values(:,9) < 0);
        else 
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2));
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2));
        end
        [ai,bi] = FindSync(DATA.AllSpikes{sa}.times(aid),...
                DATA.AllSpikes{sb}.times(bid),dt);
            dcs(1,:,:) = DATA.AllSpikes{sa}.values(aid(ai),:);
            dcs(2,:,:) = DATA.AllSpikes{sb}.values(bid(bi),:);
    end
    if ~isempty(pts)
        a = mean(dcs(1,:,pts),3);
        b = mean(dcs(2,:,pts),3);
        scatter(a(:),b(:),'.');
    else
        scatter(dcs(1,:),dcs(2,:),'.');
    end
    

    
    function xc = CalcXcorr(DATA, eids, sa, sb)
    dt = [-1000:5:1000];
    sign = DATA.syncsign;
    colors = 'bmr';
    if DATA.spikelist(1) >= 0
        spikelist = DATA.spikelist;
    else
       spikelist = [0:4];
    end
    xc = [];
    fprintf('Calculating Cross Correlation %d,%d\n',sa,sb);
    tic;
    for j = 1:length(eids)
        Expt =  DATA.Expts{eids(j)};
        trange = [Expt.Trials(1).Start(1) Expt.Trials(end).End(end)];
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2));
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2));

        isspk = [0 0];
        if sum(DATA.AllSpikes{sa}.codes(aid,2)) == 0 %% no clusters cut
            aspikelist = 0;
        else
            aspikelist = spikelist;
            isspk(1) = 1;
        end
        if sum(DATA.AllSpikes{sb}.codes(bid,2)) == 0 %% no clusters cut
            bspikelist = 0;
        else
            bspikelist = spikelist;
            isspk(2) = 1;
        end
        
        for t = 1:length(Expt.Trials)
            durs(t) = Expt.Trials(t).End(end) - Expt.Trials(t).Start(1);
        end
        ptimes = 0:10:mean(durs);
        for t = 1:length(Expt.Trials)
        trange = [Expt.Trials(t).Start(1)-0.1 Expt.Trials(t).End(end)+0.1];
        if sign == 1 %use only positive triggers
        aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
            DATA.AllSpikes{sa}.values(:,9) > 0 & ismember(DATA.AllSpikes{sa}.codes(:,2),aspikelist));
        at = DATA.AllSpikes{sa}.times(aid);
        bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
            DATA.AllSpikes{sb}.values(:,9) > 0 & ismember(DATA.AllSpikes{sb}.codes(:,2),bspikelist));
        bt = DATA.AllSpikes{sb}.times(bid);
        elseif sign < 0 %use only negative triggers
            aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ...
                DATA.AllSpikes{sa}.values(:,9) < 0 & ismember(DATA.AllSpikes{sa}.codes(:,2),aspikelist));
            at = DATA.AllSpikes{sa}.times(aid);
            bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ...
                DATA.AllSpikes{sb}.values(:,9) < 0 & ismember(DATA.AllSpikes{sb}.codes(:,2),aspikelist));
            bt = DATA.AllSpikes{sb}.times(bid);
        else
            aid = find(DATA.AllSpikes{sa}.times > trange(1) & DATA.AllSpikes{sa}.times < trange(2) & ismember(DATA.AllSpikes{sa}.codes(:,2),aspikelist));
            at = DATA.AllSpikes{sa}.times(aid);
            bid = find(DATA.AllSpikes{sb}.times > trange(1) & DATA.AllSpikes{sb}.times < trange(2) & ismember(DATA.AllSpikes{sb}.codes(:,2),aspikelist));
            bt = DATA.AllSpikes{sb}.times(bid);
        end
        dts = [];
        if length(at) & length(bt)
          for k = 1:length(at)
              dts(k,:) = hist(at(k) - bt,dt);
          end
          xc(t,:) = sum(dts,1);
          psth{1}(t,:) = hist(at-Expt.Trials(t).Start(1),ptimes);
          psth{2}(t,:) = hist(bt-Expt.Trials(t).Start(1),ptimes);
        end
        counts(t,:) = [length(at) length(bt)];
        end
    end
    if length(xc)
    xc = sum(xc);
    psth{1} = sum(psth{1},1);
    psth{2} = sum(psth{2},1);
    fprintf('%d,%d spikes, took %.2f (%s)\n',sum(psth{1}),sum(psth{2}),toc,datestr(now));
    if sa == sb
        xc(find(dt == 0)) = 0;
    end
    bar(dt(2:end-1)./10,xc(2:end-1),1,colors(sum(isspk)+1));
fgf
    end

function  TrackTemplates(DATA)

    DATA.TemplateScores = [];
probelist = DATA.probelist;
for p = 1:length(probelist)
    fprintf('Probe %d..',p);
    Spks = GetProbeFiles(DATA,p,DATA.subprobe);
    scores = [];
    for j = 1:size(DATA.Templates,1)
        scores(j,:) = Spks.values * DATA.Templates(j,:)';
    end
    for eid = DATA.exabsid;
        trials = [DATA.Expts{eid}.Trials.Trial];
        espk = find(Spks.times> DATA.Expts{eid}.Trials(1).Start(1) & ...
                Spks.times < DATA.Expts{eid}.Trials(end).End(end));
        for j = 1:length(DATA.Expts{eid}.Trials)
            times = [DATA.Expts{eid}.Trials(j).Start(1) DATA.Expts{eid}.Trials(j).End(end)];
            tspk = find(Spks.times(espk) > times(1) & ...
                Spks.times(espk) < times(2));
            if length(tspk) && ~isempty(scores)
                DATA.TemplateScores(p,:,trials(j)) = max(scores(:,espk(tspk)),[],2);
            end
        end
    end
end
GetFigure('TemplateScores');
set(DATA.toplevel,'UserData',DATA);

    
function lProbeClusters(a,b)
    if isstruct(a)
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end
    for j = 1:length(DATA.Expts)
        DATA.allexp = j;
        PlotAllProbeXY(DATA);
        drawnow;
    end
    
function DATA = PlotAllProbeXY(a,b)

    if isstruct(a)
        DATA = a;
    else
    DATA = GetDataFromFig(a);
    end
    cx = [];
    cy = [];
    DATA = GetAllProbeFig(DATA);
    set(gcf, 'KeyPressFcn',@KeyPressed);
set(gcf, 'WindowButtonDownFcn',@ButtonPressed);
set(gcf, 'WindowButtonMotionFcn',@ButtonDragged);
set(gcf, 'WindowButtonUpFcn',@ButtonReleased);
    densityplot = GetCheck('AllDensity',DATA.figs.allprobes);
    DATA.alldensityplot = densityplot;
    if isfield(DATA.plot,'useprobe') & sum(DATA.plot.useprobe) > 1
        probelist = find(DATA.plot.useprobe);
    else
        probelist = DATA.probe;
    end
    [nr,nc] = Nsubplots(length(probelist));
    if length(probelist) == 24
        nr = 6;
        nc = 4;
    end
    times = [DATA.Expts{DATA.allexp}.Header.trange];
    DATA.ptsize = 1;
    oldprobe = DATA.probe;
    probes = find(DATA.plot.useprobe);
    for j = 1:length(probes)
        DATA.probe = probes(j);
        p = DATA.probe;
        if isfield(DATA,'AllSpikes') 
            ispk = FindSpikes(DATA, times, p,[]);
            if length(ispk)
            DATA = CalcClusterVars(DATA,  ispk,'probe',p);
%            DATA.AllSpikes{p}.cx = DATA.Spikes.cx;
%            DATA.AllSpikes{p}.cy = DATA.Spikes.cy;
            cx = DATA.Spikes.cx(ispk);
            cy = DATA.Spikes.cy(ispk);
            DATA.clustervals(p).spkrange = [ispk(1) ispk(end)];
            DATA.AllSpikes{p}.spklist = ispk;
            end
        elseif isfield(DATA,'AllClusters')
            ispk = find(DATA.AllClusters(j).times > times(1) &...
                DATA.AllClusters(j).times < times(2));
            cx = DATA.AllClusters(j).cx(ispk);
            cy = DATA.AllClusters(j).cy(ispk);
            DATA.AllClusters(j).spklist = ispk;
        else
            ispk = DATA.clustervals(j).spkrange(1):DATA.clustervals(j).spkrange(2);
            cx = DATA.clustervals(j).x;
            cy = DATA.clustervals(j).y;
        end
        subplot(nr,nc,j);
        hold off;
        dprime = 0;
        if densityplot
            PlotXYDensity(cx,cy);
        else
            [a, dprime] = SetSpkCodes(DATA,ispk,DATA.probe,2);
%            plot(cx,cy,'.','markersize',1);
        end
        set(gca,'UserData',DATA.probelist(j));
        hold on;
        if iscluster(DATA.cluster,1,p)
            DATA = DrawClusters(DATA, DATA.cluster,0);
            [xr, yr] = ClusterRange(DATA.cluster,p);
        else
            xr = [0 0];
            yr = [0 0];
        end
        if DATA.plot.autoscale && length(cx)
            %or max cluster
            DATA = SetXYRanges(DATA,cx, cy);
        end
            set(gca,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
       if DATA.plot.prettyfigs
        title(sprintf('Probe %d',p));
       else
        title(sprintf('%d: %d,%.2f',p,length(ispk),dprime));
       end
    end
    DATA.probe = oldprobe;
    subplot(nr,nc,1);
    text(2,1.4,sprintf('Expt %d ids %d-%d, Trials %d-%d',DATA.allexp,DATA.Expts{DATA.allexp}.Trials(1).id,...
        DATA.Expts{DATA.allexp}.Trials(end).id,DATA.Expts{DATA.allexp}.Trials(1).Trial,...
        DATA.Expts{DATA.allexp}.Trials(end).Trial),'units','norm');

function DATA = SetXYRanges(DATA, cx, cy)

    
if nargin == 1
    cx = DATA.Spikes.cx(DATA.spklist);
    cy = DATA.Spikes.cy(DATA.spklist);
end

if DATA.plot.autoscale == 0
    return;
elseif DATA.plot.autoscale == 1 || length(cx) == 0
    DATA.plot.clusterXrange  = get(gca,'Xlim');
    DATA.plot.clusterYrange  = get(gca,'Ylim');
elseif ismember(DATA.plot.autoscale,[2 3 4])
    if min(cx) > 0
        minx = 0;
    else
        minx = prctile(cx,1 * 1.1);
    end
    if min(cy) > 0
        miny = 0;
    else
        miny = prctile(cy,1 * 1.1);
    end
    if DATA.plot.autoscale == 2
    DATA.plot.clusterYrange = [miny prctile(cy,99.9).*1.2];
    DATA.plot.clusterXrange = [minx prctile(cx,99.9).*1.2];
    elseif DATA.plot.autoscale == 3
    DATA.plot.clusterYrange = [miny prctile(cy,99.5).*1.2];
    DATA.plot.clusterXrange = [minx prctile(cx,99.5).*1.2];
    else
    DATA.plot.clusterYrange = [miny prctile(cy,95).*2];
    DATA.plot.clusterXrange = [minx prctile(cx,95).*2];
    end
    if min(cx) < 0
        DATA.plot.clusterXrange(1) = prctile(cx,1) * 1.1;
    end
    if min(cy) < 0
        DATA.plot.clusterYrange(1) = prctile(cy,1) * 1.1;
    end
elseif DATA.plot.autoscale == 4
    DATA.plot.clusterYrange = [min(cy) max(cy)];
    DATA.plot.clusterXrange = [min(cx) max(cx)];
end


function [xr, yr] = ClusterRange(C, p)
    
    for j = 1:size(C,1)
        if isfield(C{j,p},'x') & C{j,p}.x(2) > 0
        xh(j) = C{j,p}.x(1)+C{j,p}.x(2);
        xl(j) = C{j,p}.x(1)-C{j,p}.x(2);
        yh(j) = C{j,p}.y(1)+C{j,p}.y(2);
        yl(j) = C{j,p}.y(1)-C{j,p}.y(2);
        else
            xh(j) = 0;
            xl(j) = 0;
            yl(j) = 0;
            yh(j) = 0;
        end
    end
    xr = [min(xl) max(xh)];
    yr = [min(yl) max(yh)];
        
function DATA = AddMultiProbeGUI(DATA)

if ~isfield(DATA,'toplevel') %no GUI 
    return;
end

it = findobj(DATA.toplevel,'Tag','ProbeId');
set(it,'string',DATA.probenames,'value',DATA.probe);

if ~isempty(findobj(DATA.toplevel,'Tag','MultiProbeMenu'))
    return;
end
DATA.show.ed = 1;
if 0
cw = DATA.plot.cw;
    ch = DATA.plot.ch;

    SPACE = cw;

bp(1) = SPACE;
bp(2) = DATA.gui.toprow+ch+SPACE;
bp(3) = cw.*10;
bp(4) = ch * 1.5;
uicontrol(DATA.toplevel,'Style', 'pushbutton', 'Callback', @CombineAll, ...
'String', 'Combine All', 'Tag','CombineAll','Position', bp);
bp(1) = bp(1)+bp(3)+SPACE;
uicontrol(DATA.toplevel,'Style', 'pushbutton', 'Callback', @ReCombineAll, ...
'String', 'ReCombine All', 'Tag','ReCombineAll','Position', bp);
end

    it = findobj(DATA.toplevel,'Tag','ProbeId');
    bp = get(it,'Position');
    
    if isfield(DATA.probes,'traces')
        np = max([DATA.probes.traces]);
        DATA.subprobes = np;
        if np > 1
            bp(1) = bp(1)+bp(3);
            uicontrol(gcf,'Style', 'pop','String',num2str([0:np]'),'Units', 'norm','Position', bp,...
                'Tag','SubprobeId','Callback',@SetProbeHit,'value',1);
%Add extra measures needed for subprobes
            [DATA.spkvarnames, DATA.spkvarorder] = GetSpikeVals(DATA,NaN,NaN,NaN,[]);
        end
    end
    bp(1) = bp(1)+bp(3);
      uicontrol(gcf,'style','pushbutton','string','+','units','norm','Position',bp, 'Callback',{@SetProbeHit,'next'});

hm = uimenu(DATA.toplevel,'Label','Nchannel','Tag','MultiProbeMenu');
uimenu(hm,'Label','&Load All Spks','Callback',{@LoadAllProbeSpikes, 'allselect'});
uimenu(hm,'Label','&Load Selected Spks','Callback',{@LoadAllProbeSpikes,'select'});
uimenu(hm,'Label','&Load All Cluster Params','Callback',@LoadAllProbeParams);
uimenu(hm,'Label','&Load Spike Times only','Callback',@LoadSpikeTimes);
uimenu(hm,'Label','&Plot All Spks','Callback',@LoadAllProbes);
uimenu(hm,'Label','&RePlot All Spks','Callback',@PlotAllProbeXY);
uimenu(hm,'Label','&Combine All','Callback',@CombineAll);
uimenu(hm,'Label','&Combine All+LFP','Callback',{@CombineAll, 'lfp'});
uimenu(hm,'Label','&Recombine All','Callback',@ReCombineAll);
uimenu(hm,'Label','&Recombine All+','Callback',{@ReCombineAll 'mkall'});
uimenu(hm,'Label','&Recombine All LFPs','Callback',{@ReCombineAll 'lfponly'});
uimenu(hm,'Label','&Recombine All spks','Callback',{@ReCombineAll 'spkonly'});
uimenu(hm,'Label','&Recombine One probe','Callback',{@ReCombineAll 'oneprobe'});
uimenu(hm,'Label','&Recombine Check List','Callback',{@ReCombineAll 'listonly'});
uimenu(hm,'Label','&Spike Shape','Callback', @CalcSpikeShapes);
uimenu(hm,'Label','MeanSpike this probe','Callback', {@combine 'meanspike'});
uimenu(hm,'Label','->One Probe','Callback', {@SetProbeHit 'ReloadProbe'});
uimenu(hm,'Label','->CellList','Callback', @MakeCellList);
uimenu(hm,'Label','CrossCorr','Callback', @xcorrhit);
hm = uimenu(DATA.toplevel,'Label','Plot','Tag','PlotMenu');
uimenu(hm,'Label','All Expts - spkxy','Callback',@PlotAllExpts);
uimenu(hm,'Label','All Expts and probes','Callback',@PlotAllExptsAndProbes);
uimenu(hm,'Label','Plot Subspace','Callback',@PlotSubspace);
uimenu(hm,'Label','LFP Plot','Callback',{@combine, 'maklfp'});
    


function MakeCellList(a,b, varargin)
    DATA = GetDataFromFig(a);
    
wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
ch = DATA.plot.ch;
cw = DATA.plot.cw;
tag = DATA.tag.celllist;

cntrl_box = findobj('Tag',tag,'Name','Cell List');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end
if ~isfield(DATA,'figpos') | isempty(DATA.figpos{1})
   bp = get(DATA.toplevel,'Position');    
   DATA.figpos{1} =  [bp(1)+bp(3) bp(2)-ch*50 cw*40 ch*5]
end
dat.parentfigtag = DATA.tag.top;
% standard figure. Will Plot Map
cntrl_box = figure('Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Cell List','UserData',dat);

bp = [10 10 cw*4 ch*2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', {@AddCellToList, 'Tag' DATA.tag.top},...
'String', 'Add', 'Position', bp);

bbp = [10 bp(2)+ch*2 cw*4 ch*2];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @NextList,...
'String', 'Next', 'Position', bbp);
bbp(2) = bbp(2)+ch*2;

uicontrol(gcf,'Style', 'pushbutton', 'Callback', @MakeCellTemplate,...
'String', 'Tmplt', 'Position', bbp);

bbp(2) = bbp(2)+ch*2;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @SaveCellList,...
'String', 'Save', 'Position', bbp);

bbp(2) = bbp(2)+ch*2;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @AutoFillCellList,...
'String', 'Auto', 'Position', bbp);


np = length(DATA.probelist);
for j = 1:np
    pstrs{j} = num2str(DATA.probelist(j));
    pstrs{j+np} = sprintf('%d/2',DATA.probelist(j));
end
bp(1) = bp(1)+bp(3)+SPACE;
cellx = bp(1);

bp(1) = cellx;
bp(3) = cw*5;
    uicontrol(gcf,'Style', 'Text','String', 'Quality','Position', bp);
    for j = 1:DATA.state.listlen
        bp(1) = bp(1)+bp(3)+SPACE;
        uicontrol(gcf,'Style', 'pop',...
            'String', 'Nothing|MU|MU+|Poor|OK|Good|VGood|Excellent|Automatic', 'Tag', sprintf('CellQuality%d',j), 'Position', bp,'value',1);
    end
    
        bp(1) = bp(1)+bp(3)+SPACE;
    uicontrol(gcf,'Style', 'Text','String', 'Plot','Position', bp);
    bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = cw*5;
    uicontrol(gcf,'Style', 'pop',...
        'String', 'Quality|Number|One Cell|Cells|CQuality|Both|dprime|Q-Dprime|CellDprime', 'Tag', 'PlotType', 'Position',bp, 'callback', @RePlotCellList);
bp(1) = bp(1)+bp(3)+SPACE;


bp(1) = cellx;
bp(3) = cw*5;
bp(2) = bp(2)+bp(4)+VSPACE;
    uicontrol(gcf,'Style', 'Text','String', 'Probe','Position', bp);
    for j = 1:DATA.state.listlen
        bp(1) = bp(1)+bp(3)+SPACE;
    uicontrol(gcf,'Style', 'pop',...
        'String', pstrs, 'Tag', sprintf('CellProbeAdd%d',j), 'Position', bp,'value',DATA.probe);
    end
bp(1) = cellx;
bp(2) = bp(2)+bp(4)+VSPACE;
bp(3) = cw*5;
    uicontrol(gcf,'Style', 'Text','String', 'Cell','Position', bp);

    for j = 1:DATA.state.listlen
        bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = cw*5;
        'String', num2str([0:40]), 'Tag', sprintf('CellNumber%d',j), 'Position', bp,'value',1,...
    uicontrol(gcf,'Style', 'pop',...
        'String', '0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25', 'Tag', sprintf('CellNumber%d',j), 'Position', bp,'value',1,...
        'callback',@SetCellNum);
    end
    PlotCellList(DATA);
    set(gca,'position',[0.1 0.2 0.85 0.78])

function SetCellNum(a,b)
cell = get(a,'value');
RePlotCellList(a,b);



function MakeCellTemplate(a,b)
DATA = GetDataFromFig(a); 
cn = 1;
eid = get(DATA.elst,'value');
exps = DATA.expid(eid);
t(1) = DATA.Expts{exps(1)}.Trials(1).Start(1);
t(2) = DATA.Expts{exps(end)}.Trials(end).End(end);
spkt = DATA.AllData.Spikes.times;
ispk = find(DATA.AllData.Spikes.times > t(1) & DATA.AllData.Spikes.times < t(2)...
    & DATA.AllData.Spikes.codes(:,2) == cn);
mnspk = mean(DATA.AllData.Spikes.values(ispk,:));
subplot(2,1,1);
plot(mnspk);
subplot(2,1,2);

it = findobj(get(a,'parent'),'Tag','CellNumber');
if ~isempty(it)
    cell = get(it(1),'value');
else 
    cell = 1;
end
DATA.Templates(cell,:) = mnspk;
tmpl = DATA.AllData.Spikes.values * mnspk';
for j = length(DATA.Expts):-1:1
    for k = length(DATA.Expts{j}.Trials):-1:1
     ispk = find(spkt > DATA.Expts{j}.Trials(k).Start(1) & spkt < DATA.Expts{j}.Trials(k).End(end));
     if length(ispk)
         score(DATA.Expts{j}.Trials(k).Trial) = max(tmpl(ispk));
     end
    end
end
set(DATA.toplevel,'UserData',DATA);

%SaveCellList(DATA);
plot(score(find(score > 0)));

function RePlotCellList(a,b)
if isstruct(a)
    DATA = a;
    if ~isempty(findobj('Tag',DATA.tag.celllist))
        GetFigure(DATA.tag.celllist);
    end
else
    DATA = GetDataFromFig(a);
    itype =  1;
    it = findobj(get(a,'parent'),'Tag','PlotType');
    if ~isempty(it)
        itype = get(it(1),'value');
    end
    DATA.plot.cellplot = itype;
end

if isempty(findobj('Tag',DATA.tag.celllist))
    return;
end
PlotCellList(DATA);
set(DATA.toplevel,'UserData',DATA);


function AddOneCellToList(caller,b,varargin)
DATA = GetDataFromFig(caller); 
p = DATA.probe;
cl = DATA.currentcluster;
dp = DATA.cluster{cl,p}.dprime;
it = findobj(get(caller,'Parent'),'Tag','ClusterQuality');
q = get(it,'value');
if q == 10
if dp > 5
    q = 6;
elseif dp > 3.5
    q = 5;
elseif dp > 2.8
    q = 4;
elseif dp > 2.2
    q = 3;
else
    q = 2;
end
end
it = findobj(get(caller,'Parent'),'Tag','AddOneCellToList');
cell = get(it,'value');
    [a,b] = TrialRange(DATA);
if DATA.state.online
    a = a + DATA.Expts{DATA.currentexpt}.gui.firsttrial-1;
    b = b + DATA.Expts{DATA.currentexpt}.gui.firsttrial-1;
end
DATA.CellList(cell,a:b) = p;
DATA.CellQuality(cell,a:b) = q;
DATA.CellListCluster(cell,a:b) = DATA.currentcluster;
DATA.cellid = cell(1);
set(DATA.toplevel,'UserData',DATA);
%SaveCellList(DATA);
GetFigure(DATA.tag.celllist);
PlotCellList(DATA,'qual');
SaveCellList(DATA,0,'temp');
    

function AutoFillCellList(caller,b, varargin)

DATA = GetDataFromFig(caller);
eid = DATA.exabsid;
trange(1) = DATA.Expts{eid(1)}.Trials(1).Trial;
trange(2) = DATA.Expts{eid(end)}.Trials(end).Trial;
ncells = size(DATA.CellList,1);
for p = 1:length(DATA.probelist)
    ts = [];
    for j = 1:length(eid)
        E = DATA.Expts{eid(j)};
        if iscluster(E.Cluster,1,p) & isfield(E.Cluster{1,p},'quality')
            Q(p,j) = E.Cluster{1,p}.quality;
        else
            Q(p,j) = NaN;
        end
        ts(j,:) = [E.Trials(1).Trial E.Trials(end).Trial];
    end
    for c = 1:size(DATA.CellList,1)
        id = find(DATA.CellList(c,:) == p);
        nc(c) = length(id);
        if nc(c)
            cc(c) = c;
        else
            cc(c) = 0;
        end
    end
    [a,b] = max(nc);
    id = find(~isnan(Q(p,:)));
    qu = mean(Q(p,id));
    did = find(diff(Q(p,id)));
    if qu > 0
        did(length(did)+1) = length(id);
        first  = 1;
        if a > 0
            c = cc(b);
        else qu > 0;
            ncells = ncells+1;
            c = ncells;
        end
        cellid(p) = c;
        for k = 1:length(did)
            qu = mean(Q(p,id(first):id(did(k))));
            trange = ts(id(first),1):ts(id(did(k)),2);
            DATA.CellList(c,trange) = p;
            DATA.CellQuality(c,trange) = qu;
            DATA.CellListCluster(c,trange) = 1;
            first = did(k)+1;
        end
    else
        cellid(p) = 0;
    end
end
set(DATA.toplevel,'UserData',DATA);
RePlotCellList(caller,b);

function [cells, probes] = CellsSelected(DATA)
    fign = findobj('Tag',DATA.tag.celllist,'Type','Figure');

np = length(DATA.probelist);
    for j = 1:DATA.state.listlen
    cluster(j) = 1;
it = findobj(fign,'Tag',sprintf('CellProbeAdd%d',j));
if ~isempty(it)
    probes(j) = get(it(1),'value');
    if probes(j) > np
        probes(j) = probes(j)-np;
        cluster(j) = 2;
    end
end

it = findobj(fign,'Tag',sprintf('CellNumber%d',j));
if ~isempty(it)
    cells(j) = get(it(1),'value')-1;
end
end
id = find(cells > 0);
cells = cells(id);
probes = probes(id);

function AddCellToList(a,b, varargin)

    
DATA = GetDataFromFig(a); 
eid = DATA.exabsid;

itype =  1;
np = length(DATA.probelist);

for j = 1:DATA.state.listlen
    cluster(j) = 1;
it = findobj(get(a,'parent'),'Tag',sprintf('CellProbeAdd%d',j));
if ~isempty(it)
    probe(j) = get(it(1),'value');
    if probe(j) > np
        probe(j) = probe(j)-np;
        cluster(j) = 2;
    end
end
it = findobj(get(a,'parent'),'Tag',sprintf('CellNumber%d',j));
if ~isempty(it)
    cell(j) = get(it(1),'value')-1;
end
it = findobj(get(a,'parent'),'Tag',sprintf('CellQuality%d',j));
if ~isempty(it)
    quality(j) = get(it(1),'value');
end
end

it = findobj(get(a,'parent'),'Tag','PlotType');
if ~isempty(it)
   itype = get(it(1),'value');
end

[a,b] = TrialRange(DATA);


if max(Counts(cell(cell>0))) > 1  %If list same cell twice, and probes are differnt, use templace scores to assign
    [n,c] = Counts(cell);
    cid = find(n > 1 & c > 0);
    for j = 1:length(cid)
        id = find(cell == c(cid(j)));
        SetCellByTemplateScore(DATA, c(cid(j)), probe(id), [a:b],3);
    end
    return;
end
for j = 1:DATA.state.listlen
    if cell(j) > 0
if quality(j) <= 1
DATA.CellList(cell(j),a:b) = 0;
else
DATA.CellList(cell(j),a:b) = probe(j);
end
if quality(j) > 8  %%Use recorded cluster quality
    for e = 1:length(eid)
        c = DATA.Expts{eid(e)}.Trials(1).Trial;
        if eid(e) < length(DATA.Expts)
            d = DATA.Expts{eid(e)+1}.Trials(1).Trial-1; %fill in gaps
        else
            d = DATA.Expts{eid(e)}.Trials(end).Trial; %fill in gaps
        end
        starts(e) = c;
        ends(e) = d;
        C = DATA.Expts{eid(e)}.Cluster{1,probe(j)};
        if isfield(C,'quality') & C.quality > 0;
            Q(e) = C.quality;
        elseif isfield(C,'autocut') & C.autocut == 0
             Q(e) =  10; %tag. Sort later
        else
            Q(e) = 0;
        end
    end
    id = find(Q < 10);
    if length(id) % at least one is set
        last = id(end);
        id = find(Q == 10)
        isset = 0;

        for e = 1:length(eid)
            %if have set qual in this block, and it is set again later, fill in with the
            %last set value
            if Q(e) == 10 & isset & e < last
                Q(e) = Q(e-1);
            else
                isset = 1;
            end
            if Q(e) == 0
                DATA.CellQuality(cell(j),starts(e):ends(e)) =  Q(e);
                DATA.CellList(cell(j),starts(e):ends(e)) =  0;

            elseif Q(e) < 10
                DATA.CellQuality(cell(j),starts(e):ends(e)) =  Q(e);
            end
        end
    else
        fprintf('No Cluster Quality data for Cell%d, Probe %d\n',cell(j),probe(j));
    end
else
    DATA.CellQuality(cell(j),a:b) = quality(j);
end
DATA.CellListCluster(cell(j),a:b) = cluster(j);
    end
end

len = max([length(DATA.CellQuality)  length(DATA.CellList) length(DATA.CellListCluster)]);
if length(DATA.CellQuality) < len
    lq = length(DATA.CellQuality);
    DATA.CellQuality(:,lq+1:len)= 0;
end
if length(DATA.CellList) < len
    lq = length(DATA.CellList);
    DATA.CellList(:,lq+1:len) = 0;
end
if length(DATA.CellListCluster) < len
    DATA.CellListCluster(:,len) = 0;
end

DATA.cellid = cell(1);
set(DATA.toplevel,'UserData',DATA);
%SaveCellList(DATA);
PlotCellList(DATA,'plotshapes');
SaveCellList(DATA,0,'temp');


function [a,b] = TrialRange(DATA)
    eid = get(DATA.elst,'value');
    if isfield(DATA,'expid')
    exps = DATA.expid(eid);
    else
        exps = 1;
    end
    if DATA.firsttrial > 1
        a = DATA.Expts{exps(1)}.Trials(DATA.firsttrial).Trial;
    else
        a = DATA.Expts{exps(1)}.Trials(1).Trial;
    end
    if DATA.lasttrial >= DATA.firsttrial && DATA.lasttrial > 0
        b = DATA.Expts{exps(end)}.Trials(DATA.lasttrial).Trial;
    else
        b = DATA.Expts{exps(end)}.Trials(end).Trial;
    end
    
    
function SaveCellList(a,b, varargin)
    
            
    if isstruct(a)
        DATA = a;
    else
        DATA = GetDataFromFig(a);
    end

    if isfield(DATA,'CellList')
    CellList = DATA.CellList;
    CellQuality = DATA.CellQuality;
    CellListCluster = DATA.CellListCluster;
    else
        CellList = [];
        CellQuality = [];
        CellListCluster = [];
    end
        
    Templates = DATA.Templates;
    TemplateInfo = DATA.TemplateInfo;
    if ~isfield(DATA,'cellfile')
        DATA.cellfile = strrep(DATA.datafilename,'.mat','.cells.mat');
    end

    cellfile = DATA.cellfile;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'temp',4)
            cellfile = strrep(DATA.cellfile,'cells','celltmp');
        end
        j = j+1;
    end
    
    save(cellfile,'CellList','CellQuality', 'Templates','CellListCluster','TemplateInfo');
    set(DATA.toplevel,'UserData',DATA);

function ClearAnnotations(a)
    c = get(a,'Children')
    for j = 1:length(c)
        type = get(c(j),'Type');
        if strcmp(type,'line')
            delete(c(j));
        end
    end
        
function PlotCellList(DATA, varargin)
    
QUALITY = 1;
CELLNUM=2;
ONECELL=3;
JOINCELL=4;
CLUSTERQUALITY=5;
CLUSTERANDCELL=6;
DPRIME=7;
CELLANDDPRIME = 9;
CQDPRIME=8;
replotshape = 0;

cell = 1;

if ~isfield(DATA.plot,'showspikeshapes')
    showspikeshape = 0;
else
    showspikeshape = DATA.plot.showspikeshapes;
end
if isfield(DATA.plot,'cellplot')
    plottype = DATA.plot.cellplot;
else
plottype = QUALITY;
end

n = length(DATA.probelist);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'bynumber',4)
        plottype = CELLNUM;
    elseif strncmpi(varargin{j},'onecell',4)
        plottype = ONECELL;
    elseif strncmpi(varargin{j},'joincell',4)
        plottype = JOINCELL;
    elseif strncmpi(varargin{j},'clusterquality',7)
        plottype = CLUSTERQUALITY;
    elseif strncmpi(varargin{j},'plotshapes',8)
        replotshape = 1;
    end
    j = j+1;
end

if ~isfield(DATA,'CellList') | isempty(DATA.CellList)
    DATA.CellList = zeros(1,length(DATA.AllData.Trialids));
    DATA.CellQuality = zeros(1,length(DATA.AllData.Trialids));
    DATA.CellListCluster = zeros(1,length(DATA.AllData.Trialids));
%    return;
end
it = findobj('Tag','CombinerCellList');
it = findobj(it, 'Tag','CellNumber');
if ~isempty(it)
    cell = get(it(1),'value');
end
cfig = gcf;
if showspikeshape
fign = findobj('Tag','SpikeShapes','Type','Figure');
if (isempty(fign) | replotshape) & exist(DATA.meanspkfile,'file');
    PlotSpikeShapes(DATA.meanspkfile,'cells',DATA.CellList,DATA.CellListCluster);
end
if isfield(DATA, 'exabsid')
PlotSpikeShapes([DATA.exabsid(1) DATA.exabsid(end)]);
end
figure(cfig);
end

ClearAnnotations(gca);
Q = zeros(length(DATA.probelist),DATA.Expts{end}.Trials(end).Trial);
im = zeros(length(DATA.probelist)*2,DATA.Expts{end}.Trials(end).Trial);
dp = zeros(size(Q));
for j = 1:size(DATA.CellList,1)
    id = find(DATA.CellList(j,:) > 0);
    k = DATA.CellList(j,id);
    ind = sub2ind(size(im),k*2,id);
    aid = find(DATA.CellList(j,:) > 0 & DATA.CellListCluster(j,:) <= 1);
    if length(aid)
        k = DATA.CellList(j,aid)*2;
        indb = sub2ind(size(im),k-1,aid);
    else
        indb = [];
    end
    bid = find(DATA.CellList(j,:) > 0 & DATA.CellListCluster(j,:) == 2 & DATA.CellQuality(j,:) > 0);
    if length(bid)
        k = DATA.CellList(j,bid)*2;
        indc = sub2ind(size(im),k,bid);
    else
        indc = [];
    end
    indb = [indb indc];

    if isfield(DATA,'CellQuality')
        im(ind) = DATA.CellQuality(j,id);
        im(indb) = DATA.CellQuality(j,[aid bid]);
    else
        im(ind) = j;
        im(indb) = j;
    end
end
    np = length(DATA.probelist);
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        ap = get(gca,'position');

    for e = 1:length(DATA.Expts)
        id = [DATA.Expts{e}.Trials.Trial];
        if isfield(DATA.Expts{e}.Trials,'ed')
            eds(id) = [DATA.Expts{e}.Trials.ed];
        else
            eds(id) = DATA.Expts{e}.Stimvals.ed;
        end
        if isfield(DATA.Expts{e},'Cluster') && size(DATA.Expts{e}.Cluster,1) > 0
            for k = 1:min([size(DATA.Expts{e}.Cluster,2) np])
                if isfield(DATA.Expts{e}.Cluster{1,k},'dprime')
                    dp(k,id) = DATA.Expts{e}.Cluster{1,k}.dprime;
                end
                if isfield(DATA.Expts{e}.Cluster{1,k},'quality')
                    Q(k,id) = DATA.Expts{e}.Cluster{1,k}.quality;
                end
            end
        end
        t = [DATA.Expts{e}.Trials(1).Trial DATA.Expts{e}.Trials(end).Trial];
        if t(1) >= xl(1) && t(1) <= xl(2)
            xp = ap(1)+ap(3) * (t(1)-xl(1))/range(xl);
            hlines(j) = annotation('line',[xp xp],[ap(2)+ap(4) 1.0]);
            set(hlines(j),'color','r');
        end
    end

if plottype == QUALITY || plottype == CLUSTERANDCELL
    if plottype == CELLANDDPRIME
      for j = 1:np
          im((j*2)-1,:) = dp(j,:);
      end
    end
    if plottype == CLUSTERANDCELL
      for j = 1:np
          im((j*2)-1,:) = Q(j,:);
      end
    end
    hold off;
    imagesc([1:size(im,2)],[0.75 np+0.25],im);
end
if plottype ==CLUSTERQUALITY
    if plottype == CLUSTERANDCELL
        subplot(2,1,2);
    else
        subplot(1,1,1);
    end
    Q = zeros(length(DATA.probelist),DATA.Expts{end}.Trials(end).Trial);
    for e = 1:length(DATA.Expts)
        id = [DATA.Expts{e}.Trials.Trial];
        if isfield(DATA.Expts{e}.Trials,'ed')
            eds(id) = [DATA.Expts{e}.Trials.ed];
        else
            eds(id) = DATA.Expts{e}.Stimvals.ed;
        end
        if isfield(DATA.Expts{e},'Cluster') && size(DATA.Expts{e}.Cluster,1) > 0
            for k = 1:size(DATA.Expts{e}.Cluster,2)
                if isfield(DATA.Expts{e}.Cluster{1,k},'quality')
                    Q(k,id) = DATA.Expts{e}.Cluster{1,k}.quality;
                end
            end
        end
    end
    hold off;
    imagesc(Q);
    caxis([0 8]);
    hold on;
    cl = (DATA.CellListCluster-1.5)/2;

    plot((DATA.CellList')+cl','linewidth',2);
    legend(num2str([1:size(DATA.CellList,1)]'),'location','NorthWest');
    set(gca,'YaxisLocation','right');
    [a,b] = TrialRange(DATA);
    plot([a a],[1 length(DATA.probelist)],'w:');
    plot([b b],[1 length(DATA.probelist)],'w:');
elseif plottype == CELLANDDPRIME
      for j = 1:np
          im((j*2)-1,:) = dp(j,:);
      end
    hold off;
    imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
    caxis([0 7]);
elseif plottype == CQDPRIME
      for j = 1:np
          im((j*2)-1,:) = dp(j,:);
          im((j*2),:) = Q(j,:);
      end
    hold off;
    imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
    caxis([0 7]);
elseif plottype ==DPRIME
      for j = 1:np
          im((j*2)-1,:) = dp(j,:);
          im((j*2),:) = dp(j,:);
      end
    hold off;
    imagesc([1:size(DATA.CellList,2)],[0.75 np+0.25],im);
    caxis([0 7]);
elseif plottype == CELLNUM
    DATA.CellList(find(DATA.CellList == 0)) = NaN;
    hold off;
    plot(DATA.CellList','linewidth',2);
    set(gca,'ylim',[min(DATA.CellList(:))-1 max(DATA.CellList(:))+1]);
    legend(num2str([1:size(DATA.CellList,1)]'),'location','SouthWest')
elseif plottype == ONECELL
    trials = find(DATA.CellList(cell,:) > 0);
    hold off;
    plot(trials,DATA.CellList(cell,trials),'o-');
elseif plottype == JOINCELL
    hold off;
    colors = mycolors;
    mid = (size(DATA.CellList,1)-1)/2;
    for j = 1:size(DATA.CellList,1)
    trials = find(DATA.CellList(j,:) > 0);
    plot(trials,DATA.CellList(j,trials)+0.05*(j-mid),'o-','color',colors{j});
    hold on;
    end
    [a,b] = TrialRange(DATA);
    plot([a a],[1 length(DATA.probelist)],'k:');
    plot([b b],[1 length(DATA.probelist)],'k:');
        
    legend(num2str([1:size(DATA.CellList,1)]'),'location','SouthWest')
    hold off;
    set(gca,'ylim',[min(DATA.CellList(:))-1 max(DATA.CellList(:))+1]);
end

if ~ismember(plottype,[CELLNUM])
        hold on;
    DATA.CellList(find(DATA.CellList == 0)) = NaN;
    cl = (DATA.CellListCluster-0.5)/2;

    h = plot((DATA.CellList')+cl','linewidth',2);
    for j = 1:length(h)
        colors{j} = get(h(j),'color');
    end
    legend(num2str([1:size(DATA.CellList,1)]'),'location',[0 0.5 0.08 0.5])
     set(gca,'YaxisLocation','right');
    [a,b] = TrialRange(DATA);
    plot([a a],[1 length(DATA.probelist)],'w:');
    plot([b b],[1 length(DATA.probelist)],'w:');
    id = find(eds > 0);
%plot depth with sign inverte, so that it matches predicted movement
%of cells on probe. Ingreasing depth = cells move to lowe probe #s
    if range(eds(id)) > 8 .* 0.15;
        aid = find(eds(a:b) >0) + a -1;
        a = mean(eds(aid));
        if a-min(eds(id)) < 0.65
            a = min(eds(id)+0.65);
        end
%add 0.6 because y axis starts at 0.5;
        plot(id, 0.6+((a-eds(id))./0.15),'w');
    else
        plot(id, (max(eds(id))-eds(id))./0.15,'w');
    end
    for j = 1:size(DATA.CellList,1)
        id = find(DATA.CellList(j,:) > 0);
        if length(id)
        text(id(1),DATA.CellList(j,id(1))+0.1,num2str(j),'color',colors{j});
        did = find(abs(diff(DATA.CellList(j,id)))>0)+1;
        for k = 1:length(did)
            text(id(did(1)),DATA.CellList(j,id(did(1)))+0.1,num2str(j),'color',colors{j});
        end
        end
    end
end

function AddProbeList(DATA)
% Add A list of probes for multiple probe recording.
%
ch = DATA.plot.ch;
cw = DATA.plot.cw;
wsiz = [DATA.plot.cw*40,DATA.plot.ch*40];

ival = 1;

it =  findobj(DATA.toplevel,'Tag','UseProbe1');
if ~isempty(it)
%    delete(it); %for debugging
    return;
end
%
% if probes were selected on the command line, useprobe will exist.
% Otherwise not.
if isfield(DATA.plot,'useprobe')
    up = DATA.plot.useprobe;
else
    up = zeros(size(DATA.probelist));
end
for j = 1:length(DATA.probelist)
    bp = [0.94 (j-1)./(length(DATA.probelist)+1) 0.07 0.05];
    uicontrol(DATA.toplevel,'Style', 'checkbox',...
        'String', num2str(j), 'Tag', sprintf('UseProbe%d',j),'units','norm', 'Position', bp,'value',up(j),...
        'Callback',@Update);
end
    bp = [0.94 (j)./(length(DATA.probelist)+1) 0.06 0.04];
    uicontrol(DATA.toplevel,'Style', 'checkbox',...
        'String', 'sync', 'Tag', 'ShowSync', 'units','norm','Position', bp,'value',0,...
        'Callback',@Update);
ws = get(DATA.elst,'position');
ws(3) = 0.92;
set(DATA.elst,'position',ws);
ws = get(DATA.clst,'position');
ws(3) = 0.92;
set(DATA.clst,'position',ws);
ws = get(DATA.toplevel,'position');
ws(3) = ws(3) + cw * 4;
set(DATA.toplevel,'position',ws);


function GuiResize(a,b);
        
  DATA = GetDataFromFig(a);
  sz = get(a,'position');
  if isfield(DATA,'clst')
  end
  if isfield(DATA,'elst')
  end

function DATA = BuildGUI(DATA)

scrsz = get(0,'Screensize');
DATA.gui.scrsz = scrsz;
cw= scrsz(3)/140;
ch= scrsz(4)/60;
if scrsz(3) > 2000
    cw = 10;
elseif scrsz(3) > 1900
    cw = 10;
end
if scrsz(4) > 1200
    ch = 10;
elseif scrsz(4) == 1200
    ch = 8;
elseif scrsz(4) > 1100
    ch = 9;
else
    ch= 9;
end
if DATA.state.bigwindow(1) > 0
    cw = cw * DATA.state.bigwindow(1);
end
if DATA.state.bigwindow(2) > 0
    ch = ch * DATA.state.bigwindow(2);
end
DATA.plot.cw = cw;
DATA.plot.ch = ch;
DATA.user = getenv('username');
if isempty(DATA.user) %for MAC/Linux
DATA.user = getenv('username');
end
%fprintf('User is %s\n',DATA.user);
wsiz = [cw*40,ch*44];
SPACE = 0.01;
nch = 0.05;
ncw=0.022;

cntrl_box = figure('Position', [100 scrsz(4)-(wsiz(2)+ch*4) wsiz(1) wsiz(2)],...
    'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.tag.top,'ResizeFcn',@GuiResize);
bp = [0.02 0.02 0.98 0.25];
DATA.elst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['combine(''setexpt'',''Tag'',''' DATA.tag.top ''');'],'Tag','subexptlist',...
    'Units','norm','Position',bp,'Max',3,'Min',1);
bp = [0.02 bp(2)+bp(4)+0.01 ncw*6 nch];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['combine(''combine'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'Combine','Units','norm', 'Position', bp);

bp(1) = bp(1)+bp(3)+SPACE;
bp(3) = ncw * 4;
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['combine(''Save'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'Save', 'Units','norm','Position', bp);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*7 nch];
uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['combine(''SaveLFP'',''Tag'',''' DATA.tag.top ''');'],...
'String', 'SaveLFP', 'Units', 'norm', 'Position', bp);


bp = [SPACE bp(2)+bp(4)+SPACE ncw*30 nch];
DATA.saveitem = uicontrol(gcf,'Style','Edit','String','Save As','Units', 'norm','Position',bp);

bp = [0.01 bp(2)+bp(4)+SPACE 0.98 nch*6];
DATA.clst = uicontrol(gcf, 'Style','listbox',...
    'Callback', ['combine(''listexps'',''Tag'',''' DATA.tag.top ''');'],'Tag','explist',...
    'Units', 'norm','Position',bp,'Max',3,'Min',1);

bp = [0.01 bp(2)+bp(4)+SPACE ncw*9 nch];
uicontrol(gcf,'Style','Text','String','Data File','Units', 'norm','Position',bp);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*30 nch];
uicontrol(gcf,'Style','Edit','String',DATA.datafilename,'Units', 'norm','Position',bp,'Callback',['combine(''newfile''''Tag'',''' DATA.tag.top ''');'],...
    'Tag','FileName');

bp = [0.01 bp(2)+bp(4)+SPACE ncw*6 nch];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'ShowN', 'Tag', 'ShowN', 'Units', 'norm','Position', bp, 'value',DATA.plot.showN);

bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Spikes', 'Tag', 'ShowSpikes', 'Units', 'norm','Position', bp,'value',DATA.state.showspikes,'Callback',@Update);
bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'SpkXY', 'Tag', 'SpkXY', 'Units', 'norm','Position', bp,'value',DATA.state.showspkxy,...
'Callback',@Update);

bp(1) = bp(1) + bp(3) + SPACE;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Clean', 'Tag', 'ResetClusters', 'Units', 'norm','Position', bp,'value',DATA.state.resetclusters,...
'Callback',@Update);

bp(1) = bp(1) + bp(3) + SPACE;
bp(3) = ncw*7;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Auto', 'Tag', 'AutoPlot','Units', 'norm', 'Position', bp,'value',DATA.state.autoplot,...
'Callback',@Update);

bp = [0.01 bp(2)+bp(4)+SPACE ncw*6 nch];
uicontrol(gcf,'Style', 'pop',...
'String', 'Mean|Seq:time|Seq:Trials|Seq:Id|Seq:Edit|Psych Only|Collapse Y|Cp Trials|CP Hist|sdf|ACloop|ACresp|sdf for 0|Sdf collapse X|Sdf collapse Y|Counts|Subspace', ...
'Tag', 'PlotSeq','Units', 'norm', 'Position', bp,'value',DATA.state.plotseq+1,...
'Callback',@Update);

bp = [bp(1)+bp(3)+SPACE bp(2) ncw*6 nch];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Recount', 'Tag', 'Recount', 'Units', 'norm','Position', bp,'value',DATA.state.recount,...
'Callback',@Update);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*5 nch];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Psych', 'Tag', 'PlotPsych', 'Units', 'norm','Position', bp,'value',DATA.state.plotpsych,...
'Callback',@Update);
bp = [bp(1)+bp(3)+SPACE bp(2) ncw*8 nch];
uicontrol(gcf,'Style', 'checkbox',...
'String', 'Combined', 'Tag', 'PlotCombined', 'Units', 'norm','Position', bp,'value',DATA.state.plotcombined,...
'Callback',@Update);



bp = [0.01 bp(2)+bp(4)+SPACE ncw*3 nch];
uicontrol(gcf,'Style', 'checkbox',...
'String', '0', 'Tag', 'UseCluster0', 'Units', 'norm','Position', bp,'Callback',{@SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '1', 'Tag', 'UseCluster1', 'Units', 'norm','Position', bp,'Callback',{@SetClusters, DATA.tag.top});

bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '2', 'Tag', 'UseCluster2', 'Units', 'norm','Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '3', 'Tag', 'UseCluster3', 'Units', 'norm','Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'checkbox',...
'String', '4', 'Tag', 'UseCluster4', 'Units', 'norm','Position', bp,'Callback',{@SetClusters, DATA.tag.top});
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', @FitButton, ...
'String', 'Fit', 'Tag','FitButton','Units', 'norm','Position', bp);
bp(1) = bp(1)+bp(3);
uicontrol(gcf,'Style', 'pushbutton', 'Callback', 'combine(''setexpplot'')', ...
'String', 'Plot', 'Tag','PlotButton','Units', 'norm','Position', bp);
bp(1) = bp(1)+bp(3);
bp(3) = ncw*5;
uicontrol(gcf,'Style', 'checkbox',...
'String', 'LFP', 'Tag', 'UseLFP', 'Units', 'norm','Position', bp,'Callback',@Update);

    bp(1) = bp(1)+bp(3);
h = uicontrol(gcf,'Style', 'pop','String',DATA.probenames,'Units', 'norm','Position', bp,...
      'Tag','ProbeId','Callback',@SetProbeHit,'value',1);
if length(DATA.probelist) > 1
    bp(1) = bp(1)+bp(3);
      uicontrol(gcf,'style','pushbutton','string','+','Position',cp, 'Callback',{@SetProbeHit,'next'});
end

bp = [0.01 bp(2)+bp(4)+SPACE 0.98 nch];
uicontrol(gcf,'Style','Edit','String','Comments','Units', 'norm','Position',bp,'callback',@AddComment);

DATA.gui.toprow = bp(2);

  hm = uimenu(gcf,'Label','Mark');
  uimenu(hm,'Label','&Mark','Callback','combine(''Mark'')');
  uimenu(hm,'Label','&relist','Callback','combine(''relist'')');
  uimenu(hm,'Label','&Options','Callback','combine(''options'')');
  uimenu(hm,'Label','&ShowVals','Callback','combine(''showvals'')');
  uimenu(hm,'Label','&To Front','Callback','combine(''winfront'')');
  uimenu(hm,'Label','&Update RF','Callback','combine(''rfupdate'')');
  uimenu(hm,'Label','&Next CEell','Callback','combine(''nextcell'')');
  uimenu(hm,'Label','&Update Lists','Callback','combine(''checklists'')');
  uimenu(hm,'Label','&Comments','Callback','combine(''Comments'')');
  uimenu(hm,'Label','Save Comments','Callback','combine(''SaveComments'')');
  uimenu(hm,'Label','&Close','Callback',['combine(''Close'',''Tag'',''' DATA.tag.top ''');']);
  cm = uimenu(gcf,'Label','&Cluster');
  uimenu(cm,'Label','&AutoCut','Callback','combine(''autocut'')');
  uimenu(cm,'Label','&AutoCluster','Callback','combine(''autocluster'')');
  uimenu(cm,'Label','&AutoCluster All','Callback','combine(''allautocluster'')');
  uimenu(cm,'Label','&TrackCluster','Callback','combine(''trackcluster'')');
  uimenu(cm,'Label','&DDF','Callback','combine(''plotddf'')');
  uimenu(cm,'Label','&FixMains','Callback','combine(''fixmains'')');
  uimenu(cm,'Label','&ClearAll (one probe)','Callback','combine(''clearallclusters'')');
  uimenu(cm,'Label','&ClearOnline','Callback','combine(''clearonlineclusters'')');
  uimenu(cm,'Label','&Reload All','Callback',{@ReloadClusters, 'all'});
  uimenu(cm,'Label','&Reload One Probe','Callback',{@ReloadClusters});
  uimenu(cm,'Label','&Reload and Classify','Callback',{@ReloadClusters,'relcassify'});
  uimenu(cm,'Label','&Make Templates','Callback',@MakeTemplates);
  uimenu(cm,'Label','&Add Templates','Callback',@AddTemplate);
  uimenu(cm,'Label','&Plot Templates','Callback',@PlotTemplates);
  if DATA.state.online
  cm = uimenu(gcf,'Label','&Relist','Callback','combine(''relist'')');
  end
  cm = uimenu(gcf,'Label','&Options','Callback','combine(''options'')');
  set(gcf,'Menubar','none');
DATA.toplevel = cntrl_box;

function MakeTemplates(a,b, varargin)
DATA = GetDataFromFig(a);
eid = DATA.currentexpt;
if isfield(DATA.Expts{eid}.gui,'spks')
ispk = DATA.Expts{eid}.gui.spks;
else
end
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
Spks = DATA.AllData.Spikes;
end
nc = length(unique(Spks.codes(ispk,2)));
for j = 1:nc
    id = find(Spks.codes(ispk,2) == j-1);
    Template(j,:) = mean(Spks.values(ispk(id),:));
    DATA.TemplateInfo(j).exid = eid;
    DATA.TemplateInfo(j).probe = DATA.probe;
    DATA.TemplateInfo(j).cluster = DATA.currentcluster-1;
    DATA.TemplateInfo(j).cellid = 1;
end
DATA.Templates = Template;
GetFigure('CombinerSpikeV');
plot(Template');
SaveCellList(DATA);

set(DATA.toplevel,'UserData',DATA);
NotBusy(DATA);

function PlotTemplates(a,b, varargin)
DATA = GetDataFromFig(a);
GetFigure('CombinerSpikeV');
plot(DATA.Templates');

function AddTemplate(a,b, varargin)
DATA = GetDataFromFig(a);
eid = DATA.currentexpt;
if isfield(DATA.Expts{eid}.gui,'spks')
ispk = DATA.Expts{eid}.gui.spks;
else
end
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
Spks = DATA.AllData.Spikes;
end
id = find(Spks.codes(ispk,2) == DATA.currentcluster);
Template(1,:) = mean(Spks.values(ispk(id),:));
DATA.Templates = cat(1, DATA.Templates, Template);
GetFigure('CombinerSpikeV');
plot(DATA.Templates');
nt = size(DATA.Templates,1);
DATA.TemplateInfo(nt).exid = eid;
DATA.TemplateInfo(nt).probe = DATA.probe;
DATA.TemplateInfo(nt).cluster = DATA.currentcluster;
DATA.TemplateInfo(nt).cellid = nt;
SaveCellList(DATA);
set(DATA.toplevel,'UserData',DATA);


function ReloadClusters(a,b, varargin)
DATA = GetDataFromFig(a);
probes = DATA.probe;
probe = DATA.probe;
reclassify = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'all',3)
        probes = DATA.probelist;
    elseif strncmpi(varargin{j},'relcassify',5)
        reclassify = 1;
    end
    j = j+1;
end
for j = 1:length(probes)
    DATA.probe = probes(j);
    DATA = LoadClusters(DATA,ClusterFile(DATA));
    if reclassify
        DATA = ReClassifyAll(DATA,'probes',probes(j));
        SetExptClusters(DATA);
    end
end
DATA.probe = probe;
set(DATA.toplevel,'UserData',DATA);
        
function SetCellByTemplateScore(DATA, cell, probes, trials, tn)

    scores = squeeze(DATA.TemplateScores(probes, tn, trials));
    sid = find(sum(scores > 0) == length(probes)); %trials where all probes have scores
    [ts, id] = max(scores(:,sid));
    probe(sid) = probes(id);
    id = find(probe == 0);
    ip = interp1(sid,probe(sid),id);
    probe(id) = round(ip);
    GetFigure('TemplateScore');
    hold off;
    imagesc(trials,[1:size(DATA.TemplateScores,1)],squeeze(DATA.TemplateScores(:,tn,trials)));
    hold on;
    plot(trials,probe,'w-');
    GetFigure(DATA.tag.celllist);
    plot(trials,probe,'w-');
    a = questdlg('Apply?','popup','Cancel','OK','OK');
    if strcmp(a,'OK')
        DATA.CellList(cell,trials) = probe;
% not clear how to set quatlity yet.  So leave as defined for this cell
%        DATA.CellQuality(cell,trials) = probe;
        set(DATA.toplevel,'UserData',DATA);
    end

    
    
    
function CalcSpikeShapes(a,b, varargin)
%
% caluclate mean spike waveform for biggest spikes
% for each expt,probe
DATA = GetDataFromFig(a);
rebuild = DATA.state.forcebuild;

probelist = DATA.probelist;
TrialVar.cx = [];
outname = strrep(DATA.datafilename,'.mat','spks.mat');
j = 1;
while j < length(varargin)
    if strncmpi(varargin{j},'online',5)
        outname = strrep(DATA.datafilename,'.mat','spks.mat');
    end
end
if exist(outname,'file') & ~rebuild
    load(outname);
    DATA.MeanSpike = MeanSpike;
    DATA.TrialVar = TrialVar;
    DATA.TemplateScores = TemplateScores;
    if exist('TemplateInfo','var')
        DATA.TemplateInfo = TemplateInfo;
        clear TemplateInfo;
    end
else
    if ~isfield(DATA,'TemlateScores')
        DATA.TemplateScores = [];
    end
    oldstate = DATA.state.autoplotnewprobe;
    DATA.state.autoplotnewprobe = 1;
    for p = 1:length(probelist)
        if p ~= DATA.probe
            DATA = SetProbe(DATA, probelist(p));
        end
        DATA = CalcMeanSpike(DATA,1:length(DATA.Expts));
    end
    DATA.state.autoplotnewprobe = oldstate;
    DATA = SaveSpikeShape(DATA,outname);
end

GetFigure('SpikeShape');
subplot(2,1,1);
imagesc(std(DATA.MeanSpike.v,[],3)');
subplot(2,1,2);
id = find(sum(DATA.TrialVar.cx) > 0);
imagesc(abs(DATA.TrialVar.cx(:,id)));
set(DATA.toplevel,'UserData',DATA);
GetFigure('ClusterShape');
PlotSpikeShapes(DATA.MeanSpike,'auto');
GetFigure('TemplateScore')
PlotTemplateScores(DATA,2);


function PlotTemplateScores(DATA,ti,varargin)

    trials = [1:size(DATA.TemplateScores,3)]';
    trialrange = [1 length(trials)];
    j = 1;
  
  while j <= length(varargin)
      if strncmpi(varargin{j},'Trials',5)
          j = j+1;
          trialrange = varargin{j};
      end
      j = j+1;
  end
  if isfield(DATA.Comments,'Peninfo')
      probesep = DATA.Comments.Peninfo.probesep;
  else
      probesep = 75;
  end
  nprobes = size(DATA.TemplateScores,1);
  id = find(squeeze(sum(DATA.TemplateScores(:,ti,:),1)) > 0 & trials > trialrange(1));
  imagesc(trials(id),[1 nprobes],squeeze(DATA.TemplateScores(:,ti,id)));
  Trials = [];
  eds = []
  Expts = []; nx = 1;
  for ex = 1:length(DATA.Expts)
      [a,b] = ismember(DATA.Expts{ex}.Trials(1).Trial,id);
      if a
          Expts(nx).Start = id(b); %id # of starting trial
          nx = nx+1;
      end
      Trials = [Trials [DATA.Expts{ex}.Trials.Trial]];
      if isfield(DATA.Expts{ex}.Trials,'ed')
      eds([DATA.Expts{ex}.Trials.Trial]) = [DATA.Expts{ex}.Trials.ed];
      else
      eds([DATA.Expts{ex}.Trials.Trial]) = DATA.Expts{ex}.Stimvals.ed;
      end
  end
  hold on;
  if size(DATA.TemplateInfo,1) >= ti
  estart = GetEval(DATA.Expts{DATA.TemplateInfo(ti).exid},'ed');
  ed = DATA.TemplateInfo(ti).probe + (eds(id)-estart) .* 1000./probesep;
  plot(trials(id),ed,'w');
  for j = 1:length(Expts)
      plot([Expts(j).Start Expts(j).Start],[1 nprobes],':');
  end
  end
  
    
function DATA = SaveSpikeShape(DATA, outname)
for j = 1:length(DATA.Expts)
    DATA.MeanSpike.eds(j) = DATA.Expts{j}.Stimvals.ed;
    DATA.MeanSpike.Trials(j,:) = [DATA.Expts{j}.Trials([1 end]).Trial];
    if isfield(DATA,'Header')
        DATA.MeanSpike.Header = DATA.Header;
    else
        DATA.MeanSpike.probesep = 150; %default
    end
end
MeanSpike = DATA.MeanSpike;
TrialVar = DATA.TrialVar;
TemplateScores = DATA.TemplateScores;
TemplateInfo = DATA.TemplateInfo;
Templates = DATA.Templates;
fprintf('Saving MeanSpike to %s\n',outname);

save(outname,'MeanSpike','TrialVar','TemplateScores','Templates','TemplateInfo');
        
        
function AddComment(a,b)
    DATA = GetDataFromFig(a);
    n = 1;
    if isfield(DATA.Comments,'Offline')
        n = length(DATA.Comments.Offline)+1;
    end
    cfile = strrep(DATA.datafilename,'.mat','.txt');
    oid = fopen(cfile,'a');
    DATA.Comments.Offline{n} = get(a,'string');
    set(a,'string','');
    fprintf(oid,'%s\n',DATA.Comments.Offline{n});
    fclose(oid);
    set(DATA.toplevel,'UserData',DATA);
    
    
        
function DATA = CalcMeanSpike(DATA,expid)

    p = DATA.probe;
    Spks = GetSpikeData(DATA, p);

    if ~isfield(DATA,'TemplateScores')
        DATA.TemplateScores = [];
    end
    if size(Spks.values,2) == size(DATA.Templates,2)
    for j = 1:size(DATA.Templates,1)
        scores(j,:) = Spks.values * DATA.Templates(j,:)';
    end
    else
        scores = [];
        DATA.TemplateScores = [];
    end
    if size(scores,1) > size(DATA.TemplateScores,2) %% new templates
        DATA.TemplateScores(:,size(scores,1),:) = 0;
    end
    empties = [];  %trial with no spikes, not missing trials
    for eid = expid;
        [DATA, ispk] = SetExptSpikes(DATA, eid,0);
        if length(ispk) > 100
            if length(ispk) > 10000
                xid = find(DATA.Spikes.cx(ispk) > prctile(DATA.Spikes.cx(ispk),99.9));
            else
                v = sort(DATA.Spikes.cx(ispk),'descend');
                xid = find(DATA.Spikes.cx(ispk) > v(100));
            end
            DATA.MeanSpike.v(eid,p,:) = mean(Spks.values(ispk(xid),:));
            DATA.MeanSpike.sd(eid,p,:) = std(Spks.values(ispk(xid),:));
        end
        trials = [DATA.Expts{eid}.Trials.Trial];
        for j = 1:length(DATA.Expts{eid}.Trials)
            times = [DATA.Expts{eid}.Trials(j).Start(1) DATA.Expts{eid}.Trials(j).End(end)];
            tspk = find(Spks.times(ispk) > times(1) & ...
                Spks.times(ispk) < times(2));
            cspk = find(Spks.codes(ispk(tspk),2) > 0);
            if isempty(cspk)
                sv = sort(DATA.Spikes.cx(tspk));
                if length(sv) > 5
                    DATA.TrialVar.cx(p,trials(j)) = -mean(sv(end-5:end));
                else
                    DATA.TrialVar.cx(p,trials(j)) = -mean(DATA.Spikes.cx(tspk));
                end
            else
                cspk = ispk(tspk(cspk));
                DATA.TrialVar.cx(p,trials(j)) = mean(DATA.Spikes.cx(cspk));
            end
            if length(tspk) && ~isempty(scores)
                DATA.TemplateScores(p,:,trials(j)) = max(scores(:,ispk(tspk)),[],2);
            else
               empties = [empties, trials(j)];
            end
        end
        trigt = 9;
        if isfield(DATA.Expts{eid},'Cluster') & iscluster(DATA.Expts{eid}.Cluster,1,p) ...
            & isfield(DATA.Expts{eid}.Cluster{1,p},'dprime')
            DATA.MeanSpike.dprimes(p,eid) = DATA.Expts{eid}.Cluster{1,p}.dprime;
            spks = find(Spks.codes(ispk,2) == 1);
%only average the group that triggered in the same direction
            sgns = sign(Spks.values(ispk(spks),trigt));
            if mean(sgns) > 0
                cspk = ispk(spks(find(sgns > 0)));
            else
                cspk = ispk(spks(find(sgns < 0)));
            end
            DATA.MeanSpike.Cluster(eid,p,:) = mean(Spks.values(cspk,:));
            DATA.MeanSpike.ClusterSD(eid,p,:) = std(Spks.values(cspk,:));
            if isfield(DATA.Expts{eid}.Cluster{1,p},'autocut')
            DATA.MeanSpike.autocut(p,eid) = DATA.Expts{eid}.Cluster{1,p}.autocut;
            else
            DATA.MeanSpike.autocut(p,eid) = 0;
            end
        else
            DATA = AutoCut(DATA, eid, 1,'noplot');
%codes have changed in DATA, so reset Spks
            Spks = GetSpikeData(DATA, p);
            spks = find(Spks.codes(ispk,2) == 1);
            sgns = sign(Spks.values(ispk(spks),trigt));
            if mean(sgns) > 0
                cspk = ispk(spks(find(sgns > 0)));
            else
                cspk = ispk(spks(find(sgns < 0)));
            end
            DATA.MeanSpike.Cluster(eid,p,:) = mean(Spks.values(cspk,:));
            DATA.MeanSpike.ClusterSD(eid,p,:) = std(Spks.values(cspk,:));
            if length(ispk) < 10
            DATA.MeanSpike.dprimes(p,eid) = NaN;
            else
            DATA.MeanSpike.dprimes(p,eid) = DATA.Expts{eid}.Cluster{1,p}.dprime;
            end
            DATA.MeanSpike.autocut(p,eid) = 1;
        end
        DATA.MeanSpike.nspk(p,eid) = length(cspk);
        if DATA.MeanSpike.Cluster(eid,p,trigt) > 0
            nspk = find(Spks.values(ispk,trigt) > 0 & Spks.codes(ispk,2) ~= 1);
        else
            nspk = find(Spks.values(ispk,trigt) < 0 & Spks.codes(ispk,2) ~= 1);
        end
        DATA.MeanSpike.NotCluster(eid,p,:) = mean(Spks.values(ispk(nspk),:));
    end
    id = find(DATA.TemplateScores(p,end,:) == 0);

function Spks = GetSpikeData(DATA, p)
    if isfield(DATA,'AllSpikes')
        Spks = DATA.AllSpikes{p};
    else
        Spks = DATA.AllData.Spikes;
    end
        
function ReCombineByName(DATA, name, backup)

    load(name);
    if isfield(Expt.Header,'Combineids')
        DATA.combineids = Expt.Header.Combineids;
    else
        DATA.combineids = [];
    end
    id = find(ismember(DATA.expid,DATA.combineids));
    if ~isempty(id)
        DATA.extype = 1;
        DATA.spikelist = Expt.Header.Spikelist;
        DATA.outname = name;
        SetClusterCheck(DATA);

        nstr = length(get(DATA.elst,'String'));
        if max(id) > nstr
            fprintf('%d is longer than expt list (%d)\n',max(id),nstr);
        else
            set(DATA.elst,'value',id);
        end
        [nExpt, DATA] = CombinePlot(DATA, 0,'ids',DATA.combineids);
    end
    CompareExpts(nExpt,Expt);
    if backup
        bname = strrep(name,'.mat','bak.mat');
        if strcmp(name,bname)
            fprintf('backup (%s) same as original(%s)\n');
            return;
        else
            try save(bname,'Expt'); catch fprintf('ERROR saving %s\n',bname); return; end
        end
    end
    Expt = nExpt;
    try save(name,'Expt'); catch fprintf('ERROR saving %s\n',name); end
    
function ReCombineAll(a,b, varargin)
%
% go through all Expts. If a combined expt exists, use this list, then
% reccombine

j = 1;
lfponly = 0;
listonly = 0;
oneprobe = 0;
onetype = 0; 
mkall = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'lfponly',5)
        lfponly = 1;
    elseif strncmpi(varargin{j},'listonly',5)
        listonly = 1;
    elseif strncmpi(varargin{j},'mkall',5)
        mkall = 1;
    elseif strncmpi(varargin{j},'oneprobe',5)
        oneprobe = 1;
        lfponly = -1;
    elseif strncmpi(varargin{j},'ORBW',4)
        onetype = 1;
    elseif strncmpi(varargin{j},'spkonly',5)
        lfponly = -1;
    end
    j = j+1;
end

if isstruct(a)
    DATA = a;
else
    DATA = GetDataFromFig(a);
end
if ~isfield(DATA,'Templates')
    load('StdTemplate.mat');
    DATA.Templates = Templates;
    DATA.TemplateInfo = [];
end

spikelist = WhichClusters(DATA.toplevel);
if spikelist == -1 %no clusters selected
    it = findobj(DATA.toplevel,'Tag','UseCluster1');
    set(it,'value',1);
    spikelist = WhichClusters(DATA.toplevel);
end

nprobes = length(DATA.probelist);
if length(DATA.probelist) > 1
DATA.state.includeprobename = 1;
else
DATA.state.includeprobename = 0;
end
lfpfile = strrep(DATA.datafilename,'.mat','.lfp.mat');
if listonly == 0  && lfponly >= 0 && exist(lfpfile,'file')
fprintf('Loading %s ... ',lfpfile);
tic;
load(lfpfile);
fprintf(' Took %.1f\n',toc);
end
probedone = DATA.probe;

chk = 0;

if oneprobe | lfponly > 0
    probelist = probedone;
    oneprobe = 1;
else
    probelist = DATA.probelist;
end

d = dir(fileparts(DATA.datafilename));
if nprobes > 1
str = ['.p' num2str(probedone) 'c1.'];
else
str = ['.c1.'];
end
nf = 0;
for j = 1:length(d)
    if strfind(d(j).name,str)
       nf = nf+1;
       filenames{nf} = d(j).name;
       xfile(nf) = 1;
    end
end

DATA.state.redoautocut = 1;
oldstate.autoplotnewprobe = DATA.state.autoplotnewprobe;
DATA.state.autoplotnewprobe = 0;

if nf > 0
for j = 2:length(DATA.exptypelist)
%    DATA = ListSubExpts(DATA,j);
    [dp, poutname] = fileparts(CombinedName(DATA,j,DATA.spikelist(1),'probe',probedone));
    id = strmatch(poutname,filenames);
    if length(id)
        xfile(id) = 0;
    end
end
xid = find(xfile); %files with no match
else
    xfile  = zeros(size(DATA.exptypelist));
end


if onetype == 1
    id = strmatch('orXob',DATA.exptypelist(2:end));
    if isempty(id)
        return;
    else
        onetype = id(1)+1;
    end
    nex = 2;
else
    nex = length(DATA.exptypelist)+sum(xfile);
end

for p = 1:length(probelist)
    if listonly == 0 && oneprobe == 0
        if probelist(p) ~= DATA.probe
            DATA = SetProbe(DATA, probelist(p));
        end
    end

    for j = nex:-1:2
        if onetype
            j = onetype;
        end
        if j <= length(DATA.exptypelist)
            set(DATA.clst,'value',j);
            DATA = ListSubExpts(DATA,j);
            poutname = CombinedName(DATA,j,DATA.spikelist(1),'probe',probedone);
            donename = CombinedName(DATA,j,DATA.spikelist(1),'probe',1);
            outname = CombinedName(DATA,j,DATA.spikelist(1));
            set(DATA.saveitem,'string',outname);
            if ~isempty(DATA.allcombineids{j})
                DATA.combineids = DATA.allcombineids{j};
                fprintf('%s',DATA.exptypelist{j})
                fprintf(',%d',DATA.combineids);
                fprintf('\n');
            else
                DATA.combineids = [];
            end
        else
            poutname = [dp '/' filenames{xid(j-length(DATA.exptypelist))}];
            outname = regexprep(poutname,'\.p[0-9]*c[0-9]\.',sprintf('.p%dc1.',probelist(p)));
            donename = regexprep(poutname,'\.p[0-9]*c[0-9]\.',sprintf('.p1c1.',probelist(p)));
            set(DATA.saveitem,'string',outname);
            DATA.combineids = [];
        end
        if isempty(DATA.combineids)
            if exist(poutname,'file') | exist(donename,'file')
                if exist(poutname,'file')
                load(poutname);
                else
                load(donename);
                end
                Expt.Header.Name = BuildName(Expt.Header.Name);
                DATA.Expt = Expt;
                if j > length(DATA.exptypelist)
                    id =strmatch(Expt.Header.expname,DATA.explist,'exact');
                    if length(id) ==1
                        set(DATA.clst,'value',id);
                        DATA = ListSubExpts(DATA,id);
                    end
                end
                if isfield(Expt.Header,'Combineids')
                    DATA.combineids = Expt.Header.Combineids;
                    DATA.allcombineids{j} = DATA.combineids;
                    fprintf('%s',splitpath(poutname))
                    fprintf(',%d',DATA.combineids);
                    fprintf('\n');
                elseif isfield(Expt.Header,'Combined')
                    DATA.combineids = DATA.expid(Expt.Header.Combined);
                    DATA.allcombineids{j} = DATA.combineids;
                    fprintf('%s',splitpath(poutname))
                    fprintf(',%d',DATA.combineids);
                    fprintf('\n');
                else
                    DATA.combineids = [];
                    fprintf('%s no combineids',splitpath(poutname));
                    if DATA.logfid
                        fprintf(DATA.logfid,'%s no combineids',splitpath(poutname))
                    end
                end
            else
                DATA.combineids = [];
            end
        end
        if length(DATA.combineids)
            DATA.extype = j;
            id = find(ismember(DATA.expid,DATA.combineids));
            if isempty(id) %ids from a saved file, may not match list in GUI
                DATA.expid = DATA.combineids;
                id = find(ismember(DATA.expid,DATA.combineids));
            end

            if ~isempty(id)
                if listonly
                    fprintf('Expts ');
                    fprintf('%d ',id);
                    fprintf('of %d\n',length(DATA.expid));
                else
                    nstr = length(get(DATA.elst,'String'));
                    if max(id) > nstr
                        fprintf('%d is longer than expt list (%d)\n',max(id),nstr);
                    else
                        set(DATA.elst,'value',id);
                    end
                    if lfponly > 0
                        DATA.Expt = CombinePlot(DATA, 0);
                        if j > length(DATA.exptypelist)
                        end
                        combine('savelfp',DATA,LFP);
                    else
                        [Expt, DATA] = CombinePlot(DATA, chk);
                        %             Expt.Header.SpkStats = GetSpkStats(DATA);
                        drawnow;
                        nspk = sum([Expt.Trials.count]);
                        if j > length(DATA.exptypelist)
                            file = strrep(poutname,['.p' num2str(probedone)],['.p' num2str(probelist(p))]);
                            c = '*'
                        else
                            file = CombinedName(DATA,j,1);
                            c = '';
                        end
                        save(file,'Expt');
                        fprintf('Saved %d spikes (Expts %s) to %s (%s)\n',nspk,sprintf(' %d',DATA.combineids),file,DATA.user);
                        if DATA.logfid
                            fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s%s) to %s (%s)\n',datestr(now),nspk,c,sprintf(' %d',DATA.combineids),file,DATA.user);
                        end
                        nspk = Expt.Header.nspk;
                        nc = length(Expt.Header.nspk)-1;
                        cl = 1;
                        spikelist = DATA.spikelist;
                        while cl < nc
                            cl = cl+1;
                            if nspk(cl+1) > 10
                                DATA.spikelist = cl;
                                SetClusterCheck(DATA);
                                [Expt, DATA] = CombinePlot(DATA, chk);
                                if j > length(DATA.exptypelist)
                                    file = strrep(poutname,['.p' num2str(probedone) 'c1'],['.p' num2str(probelist(p)) 'c' num2str(cl)]);
                                    c = '*'
                                else
                                    file = CombinedName(DATA,j,cl);
                                    c = '';
                                end
                                save(file,'Expt');
                            end
                        end
                        DATA.spikelist = spikelist;
                        SetClusterCheck(DATA);
                        if p ==1  && lfponly == 0 && exist('LFP')
                            DATA.Expt = Expt;
                            combine('savelfp',DATA,LFP);
                        end
                        %                CombineAll(a,DATA);
                    end
                end
            end
        end
    end
%Calculate spike shape AFTER combining, so that autocuts are
%shown
    if onetype == 0 && listonly == 0
    DATA = CalcMeanSpike(DATA,1:length(DATA.Expts));
    SetExptClusters(DATA);
    if oneprobe ==0 && lfponly < 1
        SaveSpikeShape(DATA,DATA.meanspkfile);
    end
    end
end
if mkall
    PlotAllProbes(fileparts(DATA.datafilename),'sptrig','save');
end
DATA.state.autoplotnewprobe = oldstate.autoplotnewprobe;


set(DATA.toplevel,'UserData',DATA);

function spkstats = GetSpkStats(DATA)

% summary stats for a probe by expt, to try and track drift
spkstats = [];
ids = get(DATA.elst,'value');
exid = get(DATA.clst,'value');
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
    Spks = DATA.AllData.Spikes;
end
for j = 1:length(ids)
    Expt =  DATA.Expts{DATA.expid(ids(j))};
    mins = min(Spks.values(Expt.gui.spks,:)');
    maxs = max(Spks.values(Expt.gui.spks,:)');
    if length(mins) > 10 && length(maxs) > 10
    spkstats.max(j,:) = prctile(maxs,[50 90 99]);
    spkstats.min(j,:) = prctile(maxs,[50 10 1]);
    spkstats.h(j,:) = prctile(maxs-mins,[50 90 99]);
    spkstats.k(j) = moment(maxs-mins,3);
    else
        spkstats.max(j,:) = [0 0 0];
        spkstats.min(j,:) = [0 0 0];
        spkstats.h(j,:) = [0 0 0];
        spkstats.k(j) = 0;
    end
    if isempty(Expt.gui.spks)
            spkstats.trange(j,:) = [0 0];
    else
        spkstats.trange(j,:) = [min(Spks.times(Expt.gui.spks)) ...
            max(Spks.times(Expt.gui.spks))];
        ispk = Expt.gui.spks;
        if ~isfield(DATA,'Spikes') || max(ispk) > length(DATA.Spikes.cx) || sum(DATA.Spikes.cx(ispk)) == 0
            DATA = CalcClusterVars(DATA,ispk);
        end
            if length(ispk) > 10000
                xid = find(DATA.Spikes.cx(ispk) > prctile(DATA.Spikes.cx(ispk),99.9));
            elseif length(ispk) > 100
                v = sort(DATA.Spikes.cx(ispk),'descend');
                xid = find(DATA.Spikes.cx(ispk) > v(100));
            else
                xid = 1:length(ispk);
            end
            spkstats.V = mean(Spks.values(ispk(xid),:));
            xid = find(Spks.codes(ispk,2) == 1);
            if length(xid)
                spkstats.cV = mean(Spks.values(ispk(xid),:));
            end
    end
end



function CombineAll(a,DATA, varargin)
%
%Combine all probes, for one expt
%Don't redo LFP here. Too slow.

checkonce = 1;
guicall = 0;
j = 1;

dolfp = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',3)
        checkonce = 1;
    elseif strncmpi(varargin{j},'lfp',3)
        dolfp = 1;
    elseif strncmpi(varargin{j},'nocheck',3)
        checkonce = 0;
    end
    j = j+1;
end

if ~isstruct(DATA)
    DATA = GetDataFromFig(a);
    guicall = 1;
end
DATA.state.includeprobename = 1;
    eid = get(DATA.clst,'value');
    outname = DATA.outname;
    DATA.extype = eid;
    oldauto = DATA.plot.autoclustermode;
%hitting "combine all" online is asking to autocut uncut probes
    if DATA.plot.autoclustermode == 3 && DATA.state.online 
    DATA.plot.autoclustermode = 1;
    end
for j = 1:length(DATA.probelist)
    DATA = SetProbe(DATA, DATA.probelist(j));
    if j == checkonce
        chk = 1;
    else
        chk = 0;
    end
    Expt = CombinePlot(DATA, chk);
    if isfield(DATA,'AllClusters') && ~isfield(Expt.Header,'probesep')
    Expt.Header.probesep = 50;
    end
%    Expt.Header.SpkStats = GetSpkStats(DATA);
    drawnow;
    nspk = sum([Expt.Trials.count]);
    file = regexprep(outname,'\.p[0-9]*c1\.',['.p' num2str(DATA.probelist(j)) 'c1.']);
%    file = CombinedName(DATA,eid,1);
    if DATA.state.nospikes == 0
        save(file,'Expt');
        fprintf('Saved %d spikes to %s\n',nspk,file);
    end
    Expts{j} = Expt;
  
    if DATA.logfid
        fprintf(DATA.logfid, '%s,Saved %d spikes (Expts%s) to %s\n',datestr(now),nspk,sprintf(' %d',Expt.Header.Combineids),file);
    end
end
DATA.Expt = Expt;
DATA.plot.autoclustermode = oldauto;
if DATA.state.online == 0 && dolfp
   combine('savelfp',DATA);
end
DATA.AllExpts = Expts;
DATA.plot.autoclustermode = oldauto;
if guicall
    set(DATA.toplevel,'UserData',DATA);
end

    
function FitButton(a,b)

DATA = GetDataFromFig(a);
[DATA.Expt, plotres] = PlotCombined(DATA, DATA.Expt);
fit = FitExpt(plotres,'plot');
DATA = AddFitToData(DATA, plotres, fit);
set(DATA.toplevel,'UserData',DATA);

function DATA = AddFitToData(DATA, plotres, fit)
type = strmatch(plotres.type{1},{'Op' 'Pp'});
if ~isempty(type)
    if type(1) == 1
        DATA.fitvals.Op = fit.mean;
        DATA.fitvals.Opw = fit.sd;
        DATA.fitvals.OpRo = GetEval(DATA.Expt,'Ro');
        DATA.Expt.fits.Op = fit;
    elseif type(1) == 2
        DATA.fitvals.Pp = fit.mean;
        DATA.fitvals.PpRo = GetEval(DATA.Expt,'Ro');
        DATA.fitvals.Ppw = fit.sd;
        DATA.Expt.fits.Pp = fit;
    end
    DATA.Expt.fit = fit;
end

        
        function cntrl_box = setshow(DATA, tag)
wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
h = 220 * wsc(2);
w = 350 * wsc(1);
scrsz = get(0,'Screensize');
cw = DATA.plot.cw;
ch = DATA.plot.ch;
SpkDefs;
   bp(1) = SPACE;
   bp(2) = ch+VSPACE;
   bp(3) = cw*9;
   bp(4) = ch+VSPACE;
   dat.parentfigtag = DATA.tag.top;
cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc(2) w*wsc(1) h*wsc(2)], 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Showvals','UserData',dat);
  
fn = fields(DATA.show);
for j = 1:length(fn)
    id = strmatch(fn{j},CodeNames.Codes);
    if length(id)
        str = CodeNames.Label{id};
    else
        str = fn{j};
    end
        uicontrol(gcf,'Style', 'CheckBox','String',str,'Position', bp,...
   'Tag',fn{j},'Callback',@ShowUpdate,'value',DATA.show.(fn{j}));
bp(2) = bp(2)+ch+SPACE;
end

function str = vec2str(x)
          str = [int2str(x(1)) ':' int2str(x(end))];

function cntrl_box = setoptions(DATA, tag)

wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
ch = DATA.plot.ch;
cw = DATA.plot.cw;
h = (ch+VSPACE)*18;
w = cw * 40;
scrsz = get(0,'Screensize');

cntrl_box = findobj('Tag',tag,'Name','Options');
if ~isempty(cntrl_box)
    figure(cntrl_box);
    return;
end
if ~isfield(DATA,'figpos') | isempty(DATA.figpos{1})
   bp = get(DATA.toplevel,'Position');    
   DATA.figpos{1} =  [bp(1)+bp(3) bp(2) w*wsc(1) h*wsc(2)]
end
dat.parentfigtag = DATA.tag.top;
cntrl_box = figure('Position', DATA.figpos{1} , 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Options','UserData',dat);


top = num2str(DATA.toplevel); 
HSPACE = 0.02;
VSPACE = 0.02; 
bp(1) = HSPACE;
cols(1) = HSPACE;
cols(2) = 0.26;
cols(3) = 0.5;
cols(4) = 0.74;
rows = VSPACE + [0:14]./16;
%start from bottom
bp(3) = 0.24;
colw = bp(3);
bp(4) = mean(diff(rows));
   bp(1) = cols(1,1);
nr = 1;
nc = 1;
   bp(1) = cols(nc);
   bp(2) = rows(nr);
   a = uicontrol(gcf,'Style', 'CheckBox','String','Fit Expts','units','norm','Position', bp,...
      'Tag','AutoFit','Callback',@Update,'value',DATA.state.autofit);
    bp(1) = cols(2);
  uicontrol(gcf,'Style', 'CheckBox','String','Spool if No Cluster','units','norm','Position', bp,...
      'Tag','AutoSpool','Callback',@Update,'value',DATA.state.autospool);
    bp(1) = cols(3);
  uicontrol(gcf,'Style', 'CheckBox','String','V range','units','norm','Position', bp,...
      'Tag','AutoVrange','Callback',@Update,'value',DATA.plot.autoVrange);
if length(DATA.probelist) > 1
    bp(1) = cols(4);
  uicontrol(gcf,'Style', 'CheckBox','String','Update CellList plot','units','norm','Position', bp,...
      'Tag','AutoPlotCells','Callback',@Update,'value',DATA.state.autoplotcells);
end

bp(2) = rows(2); bp(1) = cols(1);
   uicontrol(gcf,'Style', 'CheckBox','String','Plot Graph if change','units','norm','Position', bp,...
      'Tag','AutoReplotGraph','Callback',@Update,'value',DATA.state.autoreplotgraph);
    bp(1) = cols(2);
  uicontrol(gcf,'Style', 'CheckBox','String','Set List latest','units','norm','Position', bp,...
      'Tag','AutoSetList','Callback',@Update,'value',DATA.state.autosetlist);

  bp(1) = cols(3);
     uicontrol(gcf,'Style', 'CheckBox','String','Use Last Cluster','units','norm','Position', bp,...
      'Tag','ApplyLastCluster','Callback',@Update,'value',DATA.state.applylastcluster);
    bp(1) = cols(4);
  if length(DATA.probelist) > 1
   uicontrol(gcf,'Style', 'CheckBox','String','Plot All For New Probe','units','norm','Position', bp,...
      'Tag','AutoPlotNewProbe','Callback',@Update,'value',DATA.state.autoplotnewprobe);
    end

  bp(1) = cols(1);
  bp(2) = rows(3);
   uicontrol(gcf,'Style','Text','String','Automatic actions','units','norm','Position',bp);
  
    
   
  bp(2) = rows(4);
   bp(1) = cols(1);
   uicontrol(gcf,'Style', 'CheckBox','String','Limit Range','units','norm','Position', bp,...
      'Tag','FixRange','Callback',@Update,'value',DATA.state.fixrange);
  bp(2) = rows(5);
  uicontrol(gcf,'Style', 'CheckBox','String','ShowEM','units','norm','Position', bp,...
      'Tag','ShowEM','Callback',@Update,'value',DATA.plot.showem);


  bp(2) = rows(6);
   uicontrol(gcf,'Style', 'CheckBox','String','Flip','units','norm','Position', bp,...
      'Tag','Flip','Callback',@Update,'value',DATA.plot.flip);
   bp(2) = rows(7);
   uicontrol(gcf,'Style', 'CheckBox','String','Collapse 1','units','norm','Position', bp,...
      'Tag','Collapse1','Callback',@Update,'value',DATA.plot.collapse);
  bp(2) = rows(8);
  uicontrol(gcf,'Style', 'CheckBox','String','Acov','units','norm','Position', bp,...
      'Tag','Acov','Callback',@Update,'value',DATA.plot.acov);
   bp(2) = rows(8);
   bp(1) = cols(3);
   uicontrol(gcf,'Style', 'CheckBox','String','ISIH','units','norm','Position', bp,...
      'Tag','ISIH','Callback',@Update,'value',DATA.plot.showISI);

  
  bp(2) = rows(4);
  bp(1) = cols(2);
   uicontrol(gcf,'Style', 'CheckBox','String','Force rebuild','units','norm','Position', bp,...
      'Tag','ForceBuild','Callback',@Update,'value',DATA.state.forcebuild);
  bp(2) = rows(5);
   uicontrol(gcf,'Style', 'CheckBox','String','dV vs V','units','norm','Position', bp,...
      'Tag','PhasePlot','Callback',@Update,'value',(DATA.plot.dvdt == 2));
  bp(2) = rows(6);
   uicontrol(gcf,'Style', 'CheckBox','String','Remove DC','units','norm','Position', bp,...
      'Tag','RemoveDC','Callback',@Update,'value',DATA.plot.nodc);

  bp(2) = rows(7);
   uicontrol(gcf,'Style', 'CheckBox','String','Auto relist','units','norm','Position', bp,...
      'Tag','AutoList','Callback',@Update,'value',DATA.state.autolist);

  bp(2) = rows(8);
   uicontrol(gcf,'Style', 'CheckBox','String','Condense RC','units','norm','Position', bp,...
      'Tag','Condense','Callback',@Update,'value',DATA.plot.condenseRC);

  % third column of checkboxes 
  bp(2) = rows(4);
  bp(1) = cols(3);
   uicontrol(gcf,'Style', 'CheckBox','String','.pN in name','units','norm','Position', bp,...
      'Tag','NameProbe','Callback',@Update,'value',DATA.state.includeprobename);

  bp(2) = rows(5);
   uicontrol(gcf,'Style', 'CheckBox','String','Plot F1','units','norm','Position', bp,...
      'Tag','PlotMod','Callback',@Update,'value',DATA.plot.plotmod);

  bp(2) = rows(6);
   uicontrol(gcf,'Style', 'CheckBox','String','Auto Advance','units','norm','Position', bp,...
      'Tag','AutoNext','Callback',@Update,'value',DATA.state.autonext);
  
  bp(2) = rows(7);
   uicontrol(gcf,'Style', 'CheckBox','String','Hide Spikes','units','norm','Position', bp,...
      'Tag','NoSpikes','Callback',@Update,'value',DATA.state.nospikes);
  

  % fourth column of checkboxes 
  bp(1) = cols(4);
  bp(2) = rows(4);
   if ~isfield(DATA.state, 'optimizeclusters')
       DATA.state.optimizeclusters = 0;
   end
 %
 %  uicontrol(gcf,'Style', 'CheckBox','String','Optimize','units','norm','Position', bp,...
 %     'Tag','OptimizeClusters','Callback',@Update,'value',DATA.state.optimizeclusters);
   uicontrol(gcf,'Style', 'CheckBox','String','Show Artifacts','units','norm','Position', bp,...
      'Tag','ShowArtifacts','Callback',@Update,'value',DATA.plot.showartifacts);
  bp(2) = rows(5);
   uicontrol(gcf,'Style', 'CheckBox','String','xCorr','units','norm','Position', bp,...
      'Tag','ShowxCorr','Callback',@Update,'value',DATA.plot.xcorr);
  bp(2) = rows(6);
   uicontrol(gcf,'Style', 'CheckBox','String','+hash','units','norm','Position', bp,...
      'Tag','AddHash','Callback',@Update,'value',DATA.plot.addhash);
  bp(2) = rows(7);
   uicontrol(gcf,'Style', 'CheckBox','String','Wave T','units','norm','Position', bp,...
      'Tag','ShowWave','Callback',@Update,'value',DATA.plot.showwave);
  %back to left side
  bp(2) = rows(8);
   uicontrol(gcf,'Style', 'CheckBox','String','RF relative','units','norm','Position', bp,...
      'Tag','CenterRFMeasures','Callback',@Update,'value',DATA.plot.centerRFmeasures);

   
  bp(2) = rows(9);
   bp(1) = cols(1);
   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'units','norm','Position', bp,...
      'Tag','ClusterX','Callback',@Update,'value',DATA.plot.clusterX);
   
   bp(1) = cols(2);
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterXrange(1)),'units','norm','Position', bp,...
      'Tag','ClusterXmin','Callback',@Update);
  bp(1) = cols(3); 
   uicontrol(gcf,'Style', 'pushbutton','String','Plot ISIH','units','norm','Position', bp,...
      'Callback','combine(''PlotISI'')');


   bp(1) = cols(1);
   bp(2) = rows(10);
   uicontrol(gcf,'Style', 'pop','String',DATA.spkvarnames,'units','norm','Position', bp,...
      'Tag','ClusterY','Callback',@Update,'value',DATA.plot.clusterY);

   bp(1) = cols(2);
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.clusterYrange(1)),'units','norm','Position', bp,...
      'Tag','ClusterYmin','Callback',@Update);

   bp(1) = cols(3); bp(3) = 0.1;
    uicontrol(gcf,'Style', 'text','string','AutoScale','units','norm','Position',bp);

   bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
 uicontrol(gcf,'Style', 'pop','String','Matlab|99.9|%99%|95%|Tight','units','norm','Position', bp,...
      'Tag','AutoScaleMode','Callback',@Update,'value',DATA.plot.autoscalemode);

   bp(1) = cols(4); bp(3) = 0.1;
    uicontrol(gcf,'Style', 'text','string','AutoCut','units','norm','Position',bp);

   bp(1) = cols(4)+0.1; bp(3) = colw-0.1;
 uicontrol(gcf,'Style', 'pop','String','Centiles|By Density|Test|None','units','norm','Position', bp,...
      'Tag','AutoCutMode','Callback',@Update,'value',DATA.plot.autoclustermode+1);

   bp(1) = cols(1); bp(3) = 0.1;
   bp(2) = rows(11);
   uicontrol(gcf,'Style', 'text','string','Pt Range A','units','norm','Position',bp);
   bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterArange),'units','norm','Position', bp,...
      'Tag','ClusterArange','Callback',@Update);

  bp(1) = cols(2); bp(3) = 0.1;
  uicontrol(gcf,'Style', 'text','string','B','units','norm','Position',bp);
   bp(1) = cols(2)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterBrange),'units','norm','Position', bp,...
      'Tag','ClusterBrange','Callback',@Update);
   bp(1) = cols(3); bp(3) = 0.1;
  uicontrol(gcf,'Style', 'text','string','E','units','norm','Position',bp);
   bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'edit','string',vec2str(DATA.clusterErange),'units','norm','Position', bp,...
      'Tag','ClusterErange','Callback',@Update);

   bp(1) = cols(1); bp(3) = 0.1;
   bp(2) = rows(12);
   uicontrol(gcf,'Style', 'text','string','N Min','units','norm','Position',bp);
   bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.nmin),'units','norm','Position', bp,...
      'Tag','Nmin','Callback',@Update);

   bp(1) = cols(2); bp(3) = colw;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.nminrc),'units','norm','Position', bp,...
      'Tag','RCNmin','Callback',@Update);

  if 0  %all done from main menu now
   bp(1) = cols(3); bp(3) = 0.1;
   uicontrol(gcf,'Style', 'text','string','CP','units','norm','Position',bp);
   bp(1) = cols(3)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'pop','String','None|CP-time|CP trials|CP-hist|CP-EM|Psych Only|Psych Smooth','units','norm','Position', bp,...
      'Tag','ShowCP','Callback',@Update,'value',DATA.plot.showcp+1);
  end
   bp(1) = cols(4); bp(3) = 0.1;
   uicontrol(gcf,'Style', 'text','string','PtSz','units','norm','Position',bp);
   bp(1) = cols(4)+0.1; bp(3) = colw-0.1;
  uicontrol(gcf,'Style', 'pop','String','Auto|1|2|3|4|5|6|7|8', 'units','norm','Position',bp,...
      'Tag','SetPtSize','Callback',@Update,'value',DATA.plot.setptsize+1);

  
  bp(2) = rows(13);
   bp(1) = cols(1); bp(3) = 0.1;
   uicontrol(gcf,'Style', 'text','string','sdfw','units','norm','Position',bp);
   bp(1) = cols(1)+0.1; bp(3) = colw-0.1;
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.sdfw),'units','norm','Position', bp,...
      'Tag','Sdfw','Callback',@Update);

   bp(1) = cols(2); bp(3) = 0.1;
   uicontrol(gcf,'Style', 'text','string','MLFP','units','norm','Position',bp);
   if ~isfield(DATA.plot,'lfpplot');
       DATA.plot.lfpplot = 0;
   end
   bp(1) = cols(2)+0.1; bp(3) = colw-0.1;
  uicontrol(gcf,'Style', 'pop','String','None|Default|Stack|Image|Blank|RC|Movie|OneStim|monocs|Eig|Var|BlankVar|Frameresp|Trial|CSD','units','norm', 'Position',bp,...
      'Tag','LFPPlot','Callback',@Update,'value',DATA.plot.lfpplot+1);
  

   bp(1) = cols(1);
   bp(2) = rows(14); bp(3) = colw;
   uicontrol(gcf,'Style', 'text','string','Spike Display MaxV','units','norm','Position',bp);
   bp(1) = cols(2);
  if ~isfield(DATA.plot,'SpikeMaxV')
      DATA.plot.SpikeMaxV = 5;
  end
  if ~isfield(DATA.plot,'SpikeMinV')
      DATA.plot.SpikeMinV = -5;
  end
  if ~isfield(DATA.plot,'SpikeVsep')
      DATA.plot.SpikeVsep = 3;
  end
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeMinV),'units','norm','Position', bp,...
      'Tag','SpikeMaxV','Callback',@Update);
   bp(1) = cols(3);
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeMaxV),'units','norm','Position', bp,...
      'Tag','SpikeMaxV','Callback',@Update);
   bp(1) = cols(4);
   uicontrol(gcf,'Style', 'edit','string',num2str(DATA.plot.SpikeVsep),'units','norm','Position', bp,...
      'Tag','SpikeVsep','Callback',@Update);
  
   bp(1) = cols(1);
    bp(2) = rows(15);
    bp(3) = 0.1;
   %       bp(1) = SPACE;
   if isfield(DATA,'AllSpikes')
       str = 'Sync';
   else
       str = 'Trig';
   end
   uicontrol(gcf,'Style', 'text','string',str,'units','norm','Position',bp);
   bp(3) = 0.2;
   bp(1) = cols(1)+0.1;
      uicontrol(gcf,'Style', 'pop','String','Neg|Either|Pos|None|PlotbyNeg|PlotByPos|PlotByAll', 'units','norm','Position',bp,...
      'Tag','SyncSign','Callback',@Update,'value',DATA.syncsign+2);
  if isfield(DATA,'AllSpikes')
   bp(1) = cols(2);
   uicontrol(gcf,'Style', 'CheckBox','String','Overlay','units','norm','Position', bp,...
      'Tag','SyncOverlay','Callback',@Update,'value',DATA.plot.syncoverlay)

   bp(1) = cols(3);
   uicontrol(gcf,'Style', 'CheckBox','String','SpikeBySpike','units','norm','Position', bp,...
      'Tag','SpikeBySpike','Callback',@Update,'value',DATA.plot.timebyspikeprobe)
   bp(1) = cols(4);
      uicontrol(gcf,'Style', 'pop','String','None|Min|Max|Energy|PCA1-1|PCA2-2|PCA1-2|Xcorr-PCA|CX-CX|CY-CY|Sum|test', 'units','norm','Position',bp,...
      'Tag','SyncCluster','Callback',@Update,'value',DATA.plot.synccluster+1);
  end

  
  function ShowUpdate(a,b)

DATA = GetDataFromFig(a);
fn = fields(DATA.show);
for j = 1:length(fn)
    DATA.show.(fn{j}) = GetShow(DATA,fn{j});
end
set(DATA.toplevel,'UserData',DATA);

function [value, it] = GetShow(DATA,tag)
it = findobj(DATA.showid, 'Tag',tag);
if ~isempty(it) 
    value = get(it(1),'value');
else
    value = 0;
end

function SetClusterCheck(DATA)
     SetCheck('UseCluster1',ismember(1,DATA.spikelist));
     SetCheck('UseCluster2',ismember(2,DATA.spikelist));
     SetCheck('UseCluster3',ismember(3,DATA.spikelist));
     SetCheck('UseCluster4',ismember(4,DATA.spikelist));

  

function SetGui(DATA);
   %SetCheck('Recut',DATA.state.recut > 0); %% in case it is 2
    SetCheck('Recount',DATA.state.recount);
    SetCheck('AutoPlot',DATA.state.autoplot);
    SetCheck('ShowSpikes',DATA.state.showspikes);
    SetCheck('ShowN',DATA.plot.showN);
    SetCheck('PlotPsych',DATA.state.plotpsych,DATA.toplevel);
    SetCheck('ResetClusters',DATA.state.resetclusters);
    SetClusterCheck(DATA);
    if isfigure(DATA.xyfig)
    it = findobj(DATA.xyfig, 'Tag','Density');
    if it
        if DATA.densityplot
            set(it,'value',1);
        else
            set(it,'value',0);
        end
    
    ax = findobj(DATA.xyfig,'Type','axes');
    if DATA.plot.autoscale == 1
        set(ax,'Ylimmode','auto','Xlimmode','auto');
        DATA.plot.clusterXrange  = get(ax(1),'Xlim');
        DATA.plot.clusterYrange  = get(ax(1),'Ylim');
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
        SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
        SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
    elseif ismember(DATA.plot.autoscale,[2 3])
        x = get(ax(1),'Xlim');
        y = get(ax(1),'Ylim');
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    else
        set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    end
    end
    end
    if isfigure(DATA.optionfig)
    SetField(DATA.optionfig,'ClusterXmin',DATA.plot.clusterXrange(1));
    SetField(DATA.optionfig,'ClusterYmin',DATA.plot.clusterYrange(1));
    set(findobj(DATA.optionfig,'Tag','ClusterX'),'value',DATA.plot.clusterX);
    set(findobj(DATA.optionfig,'Tag','ClusterY'),'value',DATA.plot.clusterY);
    end


function SetClusters(a,b,tag)
DATA = get(findobj('Tag',tag),'UserData');
DATA.spikelist = WhichClusters(DATA.toplevel);
if DATA.state.online == 2
    DATA.Expts = CountTxtSpikes(DATA.Expts,DATA.probe,DATA.spikelist);
end
if DATA.state.autoreplotgraph
    id = union(DATA.currentexpt,DATA.exabsid);
    for j = id
        DATA = CountSpikes(DATA,j);
    end
    combine('setexpplot',DATA);
end
set(DATA.toplevel,'UserData',DATA);

function xcorrhit(a,b, varargin);
DATA = GetDataFromFig(a);
GetFigure('CrossCorrelation');
if isfield(DATA.plot,'useprobe')
id = find(DATA.plot.useprobe)
if length(id)
    PlotXcorr(DATA, id,3); %or [DATA.probes(id).probe] if not 1:nprobes
 end
end


function SetProbeHit(a,b, varargin)

next = 0;
DATA = GetDataFromFig(a);
DATA.state.savedvdt = 0;
if DATA.savedclusters == 0
    msgbox('YOu have not saved cluster parameters for last probe yet');
    return;
end
if DATA.playingspk
    set(findobj(DATA.svfig,'Tag','StopSpool'),'value',1);
    return;
end

j = 1;
while j <= length(varargin)
    if strcmp('ReloadProbe',varargin{j})
        if isfield(DATA,'AllSpikes')
        DATA = rmfield(DATA,'AllSpikes');
        end
        if isfield(DATA,'AllClusters')
            DATA = rmfield(DATA,'AllClusters');
        end
        if isfield(DATA,'sids')
            DATA = rmfield(DATA,'sids');
        end
    elseif strcmp(varargin{j},'next')
        next = 1;
    end
    j = j+1;
end

f = get(a);
if strcmp(f.Tag,'ProbeId')
    pit = a;
else
    pit = findobj(DATA.toplevel,'Tag','ProbeId');
    sit = findobj(DATA.toplevel,'Tag','SubprobeId');
end
if isfield(f,'Tag') & strcmp(f.Tag,'SubprobeId')
    DATA.subprobe = f.Value-1;
elseif isfield(f,'Value')
    id = get(pit,'value');
    if next 
        if (DATA.subprobe == 0 || isempty(sit)) && id < length(DATA.probelist)
        id = id+1;
        set(pit,'value',id);
        elseif length(sit) && DATA.subprobe > 0
            sub = get(sit,'value');
            if sub < 5
                sub = sub+1;
                set(sit,'value',sub);
                DATA.subprobe = sub-1;
            elseif id < length(DATA.probelist)
                id = id+1;
                DATA.subprobe = 1;
                set(sit,'value',DATA.subprobe+1);
                set(pit,'value',id);
            end
        end
    end
    DATA.probe= DATA.probelist(id);
end




if isfield(DATA,'AllSpikes')
    if DATA.state.applylastcluster
    [DATA, DATA.spklist] = SetExptSpikes(DATA,DATA.currentexpt,0,'useexpt');
    else
    [DATA, DATA.spklist] = SetExptSpikes(DATA,DATA.currentexpt,'setrange');
    end
    nloaded = 0;
    for j = 1:length(DATA.AllSpikes)
        if isfield(DATA.AllSpikes{j},'codes') & length(DATA.AllSpikes{j}.codes) > 10
            nloaded = nloaded+1;
        end
    end
elseif DATA.probe == 100
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif isfield(DATA,'AllClusters')
    DATA.state.nospikes =1;
else
    DATA = SetProbe(DATA, DATA.probe);
    nloaded = 1;
end
DATA.outname = regexprep(DATA.outname,'.p[0-9]*c([0-9]).',sprintf('.p%dc$1.',DATA.probe));
set(DATA.saveitem,'string',DATA.outname);
set(DATA.toplevel,'UserData',DATA);
playspk = get(findobj(DATA.toplevel,'Tag','ShowSpikes'),'value');
DATA.state.showspikes = playspk;
%
% if > 3 probes are laoded, don't spool through them every time the probe
% hit is changed.
if DATA.state.nospikes
    DATA = combine('setexp', DATA);
    plotISI(DATA);
elseif (DATA.state.autospool | playspk) & nloaded < 4
    DATA = combine('setexp', DATA);
else
    DATA = CalcClusterVars(DATA,  DATA.spklist);
    DATA.Expts{DATA.currentexpt}.gui.spks = DATA.spklist;
    if isfigure(DATA.xyfig)
    GetFigure(DATA.xyfig);
    hold off;
    DrawXYPlot(DATA, DATA.spklist);
    end
end
set(DATA.toplevel,'UserData',DATA);


function [DATA, ok] = SetProbe(DATA, probe)

    DATA.probe = probe;
    if isfield(DATA,'toplevel')
set(DATA.toplevel,'Name',sprintf('Loading Probe %d...',DATA.probe));
    end
drawnow;
fprintf('Loading Probe %d...',DATA.probe);
tstart = now;
if ~isfield(DATA,'currentexpt')
    DATA.currentexpt = 1;
end
if isfield(DATA,'AllClusters')
    DATA.Spikes.cx = DATA.AllClusters(probe).cx;
    DATA.Spikes.cy = DATA.AllClusters(probe).cy;
%    DATA.AllData.Spikes.times = DATA.AllClusters(probe).times;
%    DATA.AllData.Spikes.codes = zeros(length(DATA.AllClusters(probe).times),4);
    for j = 1:length(DATA.Expts)
        DATA.Expts{j}.gui.classified = 0;
    end
elseif isfield(DATA,'AllSpikes')
    for j = 1:length(DATA.expid)
        ts(j,:) = DATA.Expts{DATA.expid(j)}.Header.trange;
    end
    times = [min(ts(:,1)) max(ts(:,2))];
    DATA.AllSpikes{probe} = GetProbeFiles(DATA, probe,DATA.subprobe,'trange',times/10000,'nodv');
elseif sum([DATA.probes.probe] < 100) == 1 %ustim expt with just 1 channel
    DATA.probe = probe;
else
DATA = LoadSpikes(DATA, DATA.currentexpt);
DATA = LoadClusters(DATA,ClusterFile(DATA),'noclear');
end
DATA.currentcluster = 1; %else can get blank clusters ->error
for j = 1:length(DATA.Expts)
    [DATA, DATA.Expts{j}.gui.spks] = SetExptSpikes(DATA,j,'setrange');
end
mytoc(tstart);
if isfield(DATA,'toplevel') & DATA.state.autoplotnewprobe & length(DATA.probelist) > 2 & DATA.state.online == 0 && DATA.state.nospikes == 0
    PlotAllExpts(DATA);
end
NotBusy(DATA);

function filename = GetProbeFilename(DATA, eid, probe)

id = find(DATA.probelist == probe);
        if DATA.state.online == 0
            if length(id) > 1
            else
            if DATA.probelist(id) > 16
                filename = strrep(DATA.datafilename,'.mat',sprintf('A.p%s.mat',DATA.probevars{id}(3:end)));
            else
                filename = strrep(DATA.datafilename,'.mat',sprintf('.p%s.mat',DATA.probevars{id}(3:end)));
            end
            end
        else
            filename = ['C:' DATA.Expts{eid}.Header.Name];
            if DATA.probelist(id) > 16
                filename = strrep(filename,'/Expt','A/Expt');
            end
        end

        
    function DATA = LoadSpikes(DATA, eid)

    if ~isfield(DATA,'probes') || DATA.state.nospikes
        return;
    end
        id = find([DATA.probes.probe] == DATA.probe);
% with online data, where part way through there is a change in the list of
% available probes, so use DATA.probelist
      if DATA.state.online
          id = find([DATA.probelist] == DATA.probe);
      else
 %DATA.probes.probes, unlike DATA.probelist,  will have mulitple entries if the data are split across
 % multiple files (offline only)
        id = find([DATA.probes.probe] == DATA.probe);
      end
        if DATA.state.online == 0
            if length(id) > 1 || length(DATA.probes) > 2 %KLUDGE Need to detect new files with only 1 expt
                DATA.AllData.Spikes = GetProbeFiles(DATA,DATA.probe,DATA.subprobe);
            else
                if DATA.probelist(id) > 16
                    filename = strrep(DATA.datafilename,'.mat',sprintf('A.p%s.mat',DATA.probevars{id}(3:end)));
                else
                    filename = strrep(DATA.datafilename,'.mat',sprintf('.p%s.mat',DATA.probevars{id}(3:end)));
                end
                [d, a,b]  = fileparts(filename);
                filename = [d '/Spikes/' a b];
                [DATA.AllData.Spikes]= GetProbeSpikes(DATA.AllData, filename, DATA.probevars{id},DATA.subprobe);
            end
            DATA = SetSpkLists(DATA);
            DATA.spklist = DATA.Expts{DATA.currentexpt}.gui.spks;
            DATA = SetExptSpikes(DATA, eid, 0);
        else
            filename = ['C:' DATA.Expts{eid}.Header.Name];
            if DATA.probelist(id) > 16  || DATA.probesource(id) == 2
                filename = strrep(filename,'/Expt','A/Expt');
            end
            if ~exist(filename,'file')
                [a,b] = fileparts(DATA.Expts{eid}.Header.Name);
                filename = [DATA.datafilename '/'  b '.mat'];
            end
            [DATA.AllData.Spikes]= GetProbeSpikes(DATA.AllData, filename , DATA.probevars{id},DATA.subprobe);
            DATA.spklist = [];
            DATA = SetExptSpikes(DATA, eid, 0);
            DATA = CountSpikes(DATA,eid);

            %  DATA = SetSpkCodes(DATA, [1:length(DATA.AllData.Spikes.times)], 0)
        end

function DATA = SetSpkLists(DATA)
    DATA.spklist = [];
    for j = 1:length(DATA.Expts)
        DATA.Expts{j}.gui.setispk = 0;
    end

    ifile = strrep(DATA.datafilename,'.mat',sprintf('.p%dispk.mat',DATA.probe));
    if exist(ifile)
        load(ifile);
        if length(ispklen) == length(DATA.Expts)
            for j = 1:length(DATA.Expts)
                if length(expispk{j})
                    DATA.Expts{j}.gui.spks = expispk{j};
                    DATA.Expts{j}.gui.setispk = length(expispk{j});
                    DATA.Expts{j}.gui.spkrange = [expispk{j}(1) expispk{j}(end)];
                end
            end
        end
    end
    nset = 0;
    for j = 1:length(DATA.Expts)
        if DATA.Expts{j}.gui.setispk == 0
            DATA = SetExptSpikes(DATA, j, 'setrange');
            expispk{j} = DATA.Expts{j}.gui.spks;
            ispklen(j) = length(expispk{j});
            nset = nset+1;
        end
    end
    if nset > 0
        save(ifile,'expispk','ispklen');
    end

 function Spikes = GetProbeFiles(DATA, probe, subprobe, varargin)
%GetProbeFiles loads up spikes that are spli across files.
% given a time range (in sec, not tics) only loads files needed for that
% range
    trange = [];
    dvfile = [];
    nodv = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j}, 'trange',3)
        j = j+1;
        trange = varargin{j};
    elseif strncmpi(varargin{j}, 'nodv',3)
        nodv = 1;
    end
    j = j+1;
end
    id = find([DATA.probes.probe] == probe);
    [dp,pref] = fileparts(DATA.datafilename);
    Spk.times = [];
    Spk.values = [];
    Spk.codes = [];
    [a,sid] = sort([DATA.probes(id).first]);
    sid = id(sid);
    if length(trange) > 1
%fid is files that are past the point we need
        fid = find([DATA.probes(sid).first] > trange(2));
        if isempty(fid)
            fid = length(sid)+1;
        end
        lid = find([DATA.probes(sid).last] < trange(1));
        if isempty(lid) %need first file
            sid = sid(1:fid(1)-1);
        elseif fid(1) == lid(end)+1
            sid = sid(fid(1));
        else
            sid = sid(lid(end)+1:fid(1)-1);
        end
        fprintf('Spike file %d\n',sid);
    end
    if DATA.state.savedvdt
        dvfile = [dp '/Spikes/' strrep(DATA.probes(sid(1)).filename,'t0','dvdt')];
    end
    for j = 1:length(sid)
    filename = [dp '/Spikes/' DATA.probes(sid(j)).filename];
    if ~exist(filename,'file') & isfield(DATA.probes,'pathname');
        filename = [DATA.probes(sid(j)).pathname '/' DATA.probes(sid(j)).filename];
    end
    if exist(filename,'file')
        a = load(filename);
        chname = DATA.probes(sid(j)).var;
        if isempty(Spk.times)
            Spk = a.(chname);
            ns = length(Spk.times);
        else
            Spk.times = [Spk.times; a.(chname).times];
            Spk.codes = [Spk.codes; a.(chname).codes];
            Spk.values = [Spk.values; a.(chname).values];
            if size(a.(chname).times,1) > size(a.(chname).values,1)
                fprintf('Some Missing Spike values in %s\n',filename);
                if DATA.logfid
                    fprintf(DATA.logfid,'Some Missing Spike values in %s\n',filename);
                end
            end
        end
    else
        fprintf('No file %s\n',filename);
    end    
    end

if ~isempty(Spk.values) 
    if size(Spk.values,1) < length(Spk.times)
        fprintf('Some Missing Spike values in %s\n',pref);
    end
    
    if size(Spk.values,3) > 1 & subprobe > 0
        Spk.values = Spk.values(:,:,subprobe);
    end
    if nodv == 0
    Spikes = CleanSpikes(Spk,'bufl',10000,'dvfile',dvfile); 
    else
        Spikes = Spk;
    end
    Spikes.times = Spikes.times .* 10000;
    if isfield(DATA.probes, 'firsti')
        Spikes.firstspki = DATA.probes(sid(1)).firsti;
    end
    cfile = ClusterFile(DATA,'probe',probe);
    if exist(cfile,'file')
        load(cfile);
        lastspk = DATA.probes(sid(end)).firsti+DATA.probes(sid(end)).nspk-1;
        if length(clid) < lastspk
            lastspk = length(clid);
            nspk = 1+lastspk-Spikes.firstspki;
        Spikes.codes(1:nspk,2) = clid(Spikes.firstspki:lastspk);
        elseif size(Spikes.codes,1) >  lastspk
            Spikes.codes(Spikes.firstspki:lastspk,2) = clid(Spikes.firstspki:lastspk);
        else
            Spikes.codes(:,2) = clid(Spikes.firstspki:lastspk);
        end
        Spikes.codes = Spikes.codes(:,1:2);
    end
else
    fprintf('No Spike values %s\n',filename);
    Spikes = [];
end



function Spikes = GetProbeSpikes(All, filename, varname, subprobe)
    if exist(filename,'file')
    load(filename);
    else
        fprintf('No file %s\n',filename);
        Spikes = [];
        return;
    end
    if exist(varname,'var')
%        Spikes = eval(varname);
    Spikes = eval(['CleanSpikes(' varname ', ''bufl'',10000);']);        
    Spikes.times = Spikes.times .* 10000;
    else
        fprintf('No data for %s in %s\n',varname,filename);
        Spikes = [];
    end
    if size(Spikes.values,3) > 1 & subprobe > 0
        Spikes.values = Spikes.values(:,:,subprobe);
    end

function Update(a,b)

    guichange = 0;
DATA = GetDataFromFig(a);
DATA = combine('getstate');
caller = get(a,'Tag');

if strcmp(caller,'AllDensity')
    PlotAllProbeXY(DATA);
    return;
end



it = findobj('Tag',DATA.tag.options);
if ~isempty(it)
    DATA.plot.autoclustermode = get(findobj(it,'Tag','AutoCutMode'),'value')-1;
    DATA.state.fixrange = GetCheck('FixRange');
    DATA.plot.showem = GetCheck('ShowEM',it);
    %DATA.plot.showcp = get(findobj(it,'Tag','ShowCP'),'value')-1;
    DATA.plot.condenseRC = GetCheck('Condense',it);
    DATA.plot.clusterX = get(findobj(it,'Tag','ClusterX'),'value');
    DATA.plot.clusterY = get(findobj(it,'Tag','ClusterY'),'value');
    DATA.plot.clusterXrange(1) = GetField('ClusterXmin',it);
    DATA.plot.clusterYrange(1) = GetField('ClusterYmin',it);
    DATA.clusterArange = GetField('ClusterArange',it);
    DATA.clusterBrange = GetField('ClusterBrange',it);
    DATA.clusterErange = GetField('ClusterErange',it);
    DATA.plot.nmin = GetField('Nmin',it);
    DATA.plot.nminrc = GetField('RCNmin',it);
    DATA.plot.sdfw = GetField('Sdfw',it);
    DATA.plot.SpikeMaxV = GetField('SpikeMaxV',it);
    DATA.plot.SpikeVsep = GetField('SpikeVsep',it);
    DATA.plot.DensitySigma = [3 3];
    DATA.plot.addhash = GetCheck('AddHash',it);
    DATA.plot.xcorr = GetCheck('ShowxCorr',it);
    DATA.plot.autoVrange = GetCheck('AutoVrange',it);
    DATA.plot.showartifacts = GetCheck('ShowArtifacts',it);
    if length(DATA.plot.clusterX) > 1 || length(DATA.plot.clusterY) > 1
        fprintf('ClusterX/Y is too big');
    end
    DATA.plot.acov = GetCheck('Acov',it);
    DATA.state.includeprobename = GetCheck('NameProbe',it);
    DATA.state.autofit = GetCheck('AutoFit',it);
    DATA.plot.flip = GetCheck('Flip',it);
    DATA.plot.collapse = GetCheck('Collapse1',it);
    DATA.plot.showISI = GetCheck('ISIH',it);
    [DATA.state.autolist, h] = GetCheck('AutoList');
    DATA.plot.setptsize = get(findobj(it,'Tag','SetPtSize'),'value')-1;
    DATA.plot.lfpplot = get(findobj(it,'Tag','LFPPlot'),'value')-1;
    DATA.plot.centerRFmeasures = GetCheck('CenterRFMeasures',it);
    DATA.state.autoplotnewprobe = GetCheck('AutoPlotNewProbe',it);
    DATA.state.autospool = GetCheck('AutoSpool',it);
    DATA.state.autoplotcells = GetCheck('AutoPlotCells',it);
    DATA.state.autoreplotgraph = GetCheck('AutoReplotGraph',it);
    DATA.state.autosetlit = GetCheck('AutoSetList',it);
    DATA.plot.autoscalemode = get(findobj(it,'Tag','AutoScaleMode'),'value');
    DATA.state.applylastcluster = GetCheck('ApplyLastCluster',it);
    s = get(findobj(it,'Tag','SyncSign') ,'value')-2;
    if isempty(DATA.xprobes)
        DATA.TriggerSign = s;
    end
    if isfield(DATA,'AllSpikes')
    if ~isempty(s)
        DATA.syncsign = s;
        DATA.plot.synccluster = get(findobj(it,'Tag','SyncCluster') ,'value')-1;
        DATA.plot.timebyspikeprobe = GetCheck('SpikeBySpike',it);
    end
    DATA.plot.syncoverlay = GetCheck('SyncOverlay',it);
    end
end


if (DATA.xyfig & get(a, 'Parent') == DATA.xyfig) || strcmp(caller,'AutoScaleMode')
    set(0,'currentfigure',DATA.xyfig);
    DATA.plot.autoscale = GetCheck('AutoScale',DATA.xyfig) * DATA.plot.autoscalemode;
    it = findobj(DATA.xyfig,'Tag','Clusterid');
    DATA.currentcluster = get(it,'value');
        ax = findobj(DATA.xyfig,'Type','axes');
    if DATA.plot.autoscale
        set(ax,'Ylimmode','auto','Xlimmode','auto');
        DATA.plot.clusterXrange  = get(ax,'Xlim');
        DATA.plot.clusterYrange  = get(ax,'Ylim');
        DATA = SetXYRanges(DATA);
        SetField(DATA.xyfig,'ClusterXmax',DATA.plot.clusterXrange(2));
        SetField(DATA.xyfig,'ClusterYmax',DATA.plot.clusterYrange(2));
    end
    set(ax,'Ylim', DATA.plot.clusterYrange,'Xlim',DATA.plot.clusterXrange);
    DATA.state.currentclusterquality = get(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value')-1;
    if strcmp(caller,'ClusterQuality')
        DATA.Expts{DATA.currentexpt}.Cluster{DATA.currentcluster,DATA.probe}.quality = ...
            DATA.state.currentclusterquality;
    end
    it = findobj(DATA.xyfig,'Tag','ForceClusterid');
    if it
        if strcmp(caller,'Clusterid')
            set(it,'value',DATA.currentcluster);
        end
        DATA.forceclusterid = get(it,'value');
    end
elseif DATA.xyfig & strcmp(caller,'SetPtSize');
    figure(DATA.xyfig); hold off;
    DrawXYPlot(DATA,DATA.Expts{DATA.currentexpt}.gui.spks);
end
DATA.plot.dvdt = GetCheck('dVdt');
DATA.plot.nodc = GetCheck('RemoveDC');
dvv = GetCheck('PhasePlot');
if DATA.plot.dvdt && dvv
    DATA.plot.dvdt = 2;
end

DATA.state.resetclusters = GetCheck('ResetClusters',DATA.toplevel);
DATA.state.showspkxy = GetCheck('SpkXY');
DATA.state.recount = GetCheck('Recount');
DATA.state.plotpsych = GetCheck('PlotPsych');
DATA.state.plotcombined = GetCheck('PlotCombined');
DATA.plot.plotmod = GetCheck('PlotMod');
DATA.plot.showsync = GetCheck('ShowSync');
DATA.plot.showwave = GetCheck('ShowWave');
DATA.plot.showN = GetCheck('ShowN');
DATA.state.showspikes = GetCheck('ShowSpikes');
[DATA.state.autoplot, h] = GetCheck('AutoPlot');
DATA.state.autonext = GetCheck('AutoNext');
DATA.state.nospikes = GetCheck('NoSpikes');
DATA.state.optimizeclusters = GetCheck('OptimizeClusters');
DATA.spikelist = WhichClusters(DATA.toplevel);
DATA.state.uselfp = get(findobj(DATA.toplevel,'Tag','UseLFP'),'value');
DATA.state.forcebuild = GetCheck('ForceBuild');
if isfigure(DATA.xyfig)
DATA.state.currentclusterquality = get(findobj(DATA.xyfig,'Tag','ClusterQuality'),'value')-1;
end

DATA.state.plotseq = get(findobj(DATA.toplevel,'Tag','PlotSeq'),'value')-1;
strs = get(findobj(DATA.toplevel,'Tag','PlotSeq'),'string');
DATA.plot.type = deblank(strs(DATA.state.plotseq+1,:));
if isempty(DATA.plot.autoclustermode)
    DATA.plot.autoclustermode = 1;
end
for j = DATA.probelist
    DATA.plot.useprobe(j) = GetCheck(['UseProbe' num2str(j)]);
end

if DATA.state.online %%no LFP available
    DATA.state.uselfp = 0;
end
id = regexp(DATA.outname,'\.c[0-9]\.');
if id & DATA.spikelist >= 0
    DATA.outname(id+2) = num2str(DATA.spikelist(1));
    set(DATA.saveitem,'string',DATA.outname);
end

if DATA.state.plotseq == 5
    DATA.plot.showcp = 5;
end


needplot = 0;
if strmatch(caller,{'PlotSeq'})
    guichange = 1;
    s = get(a,'value');
    if ismember(s,[1 2 3])% no cp or psych 
        DATA.state.plotpsych = 0;
        DATA.plot.showcp = 0;
    elseif s == 8 %CP - trials
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 2;
    elseif s == 9 %CP - Hist
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 3;
    elseif s == 6 %Psych Only
        DATA.state.plotpsych = 1;
        DATA.plot.showcp = 5;
    else
        DATA.plot.showcp = 0;
    end
    PlotCombined(DATA,DATA.Expt);
    needplot = 0;
elseif strmatch(caller,{'Psych'})
    PlotCombined(DATA,DATA.Expt);
    needplot = 0;
elseif strmatch(caller,{'SpkXY' 'ShowSpikes'})
    if DATA.state.showspkxy || DATA.state.showspikes
        DATA.state.nospikes = 0;
    end
elseif strmatch(caller,{'LFPPlot'}) & isfield(DATA,'LFP')
    ShowLFPPlot(DATA);
    needplot = 1;
end
set(DATA.toplevel,'UserData',DATA);
    if strmatch(caller,{'ClusterX' 'ClusterY'})
        if isfield(DATA,'spklist') && ~isempty(DATA.spklist)
        expspks = DATA.spklist;
    else
        expspks = DATA.spkrange(1):DATA.spkrange(2);
        end
    GetFigure(DATA.xyfig);
    hold off;
    DATA = CalcClusterVars(DATA, expspks);
    if DATA.plot.setptsize
        DATA.ptsize = DATA.plot.setptsize;
    end
    DrawXYPlot(DATA,expspks);
    end

% used to read
%if DATA.state.autoplot & a ~= h & ~ DATA.state.showspikes
% but h is not set. What was this?

 id = get(DATA.elst,'value');

 if DATA.state.autoreplotgraph & needplot & (~DATA.state.showspikes | length(id) > 1)
    combine('setexp','flagchange');
end

if guichange
    SetGui(DATA);
end
    
function SetField(parent, tag, value, varargin)
if isfigure(parent)
    it = findobj(parent,'Tag', tag);
else
    it = findobj('Tag', tag);
end
if it
    set(it,'string',num2str(value));
end

function value = GetField(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = str2num(get(it(1),'string'));
else
    value = NaN;
end

function [value, it] = GetCheck(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = get(it(1),'value');
else
    value = 0;
end

function [success] = SetCheck(tag, value,  varargin)

if nargin == 3 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    set(it(1),'value',value);
    success = 1;
else
    success = 0;
end
    
function args = PlotArgs(DATA, Expt,varargin)

combined = 0;
if DATA.state.plotcombined == 1
    combined = 1;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'combined',6)
        combined = 1;
    end
    j = j+1;
end

if isfield(DATA.tag,'rcfiga')
    args = {'rcfiga' DATA.tag.rcfiga 'figb' DATA.tag.rcfigb};
else
    args = {};
end
if isempty(DATA.plot.showcp)
    DATA.plot.showcp = 0;
end
if DATA.plot.showcp > 0
    psych = 1;
else
    psych = get(findobj('Tag','PlotPsych','Parent',DATA.toplevel),'value');
end
if ~isfield(Expt.Trials,'RespDir') %% no psych here
    psych = 0;
end
seq = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
on = get(findobj('Tag','ShowN','Parent',DATA.toplevel),'value');

if on
    args = {args{:} 'shown'};
end
on = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
if psych & ismember(seq,[4])
        args = {args{:} 'psych' 'cpseq'};
% noplot stops PSF from being plotted, still shows seq
%        args = {args{:} 'psychnoplot' 'cpseq'};
else
    if psych
        args = {args{:} 'psych'};
    end
    if seq == 1
        args = {args{:} 'seqt'};
    elseif seq == 3
        args = {args{:} 'seqid'};
    elseif seq == 2
        args = {args{:} 'sequence'};
    elseif seq == 12 && combined == 1
        args = {args{:} 'sdfsub' 'result.x==0'};
    elseif seq == 13 && combined == 1
        args = {args{:} 'sdf' 'collapse' 1};
    elseif seq == 14 && combined == 1
        args = {args{:} 'sdf' 'collapse' 2};
    elseif seq == 15 && combined == 1
        args = {args{:} 'plotall'};
    elseif seq == 9 && combined == 1
        args = {args{:} 'sdfall'};
    elseif seq == 7 && combined == 1
        args = {args{:} 'cpseq'};
    elseif seq == 8 && combined == 1
        args = {args{:} 'cphist'};
    end
end

if ismember(DATA.plot.showcp,[1 2 3 4])
    args = {args{:},  'CPtag', DATA.tag.cptag};
end
if DATA.plot.showcp == 1
    args = {args{:} 'psych'  'cpt'};
elseif DATA.plot.showcp == 2 && combined == 0
    args = {args{:} 'cp'};
elseif DATA.plot.showcp == 3 && combined == 0
    args = {args{:} 'cphist'};
elseif DATA.plot.showcp == 4 && combined == 0
    args = {args{:} 'cpt' 'emdiff' DATA.em 'emskip' DATA.plot.emskip};
end

needsdf = 0;
for j = 1:length(args)
    if ischar(args{j}) && strncmp(args{j},'sdf',3)
        needsdf = 1;
    end
end
if needsdf
    args = {args{:} 'preperiod' DATA.state.preperiod 'postperiod' DATA.state.postperiod};
end
if DATA.plot.showem
    args = {args{:} 'eyem'};
end
if DATA.plot.addhash
    args = {args{:} 'showmu'};
end

if (DATA.plot.condenseRC || combined == 0) && isempty(strmatch(DATA.plot.type ,'Subspace'))
    args = {args{:} 'condense'};
end
if DATA.plot.centerRFmeasures
    args = {args{:} 'centerRF'};
end

if DATA.state.uselfp
    args = {args{:} 'lfpt'};
end

if DATA.plot.nminrc > 0
    args = {args{:} 'rcnmin' DATA.plot.nminrc};
end
if DATA.plot.nmin > 0
    args = {args{:} 'nmin' DATA.plot.nmin};
end
if DATA.plot.sdfw > 0 
    args = {args{:} 'sdfw' DATA.plot.sdfw};
end

if ~DATA.plot.condenseRC && Expt.Stimvals.st == 4 && strcmp(Expt.Stimvals.et,'Ol')
        args = {args{:} 'twoslice'};
end

if strfind(Expt.Header.expname,'tfXip')
    args = {args{:} 'sxcx'};
end

if DATA.plot.collapse
    args = {args{:} 'collapse' 2};
end
if DATA.plot.flip
    args = {args{:} 'reverse'};
end
if DATA.plot.acov
    args = {args{:} 'acov'};
end
if DATA.plot.plotmod
    args = {args{:} 'mod'};
end

    
function mousept = myellipse(mousept, finish)


%
%
% mousept modes
%  1 draw ellipse
%  2 re-size both dimensions
%  3 rotate
%  4 move ellipse
%  5 change left edge (move center and h radius)
if nargin > 1
    start = mousept.start;
else
    mousept.mode = 4;
end
if mousept.mode == 1
    a = (finish(1,1)-start(1,1))/2;
    b = (finish(2,2)-start(2,2))/2;
    mousept.r = abs([a b]);
    mousept.c = [finish(1,1)+start(1,1) finish(2,2)+start(2,2)]/2;
elseif mousept.mode == 2
    a = (finish(1,1)-mousept.c(1));
    b = (finish(2,2)-mousept.c(2));
    mousept.r = abs([a b]);
elseif mousept.mode == 3
    
    dy = (finish(2,2)-mousept.c(2))./mousept.yrange;
    dx = (finish(1,1)-mousept.c(1))./mousept.xrange;
    t = -dy/dx;
    mousept.angle = atan(-dy/dx);
    a = mousept.r(1);
    b = mousept.r(2);
%   [a,b] = exy(mousept.angle,mousept.r(1),mousept.r(2));
%    fprintf('%.3f %.3f %.3f\n',mousept.angle,a,b);
%    a = (mousept.r(1) * cos(mousept.angle)) + mousept.ratio * mousept.r(1) * sin(mousept.angle)
%    b = (mousept.r(2) * cos(mousept.angle)) - (mousept.r(2) * sin(mousept.angle)/mousept.ratio)
    %b = mousept.r(2);
elseif mousept.mode == 4  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
elseif mousept.mode == 11  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
    mousept.c(1) = finish(1,1)-mousept.offset(1);
    mousept.c(2) = finish(2,2)-mousept.offset(2);
elseif mousept.mode == 6 %% move R boundary
    a = (finish(1,1)-mousept.c(1));
    mousept.c(1) = mousept.c(1) + (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 7 %% move L boundary
    a = abs(mousept.c(1) - finish(1,1));
    mousept.c(1) = mousept.c(1) - (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 8 %% move Top boundary
    b = (finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) + (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 9 %% move Bottom boundary
    b = abs(finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) - (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);

else %just draw what we have
    a = mousept.r(1);
    b = mousept.r(2);
end

a = a./mousept.xrange;
b = b./mousept.yrange;
sn = sin(mousept.angle);
cn = cos(mousept.angle);
%sn = sin(pi/4);
%cn = cos(pi/4);
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = mousept.xrange .* (x .* cn + y .*sn) + mousept.c(1);
yr = mousept.yrange .* (y .* cn - x .*sn) + mousept.c(2);
mousept.lasth = plot(real(xr),real(yr),'color',mousept.color);


function ScrollWheel(src, evnt)
global mousept;

DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

mousept.angle = mousept.angle+0.02*evnt.VerticalScrollCount;
mousept.mode = 10;
if mousept.lasth & ishandle(mousept.lasth)
    delete(mousept.lasth);
end
mousept = myellipse(mousept,[0 0; 0 0 ]);


function ScrollTrial(src, evnt)

DATA = GetDataFromFig(src);

if src ~= gcf
    return;
end

DATA = GetDataFromFig(src);
PlayOneTrial(DATA,src,sign(evnt.VerticalScrollCount));

   
function KeyPressed(a, ks)
global mousept;

parent = get(a,'parent');
mousept.mode;
if strmatch(ks.Key,'delete') & mousept.mode == 5
    DeleteCluster(mousept.cluster, a);
    mousept.mode = 0;
    mousept.angle = 0;
    if mousept.lasth & ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
elseif ks.Key == 'n'
    NewCluster(a);
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
elseif ks.Key == 'r'
    mousept.angle = mousept.angle+0.02;
    mousept.mode = 10;
    if mousept.lasth & ishandle(mousept.lasth)
    delete(mousept.lasth);
    end
   mousept = myellipse(mousept,[0 0; 0 0 ]);
end

function KeyReleased(a, ks)
global mousept;

mousept.mode;
if strmatch(ks.Key,'delete') & mousept.mode == 5
elseif ks.Key == 'r'
  ClassifySpikes(mousept, a);
end

function NewCluster(a)
DATA = GetDataFromFig(a);
j = size(DATA.cluster,1);
while isempty(DATA.cluster{j,DATA.probe}) & j > 1
    j = j-1;
end
newc = j+1;
it = findobj('Tag','Clusterid');
set(it,'value',newc);

function DeleteCluster(cl,callfig, varargin)
DATA = GetDataFromFig(callfig);
redraw = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nodraw',4)
        redraw = 0;
    end
    j = j+1;
end
p = get(gca,'UserData');
if isempty(p)
p = DATA.probe;
end
if cl == 0
    for j = 1:size(DATA.cluster,1)
        DeleteCluster(j, callfig, varargin{:});
    end
    return;
end
if ~isempty(DATA.cluster{cl,p}) & isfield(DATA.cluster{cl,p},'h') &  ishandle(DATA.cluster{cl,p}.h)
    delete(DATA.cluster{cl,p}.h);
end
DATA.cluster{cl,p} = [];
DATA.cluster{cl,p}.touched = 1; %records active deletion. Will be copied to Expt
DATA.cluster{cl,p}.quality = 0; %records active deletion. Will be copied to Expt
DATA.cluster{cl,p}.deleted = 1; %records active deletion. Will be copied to Expt
if isfield(DATA.Expts{DATA.currentexpt}.Cluster{cl,p},'forceid')
    DATA.Expts{DATA.currentexpt}.Cluster{cl,p} = rmfield(DATA.Expts{DATA.currentexpt}.Cluster{cl,p},'forceid');
end
DATA.forceclusterid = 0;
DATA.Expts{DATA.currentexpt}.Cluster{cl,p}.quality = 0;
if isfield(DATA.Expts{DATA.currentexpt}.gui,'spkrange')
    ispk = DATA.Expts{DATA.currentexpt}.gui.spkrange;
    ispk = [ispk(1):ispk(2)];
    if isfield(DATA,'AllSpikes')
    spks = find(DATA.AllSpikes{DATA.probe}.codes(ispk,2) == cl);
    DATA.AllSpikes{DATA.probe}.codes(ispk(spks),2) = 0;
    else
    spks = find(DATA.AllData.Spikes.codes(ispk,2) == cl);
    DATA.AllData.Spikes.codes(ispk(spks),2) = 0;
    end
    if redraw
    PlotSpikeXY(DATA,ispk(spks),DATA.spkcolor{1});
    end
% Shouldn't need this. not deleting all clusters. And if we were would
% imply that we are setting not clusters
% DATA.Expts{DATA.currentexpt}.gui.classified = 0;
end
set(DATA.toplevel,'UserData',DATA);
if cl > 1
    newc = cl-1;
elseif iscluster(DATA.cluster,cl+1,DATA.probe) == 1
    newc = cl+1;
else
    newc = 1;
end
it = findobj('Tag','Clusterid');
set(it,'value',newc);

function FigButtonReleased(src, data)
global mousept;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 0;
mousept.lasth = 0;
DATA = GetDataFromFig(src);
ExcludeTrials(DATA, mousept);

function ExcludeTrials(DATA, mousept)
    if DATA.state.plotseq == 2
        id = find([DATA.Expt.Trials.Trial] > mousept.start(1,1) & ...
            [DATA.Expt.Trials.Trial] < mousept.finish(1,1));
        ex = [DATA.Expt.Trials(id).Trial];
        if isfield(DATA.Expt,'ExcludeCluster')
            DATA.Expt.ExcludeCluster{1} = union(DATA.Expt.ExcludeCluster{1}, ex);
        else
            DATA.Expt.ExcludeCluster{1} = ex;
        end
    end

        
    hold off;    
    DATA.Expt = PlotCombined(DATA, DATA.Expt);
    set(DATA.toplevel,'UserData',DATA);

function mousept = myrect(mousept, finish)
    start = mousept.start;
    mousept.finish = finish;
    mousept.siz = finish-start;
    if mousept.mode == 1 % horizontal line
        mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) start(2,2)],'r');
    else
        mousept.lasth = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)],'r');
    end

function FigButtonPressed(src, data)
global mousept;
set(src, 'WindowButtonMotionFcn',@FigButtonDragged);

mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
mousept.down = 1;
mousept.lasth = 0;
mousept.drags = 0;
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');
end
mousept.l = mousept.start;
mousept.siz = [0 0];

function FigButtonDragged(src, data)
global mousept;

if mousept.down
    pt = get(gca,'CurrentPoint');
    if mousept.lasth & ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
    mousept= myrect(mousept,pt);
%   plot(pt(1,1),pt(2,2),'+');
    mousept.drags = mousept.drags+1;
end

function ButtonPressed(src, data)
global mousept;

DATA = GetDataFromFig(src);


it = findobj(src,'Tag','Clusterid');
if ~isempty(it)
    mousept.cluster = get(it,'value');
else
    mousept.cluster = 1;
end

start = get(gca,'CurrentPoint');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
if start(1,1) < xl(1) || start(1,1) > xl(2) || start(1,2) < yl(1) || start(1,2) > yl(2)
    return; %out of range
end
mousept.drags = 0;

mousept.color = DATA.spkcolor{mousept.cluster+1};
mousept.lasth = 0;
lastkey = get(gcf,'CurrentCharacter');

% don't want the axis rescaling as we draw the ellipse
set(gca,'Xlimmode','Manual','Ylimmode','Manual');
mousept.down = 1;
oldmode = mousept.mode;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
if oldmode == 10  && mousept.mode ==1% a left button press after scrolling
    mousept.mode = 10;
    return;
end
p = get(gca,'UserData');
mousept.xrange = diff(get(gca,'Xlim'));
mousept.yrange = diff(get(gca,'Ylim'));

if isempty(p)
p = DATA.probe;
end
if mousept.cluster <= size(DATA.cluster,1) && p <= size(DATA.cluster,2)
    C = DATA.cluster{mousept.cluster,p};
else
    C.h = 0;
end
mousept;
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');

    if size(DATA.cluster,1) >= mousept.cluster & size(DATA.cluster,2) >= p &...
            isfield(DATA.cluster{mousept.cluster,p},'h')
        for j = 1:size(DATA.cluster,1)
            distance(j) = DistanceToCluster(DATA.cluster{j,p},mousept.start(1,1:2));
        end
%        a = DistanceToEllipse(mousept,mousept.start(1,1:2))
        if(min(distance) > 1.05) % pressed outside = start over
            if ishandle(C.h)
            delete(C.h);
            end
        else  %% pressed inside; select this cluster, move if mouse moves
            [d, cl]= min(distance);
            C = DATA.cluster{cl,p};
            if ~isempty(it)
                set(it,'value',cl);
            end
            fi = findobj(src,'Tag','Clusterid');
            if isfield(C,'forceid') && C.forceid > 0
                DATA.forceclusterid = C.forceid;
                if ~isempty(fi)
                    set(fi,'value',C.forceid);
                end
            else
                DATA.forceclusterid = 0;
                if ~isempty(fi)
                    set(fi,'value',cl);
                end
            end
            DATA.cluster{cl,p}.h = DrawCluster(DATA.cluster{cl,p}, DATA.spkcolor{cl+1});
            mousept.mode = 5;
            mousept.cluster = cl;
            mousept.color = DATA.spkcolor{mousept.cluster+1};
            mousept.c = [DATA.cluster{cl,p}.x(1) DATA.cluster{cl,p}.y(1)];
            mousept.r = [DATA.cluster{cl,p}.x(2) DATA.cluster{cl,p}.y(2)];
            DATA.cluster{cl,p}.touched = 1;
            mousept.offset = mousept.start(1,1:2) - mousept.c;
            mousept.angle = -DATA.cluster{cl,p}.angle;
            mousept.lasth = DATA.cluster{cl,p}.h;
            set(DATA.cluster{cl,p}.h,'linewidth',2);
            if oldmode == 5 %second press in ellipse - move it
                mousept.down = 1;
                mousept.mode = 11;
                mousept.lasth = DATA.cluster{cl,p}.h;
            else
                mousept.down = 1;  %%ignore drag, release
            end
            set(DATA.toplevel,'UserData',DATA);
        end
    else
        p
    end
elseif mousept.mode == 2 %R button
     if isfield(C,'h') & C.h & ishandle(C.h) 
         delete(C.h); 
     end
     cl = mousept.cluster;
     mousept.c = [DATA.cluster{cl,p}.x(1) DATA.cluster{cl,p}.y(1)];
     mousept.r = [DATA.cluster{cl,p}.x(2) DATA.cluster{cl,p}.y(2)];
     mousept.angle = -DATA.cluster{cl,p}.angle;
     mousept
    mousept.start = get(gca,'CurrentPoint');
     if mousept.start(1,1) > mousept.c(1) + mousept.r(1)
         mousept.mode = 6;
     elseif mousept.start(1,1) < mousept.c(1) - mousept.r(1)
         mousept.mode = 7;
     elseif mousept.start(2,2) > mousept.c(2) + mousept.r(2)
         mousept.mode = 8;
     elseif mousept.start(2,2) < mousept.c(2) - mousept.r(2)
         mousept.mode = 9;
     end
     fprintf('Start %.2f,%.2f C %.2f,%.2f, R%.2f %.2f, mode %d\n',...
         mousept.start(1,1),mousept.start(2,2),mousept.c(1),mousept.c(2),mousept.r(1),mousept.r(2),mousept.mode);
elseif mousept.mode == 3 % R button
     if ishandle(C.h) delete(C.h); end

end


function distance = DistanceToCluster(C, pos);
   
if isempty(C) | ~isfield(C,'x');
    distance = NaN;
    return;
end
xy = pos - [C.x(1) C.y(1)];
xy = xy ./ [C.x(3) C.y(3)];
cn = cos(-C.angle);
sn = sin(-C.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = p./[C.x(2)./C.x(3) C.y(2)./C.y(3)];
distance = sum(distance.^2);

function distance = DistanceToEllipse(E, pos);
   
if isempty(E) | ~isfield(E,'c');
    distance = NaN;
    return;
end
xy = pos - E.c;
xy = xy ./ E.r;
cn = cos(-E.angle);
sn = sin(-E.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = p./[E.r(1)./E.xrange E.r(2)./E.yrange];
distance = sum(distance.^2);


    
function PlotSpikeXY(DATA, spkid, color)
ptsize = CheckPtSize(DATA, length(DATA.spklist)); %don't set diff sizes for diff clusters!
plot(DATA.Spikes.cx(spkid),DATA.Spikes.cy(spkid),'.',...
    'color',color,'markersize',ptsize);

function PlotDprime(DATA)

  hold off;
 Cx = DATA.Spikes.cx;
 Cy = DATA.Spikes.cy;
 Spks = DATA.AllData.Spikes;
 cspks = DATA.Expts{DATA.currentexpt}.gui.spks;
 C = DATA.cluster{DATA.currentcluster,DATA.probe};
 x = (Cx(cspks) - C.x(1))./C.x(3);
 y = (Cy(cspks) - C.y(1))./C.y(3);
 xr = x .* cos(C.angle) + y .* sin(C.angle);
 yr = y .* cos(C.angle) - x .* sin(C.angle);
 id = find(Spks.codes(cspks,2) == DATA.currentcluster);
 nid = find(Spks.codes(cspks,2) ~= DATA.currentcluster);
 xc = mean(xr(id));
 yc = mean(yr(id));
 sy = ((yr-mean(yr(id)))./std(yr));
 sx = ((xr-mean(xr(id)))./std(xr));
 o =C.dprimepar(1);
 d = sx .* cos(o) + sy.* sin(o);
 sd = std(d);
 [y,x] = smhist(d,'sd',sd/10);
 h = plot(sx,sy,'.','markersize',DATA.ptsize);
 axis('equal');
 refline(tan(o));
 hold on;
 yscale = diff(get(gca,'ylim'))./max(y);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx,ry);
 [y,x] = smhist(d(id),'sd',sd/10);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx,ry,'r');
 [y,x] = smhist(d(nid),'sd',sd/10);
 rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
 ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
 plot(rx,ry,'g');
 
        
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
   expspks = DATA.AllClusters(p).spklist;
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
    fprintf('%s: %d Clusters for Probe %d\n',DATA.explabels{DATA.currentexpt},length(nc),p);
end
end
DATA.state.recut = 1;
SetGui(DATA);
if isfield(DATA,'AllClusters')
id = find(DATA.AllClusters(p).codes(expspks) == 0);
elseif isfield(DATA,'AllSpikes')
id = find(DATA.AllSpikes{p}.codes(expspks,2) == 0);
else
id = find(DATA.AllData.Spikes.codes(expspks,2) == 0);
end
if ~DATA.densityplot
% done by SetSpkCodes 
%    PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{1});
end

%set Expts.cluster before countspikes, because that can copy Expts.cluster
%to DATA.cluster
p = get(gca,'UserData');
if isempty(p)
DATA.Expts{DATA.currentexpt}.Cluster = DATA.cluster;
end

if DATA.state.autoreplotgraph
    DATA = CountSpikes(DATA, DATA.currentexpt,'replot');
else
    DATA = CountSpikes(DATA, DATA.currentexpt);
end
DATA.Expts{DATA.currentexpt}.gui.classified = 2;
if DATA.plot.showdprime
    GetFigure('DprimeCalc');
    PlotDprime(DATA);
end
plotISI(DATA);
set(0,'currentfigure',src);
set(DATA.toplevel,'UserData',DATA);
if DATA.state.verbose
DATA.cluster{DATA.currentcluster,DATA.probe}
end

function plotISI(DATA)
        
    if DATA.plot.showISI
        GetFigure('ISI');
        isis = CalcISI(DATA.Expts{DATA.currentexpt}.Trials);
        id = find(isis < 1000)
        hist(isis(id),100);
    end

function [dprime, details] = CalcDprime(DATA,cluster, varargin)
    
    trackcl = 0;
    yr = DATA.Spikes.cy(DATA.spklist);
    xr = DATA.Spikes.cx(DATA.spklist);
    ispk = DATA.spklist;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'track',4)
            trackcl = 1;
        elseif strncmpi(varargin{j},'smooth',4)
            j = j+1;
            smw = varargin{j};
        elseif strncmpi(varargin{j},'xparam',4)
            j = j+1;
            Spikes = DATA.AllData.Spikes;
            xr = GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),varargin{j}, 1,[]);
        elseif strncmpi(varargin{j},'yparam',4)
            j = j+1;
            Spikes = DATA.AllData.Spikes;
            yr = GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),varargin{j}, 1,[]);
        end
        j = j+1;
    end
    spkr = minmax(DATA.spklist);

    id = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == cluster);
    nid = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == 0);
    sy = ((yr-mean(yr(id)))./std(yr));
    sx = ((xr-mean(xr(id)))./std(xr));
    sd = abs(sx+i*sy);
    o = [0:pi/40:pi];
    for j = 1:length(o)
        d = sx .* cos(o(j)) + sy.* sin(o(j));
        ds{j} = d;
%
%        dprimes(j) = abs(mean(d(nid(tid)))-mean(d(id)))./sqrt(mean([var(d(nid(tid))) var(d(id))]));
         dprimes(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
%        dprimes(j) = abs(mean(d(nid))-mean(d(id)))./std(d(nid));
    end
    [dprime, maxi] = max(dprimes);
    detials.distances = ds{maxi};
    details.id = id;
    details.angle = o(maxi);
    if trackcl
        smw = 0;
        d = ds{maxi};
        T = DATA.Expts{1}.Trials;
        for j = 1:length(T);
            ispk = [];
            codes = [];
            for k = max([1 j-smw]):min([length(T) j+smw])
                [a, ts, b] = FindSpikes(DATA,[T(k).Start(1) T(k).End(end)+500], DATA.probe, spkr);
                ispk = [ispk; a];
                codes = [codes; b];
            end
            nid = ispk(find(codes == 0));
            did = ispk(find(codes == cluster));
            dp(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
        end
        details.dps = dp;
    end

    
function dp = CalcClusterdp(DATA, cl)  %old method using ROC of radial distances.
    expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
    if length(expspks) < 10
        dp = 0;
        return;
    end
    if isfield(DATA,'AllSpikes')
        id = find(DATA.AllSpikes{DATA.probe}.codes(expspks,2) == cl-1);
        nid = find(DATA.AllSpikes{DATA.probe}.codes(expspks,2) ~= cl-1);
    else
    id = find(DATA.AllData.Spikes.codes(expspks,2) == cl-1);
    nid = find(DATA.AllData.Spikes.codes(expspks,2) ~= cl-1);
    end
    rx = std(DATA.Spikes.cx(expspks(id)));
    cx = mean(DATA.Spikes.cx(expspks(id)));
    ry = std(DATA.Spikes.cy(expspks(id)));
    cy = mean(DATA.Spikes.cy(expspks(id)));
    dx = DATA.Spikes.cx(expspks(id))- mean(DATA.Spikes.cx(expspks(id)));
    dy = DATA.Spikes.cy(expspks(id))- mean(DATA.Spikes.cy(expspks(id)));
    dc = abs(dx+i*dy);
    dx = DATA.Spikes.cx(expspks)- mean(DATA.Spikes.cx(expspks(id)));
    dy = DATA.Spikes.cy(expspks)- mean(DATA.Spikes.cy(expspks(id)));
    d = abs(dx+i*dy);
    [s, did] = sort(d);
    dc = d;
    dc(id) = 0;
    d(nid) = 0;
% need to order these somehow first before doing fpos/drate
% ? rank order d before setting id, nid to zero.

%detection rate
    drate = cumsum(dc(did)) ./ sum(dc);
%false positive rate
fpos = cumsum(d(did)) ./ sum(d);
%area gives ROC
dp = trapz(fpos,drate);


function DATA = DrawXYPlot(DATA, expspks)
    
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
    DATA.ptsize = CheckPtSize(DATA,length(expspks));
    delete(get(gca,'children'));
    expid = DATA.currentexpt;
    if isfield(DATA,'AllSpikes')
        expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        if ~isfield(DATA,'nclusters')
            DATA.nclusters = 0;
        end
        if length(expspks)
        for j = 1:DATA.nclusters+1
            id = find(DATA.AllSpikes{DATA.probe}.codes(expspks,ctype) == j-1);
            if length(id)
                PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
                hold on;
            end
        end
        end
    else
        if length(expspks) < 10 
            expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        end
        if DATA.densityplot
            PlotXYDensity(DATA.Spikes.cx(expspks),DATA.Spikes.cy(expspks));
        else
            nc = max(DATA.AllData.Spikes.codes(expspks,ctype));
            for j = 1:nc+1
                id = find(DATA.AllData.Spikes.codes(expspks,ctype) == j-1);
                PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{j});
                hold on;
            end
            if DATA.plot.showartifacts
                id = find(DATA.AllData.Spikes.codes(expspks,ctype) == 8);
                PlotSpikeXY(DATA, expspks(id), DATA.spkcolor{9});
            end
        end
    end
    hold on;
    DATA = DrawClusters(DATA, DATA.cluster, 0);
    if isfield(DATA.Expts{expid},'OnlineCluster') & DATA.state.showonlineclusters
        DATA = DrawClusters(DATA,DATA.Expts{expid}.OnlineCluster,0);
    end

        CalcClusterdp(DATA,1);
        if showhist
            xl = get(gca,'xlim');
            [a,b] = smhist(DATA.Spikes.cx(expspks));
            plot(b,(a-xl(1)).* diff(xl)./max(a));
        end
    if ismember(DATA.plot.autoscale,[2 3]) && length(DATA.Expts{DATA.currentexpt}.gui.spks)
        expspks = DATA.Expts{DATA.currentexpt}.gui.spks;
        DATA.spklist = expspks;
        [xr, yr] = ClusterRange(DATA.cluster,DATA.probe);
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
    DATA = FinishXYPlot(DATA);
    if ~ho
        hold off;
    end
    
function PlotDDF(DATA)
    cspks = DATA.Expts{DATA.currentexpt}.gui.spks;
    plottype = 0;
    nclusters = size(DATA.cluster,1);
    p = DATA.probe;
    for cl = nclusters:-1:1
        id = find(DATA.AllData.Spikes.codes(cspks,2) == cl);
        nid = find(DATA.AllData.Spikes.codes(cspks,2) ~= cl);
        sx = std(DATA.Spikes.cx(cspks(id)));
        sy = std(DATA.Spikes.cy(cspks(id)));
        mx = mean(DATA.Spikes.cx(cspks(id)));
        my = mean(DATA.Spikes.cy(cspks(id)));
        if p <= size(DATA.cluster,2)
        C = DATA.cluster{cl,p};
        x = (DATA.Spikes.cx(cspks) - mx)./sx;
        y = (DATA.Spikes.cy(cspks) - my)./sy;
        xr = x .* cos(C.angle) + y .* sin(C.angle);
        yr = y .* cos(C.angle) - x .* sin(C.angle);
        d = sqrt(((yr.^2 + xr.^2)));
        dprime = (mean(d(id)) - mean(d(nid)))./sqrt(mean([var(d(id)) var(d(nid))]));
        x = (DATA.Spikes.cx(cspks) - C.x(1))./C.x(2);
        y = (DATA.Spikes.cy(cspks) - C.y(1))./C.y(2);
        ddf = sqrt((y).^2 + (x).^2);
        GetFigure('DDF');
        if plottype == 1
        hold off;
        plot(xr(id),yr(id),'r.');
        hold on;
       plot(xr(nid),yr(nid),'.');
       axis('image');
        elseif plottype == 2
       hist(ddf,500);
        else
            hold off;
        [y,x] = smhist(ddf,'sd',0.1,'xval',[0:0.1:10]);
        plot(x,y);
        hold on;
        [y,x] = smhist(ddf(id),'sd',0.1,'xval',[0:0.1:10]);
        plot(x,y,'r');
        ddfprime = (mean(ddf(id)) - mean(ddf(nid)))./sqrt(mean([var(ddf(id))  var(ddf(nid))]));

        x = (DATA.Spikes.cx(cspks) - C.x(1))./C.x(3);
        y = (DATA.Spikes.cy(cspks) - C.y(1))./C.y(3);
        xr = x .* cos(C.angle) + y .* sin(C.angle);
        yr = y .* cos(C.angle) - x .* sin(C.angle);
        ddf = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
        drprime = (mean(ddf(id)) - mean(ddf(nid)))./sqrt(mean([var(ddf(id))  var(ddf(nid))]));
        
        [y,x] = smhist(d,'sd',0.1,'xval',[0:0.1:10]);
        plot(x,y,':');
        hold on;
        [y,x] = smhist(d(id),'sd',0.1,'xval',[0:0.1:10]);
        plot(x,y,'r:');
        end
        title(sprintf('D = %.4f (%.4f,%.4f)',-dprime,-drprime,-ddfprime))
        id = cspks(find(d < 1));
        ddf = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
        end
    end
        
function DATA = FinishXYPlot(DATA)
    if ismember(DATA.plot.autoscale,[0 2 3 4]) % 2 means I autoscale
        DATA = SetXYRanges(DATA);
        set(gca,'Xlim',DATA.plot.clusterXrange,'Ylim',DATA.plot.clusterYrange);

    end


        
        
    xlabel(DATA.spkvarnames{DATA.plot.clusterX});
    ylabel(DATA.spkvarnames{DATA.plot.clusterY});
    it = findobj(DATA.xyfig, 'Tag','Clusterid');
    if isempty(it)
        c = 1;
    else
        c = get(it,'value');
    end
    ei = DATA.currentexpt;
    if ~isempty(DATA.explabels{DATA.currentexpt})
        expname = DATA.explabels{DATA.currentexpt};
        expname = DATA.Expts{DATA.currentexpt}.Header.expname;
     else
        expname = DATA.Expts{DATA.currentexpt}.Header.expname;
    end
    if DATA.firsttrial > 1
        expname = [expname sprintf('from %d',DATA.Expts{DATA.currentexpt}.Trials(DATA.firsttrial).Trial)];
    end
    p = DATA.probe;
    str = '';
 
    if isfield(DATA,'cluster') & iscluster(DATA.cluster,c,p)
        if isfield(DATA.cluster{c,p},'dprime')
            str = sprintf('d=%.2f',DATA.cluster{c,p}.dprime);
        end
        if isfield(DATA,'AllSpikes')
            id = find(DATA.AllSpikes{DATA.probe}.codes(DATA.spklist,2) == c);
        else
        id = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == c);
        end
        if isfield(DATA.Expts{ei}.Stimvals,'du')
            rate = length(id)./(length(DATA.Expts{ei}.Trials) .* DATA.Expts{ei}.Stimvals.du);
        else
            rate = length(id)./length(DATA.Expts{ei}.Trials);
        end
       
    title(sprintf('%s: C%d %sX%s %s P%d %.1fHz',expname,c,DATA.spkvarnames{DATA.cluster{c,p}.params(1)},DATA.spkvarnames{DATA.cluster{c,p}.params(2)},str,DATA.probe,rate)); 
    else
    title(sprintf('%s: C%d',expname,c)); 
    end
    if DATA.plot.clusterZ > 0
        Plot3DClusters(DATA,0);
    end
function res = isExptCluster(E,c,p)
    if ~isfield(E,'Cluster');
        res = 0;
    else
        res = iscluster(E.Cluster,c,p);
    end
    

    function res = iscluster(C,c,p)
    if isempty(C) | size(C,1) < c | size(C,2) < p | isempty(C{c,p})
        res = 0;
    elseif ~isfield(C{c,p},'params') | ~isfield(C{c,p},'x')
        res = 0;
    elseif ~isfield(C{c,p},'nspk') || C{c,p}.nspk == 0
        res = 2;
    else
        res = 1;
    end


function ptsize = CheckPtSize(DATA, nspk)
    if DATA.plot.setptsize
        ptsize = DATA.plot.setptsize;
    elseif ~isfield(DATA,'ptsize')
        ptsize = 1;
    elseif nspk < 2000
        ptsize = 10;
    elseif nspk < 5000
        ptsize = 6;
    else
        ptsize = 1;
    end

 function [DATA, dprime, details] = SetSpkCodes(DATA, expspks, probe, show, varargin)

     onecluster = 0;
j = 1;
while j <= length(varargin)
    j = j+1;
end
if DATA.plot.synccluster & DATA.syncsign ~= 2
    onecluster = DATA.currentcluster;
end
dprime = 0;
details.nc = 0;
if isempty(expspks)
    return;
end
if isfield(DATA,'AllClusters')
  DATA.AllClusters(probe).codes(expspks) = 0;
  Cx = DATA.AllClusters(probe).cx;
  Cy = DATA.AllClusters(probe).cy;
elseif isfield(DATA,'AllSpikes')
    if onecluster
        id = find(DATA.AllSpikes{probe}.codes(expspks,2) == onecluster);
        DATA.AllSpikes{probe}.codes(expspks(id),2) = 0;
    else
        DATA.AllSpikes{probe}.codes(expspks,2) = 0;
    end
  if iscluster(DATA.cluster,1,probe) & sum(ismember(DATA.cluster{1,probe}.params, [DATA.plot.clusterX DATA.plot.clusterY]) < 2)
      if ~isfield(DATA.AllSpikes{probe},'dVdt')
      DATA.AllSpikes{probe}.dVdt = diff(DATA.AllSpikes{probe}.values,1,2);
      end
      Spikes = DATA.AllSpikes{probe};
      if isfield(DATA.AllSpikes{probe},'pcs')
      PCs = DATA.AllSpikes{probe}.pcs;
      else
      PCs = DATA.AllSpikes{probe}.codes;
      end
  end
  if isfield(DATA.AllSpikes{probe},'cx')
%used to check 
%     && DATA.plot.synccluster == 0
%but now Allspikes has cx. AND must use this - DATA.Spikes.cx  indices only
%match one probe
      Cx = DATA.AllSpikes{probe}.cx;
      Cy = DATA.AllSpikes{probe}.cy;
  elseif isfield(DATA,'Spikes')
      Cx = DATA.Spikes.cx;
      Cy = DATA.Spikes.cy;
  else
      return;
  end
elseif isfield(DATA,'Spikes')
  Cx = DATA.Spikes.cx;
  Cy = DATA.Spikes.cy;
    DATA.AllData.Spikes.codes(expspks,2) = 0;
    Spikes = DATA.AllData.Spikes;
    if isempty(DATA.AllData.pcs) || length(DATA.AllData.pcs) < max(expspks)
        PCs = DATA.AllData.Spikes.codes;
    else
        PCs = DATA.AllData.pcs;
    end
else
    return; % No spikes
end

%
%really want to limit this to spike in the current scope;id = find(DATA.Spks.cluster == 0);


if ~isfield(DATA,'cluster')
 nclusters = 0;
 DATA.cluster = {};
end
nclusters = size(DATA.cluster,1);
p = probe;
if nclusters >1
%    fprintf('%s: %d Clusters for Probe %d\n',DATA.explabels{DATA.currentexpt},nclusters,p);
end
sumplot = 0;
details.maxx = max(Cx(expspks));
details.maxy = max(Cy(expspks));
if onecluster
    cllist = onecluster;
elseif nclusters > 7 %includes and artifact. Do this last
    cllist = [nclusters:-1:1 8];
else
    cllist = nclusters:-1:1;
end
for cl = cllist
    cspks = expspks;
    nspk = 0;
    if p > size(DATA.cluster,2)
        C = [];
    else
        C = DATA.cluster{cl,p};
    end
    if isfield(C,'Cluster') && ~isempty(C.Cluster) 
        splitlist = 1;
    else
        splitlist = 0;
    end
    if isfield(DATA.cluster{cl,p},'forceid') & DATA.cluster{cl,p}.forceid
        clid = DATA.cluster{cl,p}.forceid;
    else
        clid = cl;
    end
    while ~isempty(C) & isfield(C,'x')
        if ~isfield(C,'lastspk')  & ~isempty(DATA.spklist)
            C.lastspk = DATA.spklist(end);
        end
        if ~isfield(C,'firstspk') & ~isempty(DATA.spklist)
            C.firstspk = DATA.spklist(1);
        elseif isfield(C,'firstspk') & C.firstspk > 0 & splitlist
            cspks = expspks(find(expspks >= C.firstspk & expspks <= C.lastspk));
        end
    if C.params(1) == DATA.plot.clusterX  %%otherwise need to calc values again)     
        x = (Cx(cspks) - C.x(1))./C.x(3);
    else
        x = GetSpikeVals(DATA, cspks, Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(1), 1,PCs(cspks,:));
        x = (x- C.x(1))./C.x(3);
    end
    if C.params(2) == DATA.plot.clusterY  %%otherwise need to calc values again)     
        y = (Cy(cspks) - C.y(1))./C.y(3);
    else
        y = GetSpikeVals(DATA, cspks, Spikes.values(cspks,:), Spikes.dVdt(cspks,:),C.params(2), 1,PCs(cspks,:));
        y = (y- C.y(1))./C.y(3);
    end
    xr = x .* cos(C.angle) + y .* sin(C.angle);
    yr = y .* cos(C.angle) - x .* sin(C.angle);
    d = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
    id = find(d <= 1);
    nid = find(d > 1);
    nspk = nspk + length(id);
    if length(id) & length(nid)
    dprime(cl) = (mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
    else
        dprime(cl) = NaN;
    end
    tic;
    if length(nid) > 10000
        pcrit = 1;
    elseif length(nid) > 1000
        pcrit = 10;
    elseif length(nid) > 200
        pcrit = 20;
    else
        pcrit = 50;
    end
    %    sy = ((yr-mean(yr(id)))./std(yr(id)));
%    sx = ((xr-mean(xr(id)))./std(xr(id)));   
    sy = ((yr-mean(yr(id)))./std(yr));
    sx = ((xr-mean(xr(id)))./std(xr));
    sd = abs(sx+i*sy);
    o = [0:pi/40:pi];
    for j = 1:length(o)
        d = sx .* cos(o(j)) + sy.* sin(o(j));
        a = prctile(abs(d(nid)),pcrit);
        tid = find(abs(d(nid))< a);
        ds{j} = d;
%
%        dprimes(j) = abs(mean(d(nid(tid)))-mean(d(id)))./sqrt(mean([var(d(nid(tid))) var(d(id))]));
         dprimes(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
         bii(j) = (1+skewness(d).^2)./(kurtosis(d)+3);
%        dprimes(j) = abs(mean(d(nid))-mean(d(id)))./std(d(nid));
    end
    [dprime(cl), maxi] = max(dprimes);
    DATA.cluster{cl,p}.dprimepar = [o(maxi)];
    details.truedprime = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
    details.bii = max(bii);
    bii = max(bii);
    debug = 0;
    if debug
        [a,b] = smhist(ds{maxi}(id),'sd',0.1);
        [c,d] = smhist(ds{maxi},'sd',0.1);
        [e,f] = smhist(ds{maxi}(tid),'sd',0.1);
        of = gcf;
        GetFigure('DDF');
        hold off;
        plot(b,a);
        hold on;
        plot(d,c,'r');
        plot(f,e,'g');
        title(sprintf('Spks %d/%d/%d',length(ds{maxi}),length(id),length(tid)));
        figure(of);
    end
    details.nc(cl) = length(id);

    
    if isfield(DATA,'AllClusters')
        DATA.AllClusters(probe).codes(cspks(id)) = clid;
    elseif isfield(DATA,'AllSpikes')
        DATA.AllSpikes{probe}.codes(cspks(id),2) = clid;
    else
        DATA.AllData.Spikes.codes(cspks(id),2) = clid;
    end
    DATA.cluster{cl,p}.dprime = dprime(cl);
    DATA.cluster{cl,p}.bii = max(bii);
    if isfield(C,'Cluster')
        C = C.Cluster;
    else
        C = {};
    end
    if show &&  (~DATA.densityplot && ~DATA.alldensityplot)
        DATA.ptsize = CheckPtSize(DATA, length(expspks));
        plot(Cx(cspks(id)),Cy(cspks(id)),'.',...
            'color',DATA.spkcolor{clid+1},'markersize',DATA.ptsize);
        sumplot = sumplot + length(id);
        if sumplot
            hold on;
        end
        ei = DATA.currentexpt;
        if isfield(DATA.Expts{ei}.Stimvals,'du')
            rate = length(id)./(length(DATA.Expts{ei}.Trials) .* DATA.Expts{ei}.Stimvals.du);
        else
            rate = length(id)./length(DATA.Expts{ei}.Trials);
        end
        if isfield(DATA,'explabels') & cl == DATA.currentcluster
    title(sprintf('%s dprime %.2f(%.2f) %.1fHz:%d/%dspks',DATA.explabels{DATA.currentexpt},dprime(cl),bii,rate,length(id),length(cspks))); 
        end
    else
        if isfield(DATA,'explabels') && length(DATA.explabels) >= DATA.currentexpt & cl == DATA.currentcluster
    title(sprintf('%s dprime %.2f %d/%dspks',DATA.explabels{DATA.currentexpt},dprime(cl),length(id),length(cspks))); 
        end
    end
    end
    DATA.cluster{cl,p}.nspk = nspk;
end %for cl = cllist
if length(dprime) >= DATA.currentcluster
    dprime = dprime(DATA.currentcluster);
else
    dprime = dprime(1);
end
    if ~DATA.densityplot && show > 1 && ~DATA.alldensityplot %show unclassified also
    if isfield(DATA,'AllClusters')
        id = DATA.AllClusters(probe).codes(cspks) == 0;
    elseif isfield(DATA,'AllSpikes')
        id = DATA.AllSpikes{probe}.codes(cspks,2) == 0;
    else
        id = DATA.AllData.Spikes.codes(cspks,2) == 0;
    end
        id = cspks((find(id)));
        plot(Cx(id),Cy(id),'.','color',DATA.spkcolor{1},'markersize',DATA.ptsize);
        hold on;

    end


function ButtonReleased(src, data)
global mousept;

%mousept
if mousept.down == 0 %do nothing
    if mousept.mode == 11
        mousept.mode = 0;
    end
   % delete(mousept.lasth);
   % mousept.lasth = 0;
% Even if touching a cluster in AllProbeSpikes may need to
% classify the spikes. Will not have happened during spooling.
   if ~strcmp(get(src,'tag'),'AllProbeSpikes')
    return;
   end
end
mousept.down = 0;
pt = get(gca,'CurrentPoint');
    if mousept.lasth & ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
  mousept= myellipse(mousept,pt)
  if mousept.mode == 5 %presed inside ellise without dragging
      set(mousept.lasth,'linewidth',2);
  else
      ClassifySpikes(mousept, src);
      mousept.mode = 0;
  end
%mousept.drags

function ButtonDragged(src, data)

global mousept;

%fprintf('mode %d %.4f,%d\n',mousept.mode,mousept.angle,mousept.drags);
if mousept.down && mousept.mode == 5
    mousept.drags = mousept.drags+1;
    if mousept.drags > 2
        mousept.down = 1;
        mousept.mode = 11;
    end
elseif mousept.down
    pt = get(gca,'CurrentPoint');
    if mousept.lasth && ishandle(mousept.lasth)
        delete(mousept.lasth);
    end
    mousept= myellipse(mousept,pt);
%    plot(pt(1,1),pt(2,2),'+');
    mousept.drags = mousept.drags+1;
end


function Expts = CheckSaccades(Expts, emfile)

load(emfile);
saccth = 10;

for j = 1:length(Expts)
    crate = Expt.Header.CRrates(1) * 10000;
    for k = 1:length(Expts{j}.Trials)
        if strfind(Expts{j}.Trials(k).OptionCode,'+2a')
            bs = Expts{j}.Trials(k).TrialStart;
            es = Expts{j}.Trials(k).End(end);
            [tdiff, emt] = min(abs([Expt.Trials.Start] - bs));
            if tdiff < 160
                maxi = min([length(Expt.Trials(emt).lh) length(Expt.Trials(emt).rh) length(Expt.Trials(emt).lv) length(Expt.Trials(emt).rv)]);
                ch = (Expt.Trials(emt).lh(1:maxi)+Expt.Trials(emt).rh(1:maxi))/2;
                cv = (Expt.Trials(emt).lv(1:maxi)+Expt.Trials(emt).rv(1:maxi))/2;
                dv = smooth(diff(cv),5);
                dh = smooth(diff(ch),5);
                v = 10000 * sqrt(dv.^2+dh.^2)./crate;
                first = (bs -Expt.Trials(emt).ftime)/crate;
                last =  (es -Expt.Trials(emt).ftime)/crate;
                times = ftime + [1:length(cv)] .*crate;
                saccs = find(v > 10);
                id = find(saccs > last);
                if ~isempty(id)
                    [rate, pk] = max(v(saccs(id)));
                    pk = id(pk);
                    sv(1) = mean(cv(first:saccs(pk)-10));
                    sv(2) = mean(cv(saccs(pk)+10:end));
                    sh(1) = mean(ch(first:saccs(pk)-10));
                    sh(2) = mean(ch(saccs(pk)+10:end));
                    saccsize = abs(diff(sh) + i* diff(sv));
                    sacdir = atan2(diff(sv),diff(sh));
                end
            end
        end
    end
end


function name = BuildName(name)

if isempty([strfind(name,'/') strfind(name,'\')])
    name = [pwd '/' name];
end
id = strfind(name,':');
if id
    if isunix
        name = ['/bgc' name(id(1)+1:end)];
        name = strrep(name,'\','/');
    else
        name = name(id(1)+1:end);
    end
else
    name = name;
end


function [aid, bid, n] = FindSync(at,bt,dt,varargin)
%from two lists of times, find events synchronous witthin dt
n = 0;
k = 1;
aid = [];
bid = [];
if isempty(at) || isempty(bt)
    return;
end
if length(bt) < 101
    bts = [1:length(bt)];
else
    bts = [1:101];
end
nb = length(bt);
for j = 1:length(at)
    t = bt(bts)-at(j);
    while(t(end) < 0 & bts(end) < nb)
        if bts(end)+100 > nb
            bts = bts(end):nb;
        else
            bts = bts+100;
        end
        t = bt(bts)-at(j);
    end
    [a,b] = min(abs(t));
    if a <= dt
        n = n+1;
        aid(n) = j;
        bid(n) = bts(b);
    end
end

function [E,V, pc] = SpikePCA(spk, probes, ids)

if length(probes) == 2
    spks = [spk{probes(1)}.values(ids{1},:) spk{probes(2)}.values(ids{2},:)];
         [E,V] = eig(cov(spks));
        pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
        pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
        pc(3,:) = (E(1:32,end)'-mean(E(1:32,end-1))) * spks(:,1:32)';
        pc(4,:) = (E(33:end,end)'-mean(E(33:end,end-1))) * spks(:,33:end)';
        pc(5,:) = (E(1:32,end-1)'-mean(E(1:32,end-1))) * spks(:,1:32)';
        pc(6,:) = (E(33:end,end-1)'-mean(E(33:end,end-1))) * spks(:,33:end)';
end

function DATA = CheckForPCA(DATA, ispk, force)
    needpca = 0;
    p = DATA.probe;
    if isempty(DATA.cluster)
        needpca = 0;
    elseif sum(ismember([DATA.plot.clusterX DATA.plot.clusterY DATA.plot.clusterZ],[33 34 37 38 39])) %need PCA
       needpca = 1;
   elseif p <= size(DATA.cluster,2) & isfield(DATA.cluster{1,p},'params') & sum(ismember(DATA.cluster{1,p}.params,[33 34 37 38]))
       needpca = 1;
   end
   if needpca
       if isfield(DATA,'AllSpikes')
           for j = p
               if ~isempty(DATA.AllSpikes{j})
                   [a,b,c] = OneProbePCA(DATA,ispk);
                   DATA.AllSpikes{j}.pcs(ispk,:) = c';
                   DATA.AllSpikes{j}.EigVal = diag(b);
               end
           end
       elseif force || length(DATA.AllData.pcs) < max(ispk) || std(DATA.AllData.pcs(ispk,1)) < 0.01
           id = find(DATA.AllData.Spikes.codes(ispk,2) < 1);
        [a,b,c] = OneProbePCA(DATA,ispk,id);
        DATA.AllData.pcs(ispk,:) = c';
        DATA.AllData.EigVal = diag(b);
        DATA.AllData.EigVec = a;
       end
   end


function [E,V, pc] = OneProbePCA(DATA, ids, subid)

    if isfield(DATA,'AllSpikes')
        spks = DATA.AllSpikes{DATA.probe}.values(ids,:);
    else
        spks = DATA.AllData.Spikes.values(ids,:);
    end
    if nargin  == 3
         [E,V] = eig(cov(spks(subid,:)));
    else
         [E,V] = eig(cov(spks));
    end
        pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
        pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
        pc(3,:) = (E(:,end-2)'-mean(E(:,end-2))) * spks';
        pc(4,:) = (E(:,end-3)'-mean(E(:,end-3))) * spks';


        
function aid = PlotArtifacts(DATA)
   
aid = [];    
        Spks = DATA.AllData.Spikes;
 id = find(Spks.codes(:,2) == 8)
amean = mean(Spks.values(id,:));
scores = Spks.values * (amean-mean(amean))';
GetFigure('Artifacts');
hold off;
x = max(Spks.values(:,10:15),[],2);
plot(x,scores,'.','markersize',1);
a = questdlg('Classify?','popup','Cancel','OK','OK');
if strcmp(a,'OK')
E = AddEllipse(gcf,'wait','color','r','line','timeout',30);
if ~isempty(E)
C = XYClassify(x,scores,E);
hold off; 
plot(Spks.values(C.id,:)','color','b');
hold on;
plot(amean,'r','linewidth',2);
plot(mean(Spks.values(C.id,:)),'g','linewidth',2);

a = questdlg('Apply?','popup','Cancel','OK','OK');
if strcmp(a,'OK')
    aid = C.id;
end
end
end

