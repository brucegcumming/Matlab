function varargout = PlotClusters(name, varargin)
%PlotClusters(dir, ......)
%reads in Cluster files made by AllVPcs and plots up clusters for each
%expt/cell
%PlotClusters(dir, 'load') preloads all of the cluter files. Without this,
%files are loaded as needed.
%PlotClusters(dir, 'loadfullv') also forces loading of all the FullV files.
%This allows it to call AllVPcs quickly for any expt.  Not implemented for
%Utah arrays yet
%        BE SURE YOU HAVE ENOUGH MEMORY 

if isunix
    prefsdir = '/bgc/group/matlab/preferences/PlotClusters';
else
    if isdir('Z:/group')
        prefsdir = 'Z:/group/matlab/preferences/PlotClusters';
    else
        prefsdir = '/bgc/group/matlab/preferences/PlotClusters';
    end
end
defaultconfig = [prefsdir '/' gethostname '.' GetUserName '.config'];
defaultlayout = [prefsdir '/' gethostname '.' GetUserName '.layout.mat'];
TAGTOP = 'PlotClusters';
DATA.loadautoclusters = 0;
loadargs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
elseif strncmpi(varargin{j},'loadauto',6)
    DATA.loadautoclusters = 1;
elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        TOPTAG = varargin{j};
    end
    j = j+1;
end
initcall =0 ;
warnmode = 'nowarn';
expts = [];
runcmd = [];
if nargout >0
    varargout{1} = [];
end

if ishandle(name)  % callback
    gui = name;
    DATA = GetDataFromFig(name);
    name = varargin{2};
    if length(varargin) > 2 && isempty(varargin{1})
    varargin = varargin{3:end};
    else
        varargin = {};
    end
else
    it = findobj('Tag',TAGTOP,'type','figure');
    if isempty(it)
        initcall = 1;
        DATA.tag.top = TAGTOP;
        DATA = SetDefaults(DATA);
    else
        DATA = get(it(1),'UserData');
        Clusters = getappdata(DATA.toplevel,'Clusters');
        gui = 0;
    end
end

if isstruct(name) || iscell(name)
    DATA.name = 'ClusterArg';
    DATA.strings{1} = 'ClusterArg';
    if isstruct(name) && isfield(name,'cls') %template fits
        Clusters = ReadTemplateClusters(name.cls);
        DATA.Templates = name.Templates;
        DATA.datatype = 2;
    elseif isstruct(name) %one cluster result;
        Clusters{1} = name.Clusters;
    elseif isfield(name{1},'mahal')
        Clusters{1} = name;
    elseif initcall == 1
        Clusters = name;
    end
    if initcall
    DATA = InitInterface(DATA); 
    DATA.defaultlayoutfile = defaultlayout;
    DATA.defaultconfigfile = defaultconfig;
    ApplyLayout(defaultlayout);
    ApplyConfig(defaultconfig);
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
    end
    
elseif isdir(name)
    d = dir(name);
    strings = {};
    for j  = 1:length(d)
       if strfind(d(j).name,'ClusterTimes.mat') & isempty(strfind(d(j).name,'OldClusterTimes.mat')) & ...
           (isempty(strfind(d(j).name,'AutoClusterTimes.mat')) | DATA.loadautoclusters)
            strings = {strings{:} d(j).name};
        end
    end
    DATA.name = regexprep(name,'/$','');
    DATA.datadir = DATA.name;
    DATA.strings = strings;
    DATA.dates = zeros(size(DATA.strings));
    DATA.tag.top = TAGTOP;
    DATA = SetDefaults(DATA);
    DATA = InitInterface(DATA);
elseif strncmpi(name,'close',4) %must be in front of exist test because of close.m
    f = fields(DATA.tag);
    for j = 1:length(f)
        CloseTag(DATA.tag.(f{j}));
    end
    return;
elseif strncmpi(name,'profile',4) %must be in front of exist test because of close.m
    DATA.profiling = ~DATA.profiling;
    fprintf('Profiling %d\n',DATA.profiling);
    
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif sum(strcmp(name,{'checkclusters' 'test'})) %catch arguments that exist says are files
elseif exist(name,'file')
    DATA.tag.top = TAGTOP;
    DATA = SetDefaults(DATA);
    DATA.exptname = name;
    [name, exptname] = fileparts(name);
    d = dir(name);
    strings = {};
    for j  = 1:length(d)
        if strfind(d(j).name,'ClusterTimes.mat') & isempty(strfind(d(j).name,'OldClusterTimes.mat'))
            strings = {strings{:} d(j).name};
        end
    end
    DATA.name = name;
    if isdir(name)
        DATA.datadir = name;
    else
        DATA.datadir = fileparts(name);
    end
    DATA.strings = strings;
    DATA.dates = zeros(size(DATA.strings));

    DATA = InitInterface(DATA);
end

if strcmp(name,'getstate')
    varargout{1} = DATA;
    return;
end
doload = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'LoadAuto',8)
        DATA = LoadAll(DATA,'auto','force');
        set(DATA.lstui,'string',DATA.strings);
    elseif strncmpi(varargin{j},'check',4)
        runcmd = varargin{j}
    elseif strncmpi(varargin{j},'updatearray',10)
        runcmd = varargin{j}
    elseif strncmpi(varargin{j},'config',4)
        j = j+1;
        if regexp(varargin{j},'^[A-Z]:') | strfind(varargin{j},'/') %real path
            DATA.configfile = varargin{j};
        else
            DATA.configfile = [prefsdir '/' varargin{j} '.config'];
        end
        if exist(DATA.configfile)
            DATA = ApplyConfig(DATA);
        end
    elseif strncmpi(varargin{j},'layout',4)
        j = j+1;
        if regexp(varargin{j},'^[A-Z]:') | strfind(varargin{j},'/') %real path
            DATA.layoutfile = varargin{j};
        else
            DATA.layoutfile = [prefsdir '/' varargin{j} '.layout.mat'];
        end
        if exist(DATA.layoutfile)
        load(DATA.layoutfile);
        ApplyLayout(DATA);
        end
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'exptload',7)
        DATA = LoadExpts(DATA);
    elseif strncmpi(varargin{j},'Load',4)
        doload = varargin{j};
    elseif strncmpi(varargin{j},'mahaltype',7)
    j = j+1;
    DATA.mahaltype = varargin{j};
    elseif strncmpi(varargin{j},'rebuild',6)
        DATA.rebuild = 1;
    elseif strncmpi(varargin{j},'norebuild',6)
        DATA.rebuild = 0;
    elseif strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
%        DATA.templatesrc = 1; 
    elseif strncmpi(varargin{j},'templateid',10)
        j = j+1;
        DATA.templateid = varargin{j};
    elseif strncmpi(varargin{j},'useauto',7)
        DATA.useautoclusters = 1;
    elseif strncmpi(varargin{j},'usealltrials',7)
        loadargs = {loadargs{:} varargin{j}};
        DATA.usealltrials = 1;
    end
    j = j+1;
end

if ~isempty(doload)
        if ~isfield(DATA,'name') && ischar(name)
            fprintf('%s is not a directory\n',name);
            return;
        elseif isdir(DATA.name)
          if strncmpi(doload,'LoadSpikes',8)
            DATA = LoadAll(DATA,[],'loadspikes',warnmode);
          else
            DATA = LoadAll(DATA,[],'force', warnmode,'expts',expts,loadargs{:});
          end
            set(DATA.lstui,'string',DATA.strings);
            if strncmpi(doload,'Loadfullv',9)
                DATA = LoadFullVs(DATA);
            end
        else
            fprintf('%s is not a directory\n',DATA.name);
            return;
        end
        
end

if strncmpi(name,'close',4)
    f = fields(DATA.tag);
    for j = 1:length(f)
        CloseTag(DATA.tag.(f{j}));
    end
    return;
elseif strncmpi(name,'checkexpts',7)
    if length(varargin)
        CheckExpts(DATA,varargin{:});
    else
        CheckExpts(DATA,'errs');
    end
elseif strncmpi(name,'checktimes',8)
        PlotTimeRanges(DATA);
elseif strncmpi(name,'checkspikes',7)
    res = CheckSpikeFiles(DATA,'spkid');
    b = CheckSpikeFiles(DATA,'build');
    if isfield(res,'acts') && isfield(b,'acts')
        res.acts = cat(2,res.acts,b.acts);
    end
    varargout{1} = res;

elseif strncmpi(name,'checkclusters',7)
    GuiMenu(DATA,[],'checkclusters');
elseif strncmpi(name,'checklist',7)
    CheckExclusion(DATA);
elseif strncmpi(name,'convertlist',7)
    DATA = ConvertExclusion(DATA);
    if strcmp(name,'convertlistsave')
        SaveCellList(DATA);
    end
elseif strncmpi(name,'correlogram',7)
    PlotCorrelogram(Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)});
elseif strncmpi(name,'findcells',7)
    CellFinderTest(DATA,1);
elseif strncmpi(name,'quickspikes',7)
%    QuickSpikes(DATA, DATA.currentpoint);
    SpoolSpikes(DATA, DATA.currentpoint);
elseif strncmpi(name,'reloadtrials',7)
    DATA = LoadTrialLists(DATA);
elseif strncmpi(name,'loadextra',7)
    DATA = LoadExtra(DATA,0);
elseif strncmpi(name,'reloadcellist',9)
    DATA = LoadCellFile(DATA);
elseif strncmpi(name,'readres',7)
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
elseif strncmpi(name,'CallAllVPcs',7)


elseif strncmpi(name,'loadspikes',7)
    LoadAllSpikes(DATA);

elseif strncmpi(name,'refitxysdx',7)
    DATA = ReFitAll(DATA, 'sdindex');
elseif strncmpi(name,'refit1d',7)
    ReFit1D(DATA, 1.1);
elseif strncmpi(name,'refit3means',7)
    ratio = 1.1;
    if length(varargin) & isnumeric(varargin{1})
        ratio = varargin{1};
    end
    DATA = ReFit3means(DATA, ratio);
elseif strncmpi(name,'refitgm',6)
    DATA = ReFitGMDip(DATA);
elseif strcmp(name,'recalcxcorrcell')
    OptionMenu(DATA, [], name);
    return; %Don't set UserData below
elseif strncmpi(name,'setentry',7)
    id = get(gui,'value');
    DATA.id = id;
    Clusters = getappdata(DATA.toplevel,'Clusters');
    if id > length(Clusters) || isempty(Clusters{id})
    Clusters{id} = LoadClusters([DATA.name '/' DATA.strings{id}]);
    end
    figure(DATA.fig.clusters);
    DATA = doentry(DATA, Clusters, id);
elseif strncmpi(name,'followcorr',7)
    DATA.plot.alltype = 'followcorr';
    DATA.templateid = fliplr(DATA.currentpoint);
    DATA = PlotAllClusters(DATA,[]);
elseif strncmpi(name,'MakeCells',6)
        cellps = getappdata(DATA.toplevel,'CellPs');
        MakeCellId(cellps);
elseif strncmpi(name,'plotexpts',6)
        PlotExpts(DATA);
elseif strncmpi(name,'plottest',6) %for testing new plots
    OptionMenu(DATA,[],'xcorrallprobes');
    return;
    
    OptionMenu(DATA,[],'findconnect');
    eid = DATA.currentpoint(1);
    probes = find(DATA.selectprobe(eid,:));
    SetFigure(DATA,DATA.tag.spikes,'front');
    PlotSyncSpikes(DATA, eid, fliplr(probes),[1 1], varargin{:});
    return;
    DATA.plotcells.showmahal = 2;
    DATA.plotcells.showfitdp = 2;
    PlotCellList(DATA, 'showfig');
elseif strncmpi(name,'plotxcorrprobes',12)
    SetFigure(DATA, DATA.tag.xcorr);
    PlotXcorrs(DATA, DATA.xcorrs, 1:length(DATA.exptids), 0);
elseif strncmpi(name,'plotxcorr',6)
    SetFigure(DATA, DATA.tag.xcorr);
    PlotXcorrs(DATA, DATA.xcorrs, 1:length(DATA.exptids),1);
elseif strncmpi(name,'popplot',7)
    if length(varargin) && ischar(varargin{1})
        DATA.plot.alltype = varargin{1};
    end
    DATA = PlotAllClusters(DATA,[]);
elseif strncmpi(name,'plotmenu',7)
    PlotMenu(DATA,[], varargin{1}, varargin{2});
elseif strncmpi(name,'setplot',7)
    if gui
        strs = get(gui,'string');
    DATA.plot.onetype = deblank(strs(get(gui,'value'),:));
    end
    DATA = doentry(DATA, Clusters, DATA.id);
elseif strncmpi(name,'showdata',7)
    ShowData(DATA, DATA.currentpoint(1),DATA.currentpoint(2));
end
if initcall 
    if isfield(DATA,'exptid') %Do pop plot if data loaded
        DATA = PlotAllClusters(DATA,[]);
    else %have not loaded the clusters yets
        if ~isfield(DATA,'name') && ischar(name)
            fprintf('%s is not a directory\n',name);
            return;
        end

        listname = [DATA.name '/ExptList.mat'];
        if exist(listname)
            load(listname);
            DATA.exptid = ExptList.exptid;
            DATA.expnames = ExptList.expnames;
        end
        DATA.voffset = 1:DATA.nprobes *2;
    end
    DATA.defaultlayoutfile = defaultlayout;
    DATA.defaultconfigfile = defaultconfig;
    DATA = ApplyLayout(DATA, defaultlayout);
    DATA= ApplyConfig(DATA, defaultconfig);
    DATA = LoadCellFile(DATA);
    DATA = LoadExtra(DATA,0);
    DATA.ArrayConfig = GetArrayConfig(DATA.name);
    if isfield(DATA.ArrayConfig,'X')
        DATA.nprobes = length(DATA.ArrayConfig.X);
    elseif size(DATA.CellList,2) > 0
        DATA.nprobes = size(DATA.CellList,2);
    end
    DATA = PlotCellList(DATA,'showfig');
    if ~isfield(DATA,'exptid') && isfield(DATA.CellDetails,'exptids')
        DATA.exptid = DATA.CellDetails.exptids;
%        DATA = LoadExpts(DATA);
        for j = 1:length(DATA.exptid)
            Clusters{j} = {};
        end
        setappdata(DATA.toplevel,'Clusters',Clusters);
    end
    if ~isfield(DATA,'exptid')
        DATA.nexpts= 0;
    else
        DATA.nexpts = length(DATA.exptid);
    end
    DATA.selectprobe = zeros(DATA.nexpts,DATA.nprobes);
    DATA.user = GetUserName;

end

if strcmp(runcmd,'checkrateseq')
    varargout{1} = CheckAllRateSequences(DATA);
    varargout{1} = CopyFields(varargout{1},DATA, {'datadir', 'exptid' 'errs'});
elseif strcmp(runcmd,'updatearray')
    Clusters = getappdata(DATA.toplevel,'Clusters')
    varargout{1} = UpdateArrayConfig(DATA.datadir,Clusters);
    varargout{1} = CopyFields(varargout{1},DATA, {'datadir', 'exptid' 'errs'});
end
SetGUI(DATA);
set(DATA.toplevel,'UserData',DATA);


function out = doentry(DATA, Clusters, id)

type = DATA.plot.onetype;

done = 0;
DATA.currentpoint(1) = id;
if DATA.plotallxy
    PlotAllProbeXY(DATA);
    done = 1;
end
if DATA.plotspks
    PlotAllProbe(DATA, 'allspks');
    done = 1;
end
if done
    PlotCellList(DATA);
    out = DATA;
    return;
end
figure(DATA.fig.clusters);
if strcmp(type,'cov(shape)')
    CompareProbesShape(DATA,id,1);
elseif strcmp(type,'xc(shape)')
    CompareProbesShape(DATA,id,2);
elseif strcmp(type,'xy')
    cres = PlotSpikeTimes(Clusters{id},'xy');
else
    cres = PlotSpikeTimes(Clusters{id},type,'hitfcn',@HitPoint);
end
if strcmp(type,'xcorr')
     SetFigure(DATA,DATA.tag.xcorr);
     imagesc(squeeze(cres.synci(:,:,1)));
     DATA.synic{id} = cres;
     xc = cres;
     ofile = [DATA.name '/' strrep(DATA.strings{id},'ClusterTimes','Xcorrs')];
     save(ofile,'xc');
end
out = DATA;

function DATA = SetDefaults(DATA)
    
    DATA.profiling = 0;
    DATA.renderer = 'OpenGL'; %default for large figures
    DATA.exptname = [];
    DATA.probesperfile = 0;
    DATA.allclustering = 0;
    DATA.version = 1.3;
    DATA.progname = ['PlotClusters' DATA.version];
    DATA.markexpts = 'none';
    DATA.colors = mycolors('spkcolors');
    DATA.elmousept.down = 0;
    DATA.markh = [];
    DATA.selecth = 0;
    DATA.modified = 0;
    DATA.show.watchallcellxy = 0;
    DATA.show.watchallcellspks = 0;
    DATA.show.watchexptallspks = 0;
    DATA.show.watchallcellmean = 0;
    DATA.show.cellsummary = 0;
    DATA.show.linecontextmenus = 1;
    DATA.show.allvpcs = 0;
    DATA.clusteroffset = 0;
    DATA.currenttrial = 1;
    DATA.currentcell = 1;
    DATA.spikesloaded = 0;
    DATA.selectexpts = [];
    DATA.proberange = [];
    DATA.selectprobe = []; % not zero, else won't make it a matrix later
    DATA.currentcell = 1;
    DATA.currentcluster = 1;  %cluster to display
    DATA.currentcutcluster = 1;  %cluster for cutting
    DATA.currentcellcluster = 1; %cluster of currently setelcted cell
DATA.usesavedcodes = 0;
DATA.useautoclusters = 0;
DATA.usealltrials = 0;
DATA.Comments = [];
DATA.elmousept.shape = -1;
DATA.cellcluster = 1;  %current cluster for defining cells
DATA.checkrate.ff = 1;
DATA.checkrate.slope = 0.25;
DATA.checkrate.diff = 0.25;
DATA.tag.onecluster = [DATA.tag.top 'One'];
DATA.tag.clusters = [DATA.tag.top 'Fig'];
DATA.tag.all = [DATA.tag.top 'All'];
DATA.tag.popscatter = [DATA.tag.top 'Scatter'];
DATA.tag.templatesrc = [DATA.tag.top 'Template'];
DATA.tag.spkmean = [DATA.tag.top 'MeanSpike'];
DATA.tag.spkmeanb = [DATA.tag.top 'MeanSpike2'];
DATA.tag.celllist = [DATA.tag.top 'CellList'];
DATA.tag.expt = [DATA.tag.top 'Expt'];
DATA.tag.allxy = [DATA.tag.top 'AllXY'];
DATA.tag.allmean = [DATA.tag.top 'AllMeans'];
DATA.tag.allexpt = [DATA.tag.top 'AllExpt'];
DATA.tag.xyplot = [DATA.tag.top 'XYplot'];
DATA.tag.xcorr = [DATA.tag.top 'xCorr'];
DATA.tag.xcorrpop = [DATA.tag.top 'XCpop'];
DATA.tag.hist = [DATA.tag.top 'Histogram'];
DATA.tag.misc = [DATA.tag.top 'Misc'];
DATA.tag.spikes = [DATA.tag.top 'Spikes'];
DATA.tag.allspikes = [DATA.tag.top 'AllSpikes'];
DATA.tag.comments = [DATA.tag.top 'Comments'];
DATA.tag.xyseq = [DATA.tag.top 'XYseq'];
DATA.tag.rateseq = [DATA.tag.top 'Rateseq'];
DATA.tag.eigenvectors = [DATA.tag.top 'eigenvectors'];
DATA.tag.onecell = [DATA.tag.top 'Onecell'];
DATA.tag.smoothed = [DATA.tag.top 'smoothed'];
DATA.tag.fitgauss = [DATA.tag.top 'FitGauss'];
DATA.tag.checktimes = [DATA.tag.top 'checktimes'];

DATA.rebuild = 0;
DATA.cellbackup = 0;
DATA.spoolspikes = 0;
DATA.showspkxy=0;
DATA.showspkmean=1;
DATA.plotmeantype='both';
DATA.refitgm = 0;
DATA.plothist = 0;
DATA.nprobes = 24;
DATA.ArrayConfig = [];
DATA.plot.xcorrlabeltop = 0;
DATA.plot.xcorrlabelside = 0;
DATA.plot.xcorrexclude = [];
DATA.plot.xcorrtype = 'xcorr';
DATA.plot.xcorrpoptype = 'Shape/Efficacy';
DATA.plot.synci = 0;
DATA.plot.dprimemin = 0;

DATA.plotallxy = 0;
DATA.plotxyseq = 0;
DATA.plottrighist = 0;
DATA.plotspks = 0;
DATA.plotspk.bytrial = 1;
DATA.plotspk.showexcluded = 0;
DATA.plotspk.showcellmeans = 0;
DATA.plotspk.onescale =  0;
DATA.plotcells.showmahal = 0;
DATA.plotcells.showfitdp = 0;
DATA.plot.gmcid = 0;
DATA.plot.xcmax = 200;
DATA.plot.synctmax = 2;
DATA.plot.markauto = 0;
DATA.plot.markmanual = 0;
DATA.plot.markquick = 0;
DATA.plot.markcopy = 0;
DATA.plot.markgood = 0;
DATA.plot.markbadprobe = 0;
DATA.plot.markgoodmu = 0;
DATA.plot.showgm = 0;
DATA.plot.markexclusions = 0;
DATA.plot.scaledensity = 0;
DATA.plot.xyseqdensity = 0;

DATA.id = 1;
DATA.exclude.offpeak = 0;
DATA.exclude.onpeak = 0;
DATA.exclude.noncell = 0;
DATA.colorscheme = 'Plain';
DATA.datatype = 1;
DATA.templatesrc = 1;
DATA.show.exptno = 1;
DATA.show.exptname = 1;
DATA.show.ed = 0;
DATA.SpaceTypeLabels = {'PCs',  'ADC', 'Tmplt', 'Var-E', 'undef', 'N-D'};
DATA.SpaceChars = 'PATVUN';
DATA.NDSpaceChars = 'PATT';
DATA.mahaltype = 1;
DATA.crit.mahal = 2;
DATA.crit.synci = 0.2;
DATA.crit.dropi = 1.5;
DATA.xccrit = 0.9;
DATA.CellList = [];
DATA.CellDetails = [];
DATA.CellChanges = [];
DATA.steptype = 1;
DATA.nprobespool = 1;
DATA.currentpoint = [1 1];
DATA.NewCut.exptid = 0;
DATA.NewCut.probe = 0;
DATA.NewCut.saved = 1;%set to 0 if change
DATA.plotexpttype = 'none';

DATA.plot.onetype = 'dips';
DATA.plot.density = 0;
DATA.plot.mahalcmax = 5;
DATA.plot.alltype = 'mahal';
DATA.plot.showboundary = 0;
DATA.plot.cellgrid = 5;

DATA.markcell.candidates = 1; %mark possible cells
DATA.markcell.dropi = DATA.crit.dropi; %show cells with low drop index
DATA.markcell.ellipses = 0; %mark where eillpses were used for clsuter
DATA.markcell.mahal = 0; %show cells with low mahal distance
DATA.markcell.readmethod = 0; %mark expts using old read method
DATA.markcell.quick = 0; %show incomplet quantifications
DATA.markcell.grid = 1; %show grid lines
DATA.markcell.tagged = 0; %cells with manual tags added
DATA.markcell.expnames = 0; %show names on cell plot
DATA.markcell.goodmu=0;
DATA.markcell.goodcell=0;
DATA.allvpcsmode = 'fromspikes';


function TagMenu(a, b, fcn)
    DATA = GetDataFromFig(a);
    DATA.currentprobe = DATA.currentpoint(2);
    DATA.exptno = DATA.exptid(DATA.currentpoint(1));
    set(DATA.toplevel,'UserData',DATA);
    if strcmp(fcn,'comment')
        PlotComments(DATA.datadir,'parent',DATA.toplevel);
    return;
    end
    if strcmp(fcn,'flipline')
        FlipLineCrit(DATA);
    return;
    end
    if sum(strcmp(fcn,{'poor isolation' 'dropping spikes' 'poor stability'}))
        PlotComments(DATA.datadir,DATA,'addhidden',fcn);
        return;
    end
    reason = strmatch(fcn,{'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear' 'print'});
    if isempty(reason)
        reason = NaN;
    end
    
    DATA = LoadComments(DATA,DATA.name);
    if reason == 10 %clear
        DATA.selectprobe = zeros(size(DATA.selectprobe));
        reason = 0;
    elseif reason == 11 %print
        [a,b] = find(DATA.tagged > 0);
        for j = 1:length(a)
            fprintf('E%.1fP%d %s\n',a(j),b(j),DATA.tagstrings{DATA.tagged(a(j),b(j))});
        end
        reason = 0;
    end
    if ~isfield(DATA,'selectprobe')
        DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    end
    if ~isfield(DATA,'tagged')
        DATA.tagged = zeros(size(DATA.selectprobe));
    end
    id = find(DATA.selectprobe > 0);
    if reason < 6
        DATA.tagged(id) =  reason;
        DATA.tagged(DATA.currentpoint(1),DATA.currentpoint(2)) =  reason;
    end
    [a,b] = ind2sub(size(DATA.selectprobe),id);
    if isempty(id)
        fprintf('No Probes/Expts Selected\n');
    end
    for j = 1:length(id) 
        if  reason== 6
            GetComment(DATA,a(j),b(j));
        elseif ismember(reason, [7 8 9])
            AddComment(DATA, [], fcn);
            DATA = get(DATA.toplevel,'UserData');
        elseif reason == 1
            DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),:) = NaN;
            PlotCellList(DATA);
        end
        fprintf('Tagging E%d,P%d\n',a(j),b(j));
    end   
    if strcmp(DATA.plot.alltype,'Tagged')
        PlotAllClusters(DATA,[]);
    end
    SaveComments(DATA, DATA.name);
    set(DATA.toplevel,'UserData',DATA);

function GuiMenu(a,b, type)
    DATA = GetDataFromFig(a);
    if strcmp(type,'setmahaltype')
        DATA.mahaltype = strmatch(GetPopStr(a),{'2D' 'ND' ''  '1D' '2Gauss'});
        DATA = CalcDistances(DATA);
    elseif strcmp(type,'setmahalcrit')
        DATA.mahalcrit = str2num(get(a,'string'));
    elseif strcmp(type,'CallAllVPcs')
        CallAllVPcs(DATA, (DATA.currentpoint(1)),DATA.currentpoint(2));
    elseif strcmp(type,'checkclusters')
        CheckClusters(DATA,'exclusions');
        CheckClusters(DATA,'nclusters');
        CheckClusters(DATA,'fitspace');
        CheckClusters(DATA,'fittimes');
        DATA = CheckClusterLineSign(DATA);
    end
    set(DATA.toplevel,'UserData',DATA);
    
    
function SetGUI(DATA)
    SetCheck('SpikeXY',DATA.showspkxy,DATA.toplevel);
    SetCheck('PlotTrigHist',DATA.plottrighist,DATA.toplevel);
    SetCheck('PlotHistogram',DATA.plothist,DATA.toplevel);
    SetCheck('PlotXYseq',DATA.plotxyseq,DATA.toplevel);
    SetCheck('PlotSpks',DATA.plotspks,DATA.toplevel);
    SetCheck('LabelEd',DATA.show.ed,DATA.toplevel);
    SetCheck('LabelExptName',DATA.show.exptname,DATA.toplevel);
    SetCheck('LabelExptno',DATA.show.exptno,DATA.toplevel);
    SetCheck('AllVPcs',DATA.show.allvpcs,DATA.toplevel);

    F = findobj('Tag',DATA.tag.celllist,'type','figure');
    f = fields(DATA.markcell);
    onoff = {'off' 'on'};
    for j = 1:length(f)
        hm = findobj(F,'Tag',['mark' f{j}],'type','uimenu');
        on = DATA.markcell.(f{j}) > 0;
        set(hm,'checked',onoff{1+on});
    end
    f = {'watchallcellspks' 'watchallcellmean' 'watchallcellxy' 'watchexptallspks'};
    for j = 1:length(f)
        hm = findobj('Tag', f{j},'type','uimenu');
        if ~isempty(hm) 
        on = DATA.show.(f{j}) > 0;
        set(hm,'checked',onoff{1+on});
        end
    end    
    SetMenuCheck(DATA.toplevel, 'allvpcsmode', DATA.allvpcsmode,'exclusive');
    
function [C, ok] = GetClusterInfo(C, cl)
        ok = 1;
  if cl > 1 
      if ~isfield(C,'next')
          ok = 0;
      elseif length(C.next) > cl-1 && isfield(C.next{cl-1},'mahal')
          C = C.next{cl-1};
      elseif cl >= length(C.next)-1 || ~isfield(C.next{cl-1},'mahal')
          ok = 0;
      end
  end
  if ~isfield(C,'mahal')
      ok = 0;
  end
  
function DATA = CalcDistances(DATA)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    nc = max(DATA.nclusters(:));
    tic;
     for e = 1:length(Clusters)
         for p = 1:length(Clusters{e})
             for c = 1:nc
                 d(e,p,c) = DistanceMeasure(Clusters{e}{p},c,DATA.mahaltype);
             end
         end
     end
    toc
    DATA.mdistance = d;
    
function d = DistanceMeasure(C, cl, type, varargin)
    [C, ok] = GetClusterInfo(C, cl);
    if ok == 0
        d = NaN;
    elseif type < 5
            d = C.mahal(type);
        elseif type == 5
            if isfield(C,'fitdprime')
                d  = C.fitdprime(1);
            else
                d = NaN;
            end
        end
        
        
function GetComment(DATA,a,b)
    f = GetString(DATA.tag.comments,DATA.toplevel,{@AddComment, []});
    str = sprintf('Ex%.1f (%d), Probe%d',DATA.exptid(a),a,b);
    it = findobj(f,'tag','CommentLabel');
    if length(it) == 1
        set(it,'string',str);
    end
  
function AddComment(a,b,str)
    DATA = GetDataFromFig(a);
    DATA = LoadComments(DATA, DATA.name);
    n = length(DATA.Comments) + 1;
    DATA.Comments(n).user = DATA.user;
    DATA.Comments(n).time = now;
    DATA.Comments(n).ex = DATA.currentpoint(1);
    DATA.Comments(n).p = DATA.currentpoint(2);
    DATA.Comments(n).exptno = DATA.exptid(DATA.currentpoint(1));
   if length(str)  == 0
       x = get(a,'UserData');
       if isempty(x)
        str = get(a,'string');
       else
        str = get(x,'string');
       end
       close(get(a,'parent'));
   end
    DATA.Comments(n).comment = str;
    SaveComments(DATA, DATA.name);
    set(DATA.toplevel,'UserData',DATA);
    

function PrintComments(DATA,e,p)
    if isempty(DATA.Comments)
        return;
    end
    
    exid = sum(ismember(cat(2,DATA.Comments.ex),e),2)';
    for j = 1:length(p)
        id = find(exid > 0 & [DATA.Comments.p] == p(j));
        for k = 1:length(id)
            fprintf('E%.0f P%d %s\n',DATA.Comments(id(k)).exptno,p(j),DATA.Comments(id(k)).comment);
        end
    end

function HitTrial(data,b, cell)
    DATA = GetDataFromFig(data);
    Expts = getappdata(DATA.toplevel,'Expts');
    f = gcf;
    D = get(f,'UserData');
    for j = 1:length(D.exptlist)
        exptime(j) = Expts{D.exptlist(j)}.Header.timeoffset;
        timeoffset(j) = Expts{D.exptlist(j)}.Header.timeoffset-D.timeadjust(j);
    end
    pos = get(gca,'currentpoint');
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    x = get(data,'xdata');
    y = get(data,'ydata');
    r = ((pos(1,1)-x)./diff(xl)).^2+((pos(1,2)-y)./diff(yl)).^2;
    [a,b] = min(r);
    pos = [x(b) y(b)];
    id = find(D.expstarts <= pos(1));
    eid = D.exptlist(id(end));
    D.currentexpt = id(end);
    start = (pos(1)-timeoffset(id(end))).*10000;
    [a, id] = min(abs(start-[Expts{eid}.Trials.TrialStart])); 
    L = squeeze(DATA.CellList(eid,:,:));
    [probe,cl] = find(L == cell);

    if sum(DATA.currentpoint == [eid probe]) < 2 %Changed Expt/probe
        newpoint = 1;
    else
        newpoint = 0;
    end

    DATA = ClearSelections(DATA, 0, [eid probe]);
    DATA.currentcellcluster = cl;
    DATA.currentcluster = cl; %so xy display is right
    DATA.currentcell = cell;
    set(DATA.toplevel,'UserData',DATA);
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    if bt == 3 %shift
        id = DATA.currenttrial:id;
    end
    t = title(sprintf('Expt%d %s Cell %d P%d',eid,DATA.expnames{eid},cell,probe,cl),'color',get(data,'color'));
    AddClusterString(DATA, t, cl);
    if newpoint
        PlotCellList(DATA);
        DATA = ConditionalPlotXY(DATA, [], 0);
    end
    Expt = PlotExptCounts(DATA, eid, probe, cl);
    [a, h] = SpoolSpikes(DATA,DATA.currentpoint,'ids',id);
    AddLineContextMenu(DATA, h, eid, probe);
    figure(f);
    set(f,'UserData',D);
    fprintf('TrialHit %.1f,%.1f E%dP%dCl%d\n',pos(1),pos(2),eid,probe,cl);

    

    
function SetRateseqPlot(a,b,cellid)
    DATA = GetDataFromFig(a);
    DATA.currentcell = cellid;
PlotMenu(DATA,b,'cells','rateseqone', cellid);        
    
function PlotMenu(a, b, fcn, type, varargin)
    DATA = GetDataFromFig(a);
    onoff = {'off' 'on'};
    
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if strcmp(fcn,'expt') && strcmp(type,'combine')
        Data.Expt = PlotCombinedExpt(DATA);
    elseif strcmp(fcn,'expt')
%       PlotAllCell(DATA, type); %surely don't need this
        DATA.plotexpttype = type;
        m = get(a,'Parent');
        c = get(m,'Children');
        for j = 1:length(c)
            set(c(j),'Checked','off');
        end
        set(a,'checked','on');
        PlotExptCounts(DATA,DATA.currentpoint(1),DATA.currentpoint(2),DATA.currentcluster);
    elseif strcmp(fcn,'spikes')
        if strmatch(type,'Spool')
            SpoolCurrentSpikes(a,b);
        end
    elseif strcmp(fcn,'xyseq')
        if strmatch(type,'density')
            
            Expts = getappdata(DATA.toplevel,'Expts');
            DATA.plot.xyseqdensity = ~DATA.plot.xyseqdensity;
            if DATA.plot.xyseqdensity
                PlotXYSequence(DATA, DATA.currentpoint,Expts{e},'density');
            else
                PlotXYSequence(DATA, DATA.currentpoint,Expts{e});
            end
        end
    elseif strcmp(fcn,'xcorr')
        if strmatch(type,'zoom')
            DATA.plot.xcmax = varargin{1};
        elseif strmatch(type,'spkzoom')
            DATA.plot.synctmax = varargin{1};
        elseif strmatch(type,'recalcallcells')
            DATA.plot.synctmax = varargin{1};
        elseif strmatch(type,'replot')
            PlotXcorrs(DATA,DATA.xcorrs,1:length(DATA.exptid),1);
        elseif strmatch(type,'labelcellleft')
            DATA.plot.xcorrlabelleft = ~DATA.plot.xcorrlabelleft;
            set(a,'checked',onoff{1+DATA.plot.xcorrlabelleft});
        elseif strmatch(type,'labelcelltop')
            DATA.plot.xcorrlabeltop = ~DATA.plot.xcorrlabeltop; 
            set(a,'checked',onoff{1+DATA.plot.xcorrlabeltop});
        elseif strmatch(type,'exclude')
            if length(DATA.plot.xcorrexclude) < varargin{1} || DATA.plot.xcorrexclude(varargin{1}) == 0
                DATA.plot.xcorrexclude(varargin{1}) = 1;
                set(a,'checked','on');
            else
                DATA.plot.xcorrexclude(varargin{1}) = 0;
                set(a,'checked','off');
            end
        end
    elseif strcmp(fcn,'clusters')
       if strmatch(type,{'spksxy'})
           PlotAllProbeXY(DATA,'select');
           DATA = PlotAllProbe(DATA, 'selectspks');
       end
    elseif strcmp(fcn,'selectexpts')
        SetCheckExclusive(a);
        DATA.markexpts = type;
        PlotCellList(DATA);
        PlotClusterRates(DATA, 'rateseqall');
%        PlotClusterRates(DATA, 'rateseqone'); Migh tneed to make this
%        depend
    elseif strcmp(fcn,'markexpts')
        if strcmp(type,'relist')
            sms = get(get(a,'parent'),'children');
            delete(sms);
            AddExptList(get(a,'parent'),'markexpts',DATA)
            return;
        end
        SetCheckExclusive(a);
        DATA.markexpts = type;
        PlotCellList(DATA);
    elseif strcmp(fcn,'probes')
       if strmatch(type,{'spoolall' 'spoolcells' 'allspks' 'allcellspks' 'spooleverything' 'spooleverycell' 'selectspks' 'allprobespks' 'AllProbeMean'})
            DATA = PlotAllProbe(DATA, type);
       elseif strmatch(type,{'exptall'})
            DATA = PlotAllProbe(DATA, 'allspks');
            PlotAllProbeXY(DATA);
       elseif strmatch(type,{'onescale'})
           DATA.plotspk.onescale = ~DATA.plotspk.onescale;
           set(a,'Checked',onoff{1+DATA.plotspk.onescale});
       elseif strmatch(type,{'probeall'})
            DATA = PlotAllProbe(DATA, 'allprobespks');
            DATA = PlotExptsProbe(DATA, 'AllprobeXY');
           
       elseif strmatch(type,{'allprobespks' 'alltemplatespks' 'AllprobeXY' 'AllTemplateXY'})
            DATA = PlotExptsProbe(DATA, type);
       elseif strmatch(type,{'AllMean' 'AllMeanIm' 'AllExptIm'})
           PlotAllProbeMean(DATA, type);
       elseif strmatch(type,{'AllExptMeanIm'})
           PlotAllExptProbeMean(DATA, type);
       elseif strncmp(type,'spoolcurrent',8)
           SpoolCurrentSpikes(a,b);
       elseif strncmp(type,'plotsync',8)
           eid = DATA.xcorrs(DATA.xcid).eid;
           probes = find(DATA.selectprobe(eid,:));
           if strcmp(type,'plotsyncb')
               PlotSyncSpikes(DATA, eid, DATA.xcorrs(DATA.xcid).cells, DATA.xcorrs(DATA.xcid).clnum);
           else
               PlotSyncSpikes(DATA, eid, fliplr(DATA.xcorrs(DATA.xcid).cells), fliplr(DATA.xcorrs(DATA.xcid).clnum));
           end
       elseif strcmp(type,'spoolsync')
           Clusters = getappdata(DATA.toplevel,'Clusters');
           eid = DATA.xcorrs(DATA.xcid).eid;
           probes = find(DATA.selectprobe(eid,:));
           AllSpikes = CheckAllSpikes(DATA, eid, DATA.xcorrs(DATA.xcid).cells);
           if isempty(AllSpikes)
               fprintf('Spikes have not been loaded\n')
               return;
           end
           SetFigure(DATA,DATA.tag.spikes);
           SpoolAllProbes(DATA, eid, AllSpikes(eid,DATA.xcorrs(DATA.xcid).cells), Clusters{eid});

       elseif strcmp(type,'SelectXY')
           PlotAllProbeXY(DATA,'select');
       elseif strcmp(type,'AllXY')
           PlotAllProbeXY(DATA);
       else
            PlotAllCell(DATA, type);
        end
    elseif strncmp(fcn,'allcluster',8)
        if sum(strcmp(type,{'markmanual' 'markauto' 'markcopy' 'markexclusions' 'markgood' 'markgoodmu'}))
            DATA.plot.(type) = ~DATA.plot.(type);
            set(a,'Checked',onoff{DATA.plot.(type)+1});
            DATA = PlotAllClusters(DATA,[]);
        end
            
    elseif strncmp(fcn,'cells',4)
        if strcmpi(type,'allplots') || strcmpi(type,'allnext')
            if strcmpi(type,'allnext')
                DATA.currentcell = DATA.currentcell+1;
                F = findobj('type','figure','tag',DATA.tag.celllist);
                it = findobj(F,'tag','CellNumberId');
                set(it,'value',DATA.currentcell);
            end
            PlotAllCellXY(DATA);
%            PlotAllCellMean(DATA,'lineonly');
            PlotCellRates(DATA,'rates');
            PlotAllCell(DATA, 'allspkswithmean');
            [a,b,c] = FindCell(DATA, DATA.currentcell);
            X = zeros(size(DATA.CellList,1),size(DATA.CellList,2));
            id = sub2ind(size(X),a,b);
            X(id) = 1;
            if strmatch(DATA.plot.alltype,{'Fit-mahal'})
                PlotAllClusters(DATA, [], 'includelist', X);
            end
        elseif strcmpi(type,'addcontext')
            DATA.show.linecontextmenus = ~DATA.show.linecontextmenus; 
            set(a,'checked',onoff{1+DATA.show.linecontextmenus});
        elseif strcmpi(type,'allexpttune')
            Expts = getappdata(DATA.toplevel,'Expts');
            cname = Expt2Name(Expts{e},'addsuff');
            eid = find(strcmp(cname,DATA.expnames));
            cellid = unique(DATA.CellList(eid,:,:));
            cellid = cellid(cellid > 0);
            for j = 1:length(cellid)
                id = find(DATA.CellList == cellid(j));
                [exid, cid, clid] = ind2sub(size(DATA.CellList),id);
                depths(j) = mean(cid);
            end
            [depths,b] = sort(depths);
            cellid = cellid(b);
            
            if DATA.profiling
                profile on;
            end
            [nr,nc] = Nsubplots(length(cellid));
            SetFigure(DATA,DATA.tag.allexpt);
            colors = mycolors;
            for j = 1:length(cellid)
                DATA.selectprobe= zeros(size(DATA.CellList,1),size(DATA.CellList,2));
                mysubplot(nr,nc,j);
                id = find(DATA.CellList == cellid(j));
                [exid, cid, clid] = ind2sub(size(DATA.CellList),id);
                [b, id] = ismember(eid,exid);
                icolor = find(b);
                id = id(id > 0);
                [a,b] = sort(exid(id));
                id = id(b);
                b = sub2ind(size(DATA.selectprobe),exid(id),cid(id));
                fprintf('Building Expt for cell %d - %d expts\n',cellid(j),length(b));
                DATA.selectprobe(b) = 1;
                hold off; 
                [a, AllExpt] = PlotSelectedExpts(DATA,'expts', exid(id), cid(id), 'cluster',clid(id),'colors',colors(icolor),'cell',cellid(j));
                title(sprintf('Cell %d P%.1f',cellid(j),depths(j)),'verticalalignment','middle');
                set(gca,'xtick',[]);
            end
            if DATA.profiling
                profile viewer;
            end

        elseif strcmpi(type,'allcelltune')
            Expts = getappdata(DATA.toplevel,'Expts');
            extypes = unique(DATA.exptnames);

            cname = Expt2Name(Expts{e},'addsuff');
            eid = find(strcmp(cname,DATA.expnames));
            cellid = DATA.currentcell;
            
            id = find(DATA.CellList == cellid);
            [exid, cid, clid] = ind2sub(size(DATA.CellList),id);
            [nr,nc] = Nsubplots(length(extypes));
            SetFigure(DATA,DATA.tag.allexpt);
            for j = 1:length(extypes)
                DATA.selectprobe= zeros(size(DATA.CellList,1),size(DATA.CellList,2));
                mysubplot(nr,nc,j);
                eid = find(strcmp(extypes{j},DATA.expnames));
                [b, id] = ismember(eid,exid);
                id = id(id > 0);
                b = sub2ind(size(DATA.selectprobe),exid(id),cid(id));
                DATA.selectprobe(b) = 1;
                hold off; 
                [a, AllExpt] = PlotSelectedExpts(DATA);
                title(sprintf('%s',extypes{j}),'verticalalignment','top');
            end
                
        elseif strcmpi(type,'driftestimate')
           DATA = CompareShapes(DATA,'exptimage');
           %Calls FitDriftMatrix
           set(DATA.toplevel,'UserData',DATA);
        elseif strmatch(type,{'rateseqall' 'rateseqone'})
            PlotClusterRates(DATA, type);
        elseif strcmpi(type,'allxy')
            PlotAllCellXY(DATA);
        elseif strcmpi(type,'watchallcellmean')
            DATA.show.watchallcellmean = ~DATA.show.watchallcellmean;
            set(a,'checked',onoff{1+DATA.show.watchallcellmean});
        elseif strcmpi(type,'watchallcellxy')
            DATA.show.watchallcellxy = ~DATA.show.watchallcellxy;
            set(a,'checked',onoff{1+DATA.show.watchallcellxy});
        elseif strcmpi(type,'watchallcellspks')
            DATA.show.watchallcellspks = ~DATA.show.watchallcellspks;
            set(a,'checked',onoff{1+DATA.show.watchallcellspks});
        elseif strcmpi(type,'summary')
            DATA.show.cellsummary = ~DATA.show.cellsummary;
            set(a,'checked',onoff{1+DATA.show.cellsummary});
        elseif strcmpi(type,'markdropi')
            if DATA.markcell.dropi > 0
                DATA.markcell.dropi = 0;
                set(a,'Checked','off');
            else
                DATA.markcell.dropi = DATA.crit.dropi;
                set(a,'Checked','on');
            end
            PlotCellList(DATA);
        elseif strcmpi(type,'markmahal')
            if DATA.markcell.mahal > 0
                DATA.markcell.mahal = 0;
                set(a,'Checked','off');
            else
                DATA.markcell.mahal = DATA.crit.mahal;
                set(a,'Checked','on');
            end
            PlotCellList(DATA);
        elseif strncmpi(type,'mark',4)
                f = type(5:end);
                DATA.markcell.(f) = ~DATA.markcell.(f) ;
                set(a,'Checked',onoff{DATA.markcell.(f)+1});
            PlotCellList(DATA);
        elseif strncmpi(type,'markellipses',10)
            DATA.markcell.ellipses = ~DATA.markcell.ellipses;
            set(a,'Checked',onoff{1+DATA.markcell.ellipses});
            PlotCellList(DATA);
        elseif strcmpi(type,'markcandidates')
            DATA.markcell.candidates = ~DATA.markcell.candidates;
            set(a,'Checked',onoff{1+DATA.markcell.candidates});
            PlotCellList(DATA);
        elseif strcmpi(type,'spkrate+xy')
            PlotCellRates(DATA,'both');
        elseif strcmpi(type,'spkrate')
            PlotCellRates(DATA,'rates');
        elseif strcmpi(type,'allmean')
            PlotAllCellMean(DATA,'meanlines');
        elseif strcmpi(type,'allmeanim')
            PlotAllCellMean(DATA,'imageonly');
        elseif strcmp(type,'trighist')
            PlotAllCell(DATA, type);
        elseif strcmp(type,'allspks')
            PlotAllCell(DATA, type);
        elseif strcmp(type,'spoolall')
            PlotAllProbe(DATA, 'spoolcells');
        elseif strcmp(type,'spoolone')
            PlotAllCellSpikes(DATA);
        elseif strcmp(type,'spoolone')
            PlotAllCell(DATA, type);
        elseif strncmpi(type,'shapesize',5)
            res = getappdata(DATA.toplevel,'CellSizes');
            if isempty(res) || 1
                C = getappdata(DATA.toplevel,'Clusters');
                res = CalcCellSize(C, 'CellList', DATA.CellList);
                setappdata(DATA.toplevel,'CellSizes',res);
            end
            F = SetFigure(DATA, 'CellSize');
            if strncmpi(type,'shapesize',8)
                CalcCellSize(res,'callback',@HitShapePlot);
            elseif strncmpi(type,'shapecell',8)
                CalcCellSize(res,'onepercell');
            end
            D = get(F,'UserData');
            D.parent = DATA.toplevel;
            set(F,'UserData',D);
            DATA.figs.cellshape = F;
        end
    elseif strncmp(fcn,'means',4)
        DATA.plotmeantype = type;
        C = getappdata(DATA.toplevel,'Clusters');
        SetFigure(DATA,DATA.tag.spkmean);
        set(a,'Checked','on');
        PlotMeanSpike(C{e}{p}, p,0,type);
    end
    
    set(DATA.toplevel,'UserData',DATA);

 function HitShapePlot(a,b,id, c)
     D = GetDataFromFig(a);
     DATA = get(D.parent,'UserData');
     DATA.currentpoint(1) = D.cells{id}.eid;
     DATA.currentpoint(2) = D.cells{id}.probe;
     ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
     
     
function PlotClusterRates(DATA, type,varargin)
    markexpts = DATA.markexpts;
    currentcell = DATA.currentcell;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cell',4)
            j = j+1;
            currentcell = varargin{j};
        elseif strncmpi(varargin{j},'expt',4)
            j = j+1;
            markexpts = varargin{j};
        end
        j = j+1;
    end
    
    Expts = getappdata(DATA.toplevel,'Expts');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    cellids = unique(DATA.CellList(:));
    cellids = cellids(cellids > 0);
    np = DATA.nprobes;
    f = SetFigure(DATA,DATA.tag.rateseq);
    X = get(f,'UserData');
    hold off;
    colors = mycolors;
    AllExpts = [];
    AllProbes = [];
    Allclid = [];
    
    ts = now;
    if strcmp(type,'rateseqone')
        cellids = currentcell;
        scaling = 'rates';
    else
        scaling = 'normalize';
    end
    if strcmp(markexpts,'none')
        exid = 1:DATA.nexpts;
        exstr = 'All Expts';
    else
        exid = find(strcmp(markexpts, DATA.expnames));
        exstr = markexpts;
    end
    
%If have not loaded allexpts yet, could either force a load here, or only plot
%what is loaded.  For now, only plot what is loaded.
badex = [];
for j = 1:length(exid)
    if exid(j) > length(Expts) || isempty(Expts{j})
       badex = [badex exid(j)];
    end
end
exid = setdiff(exid,badex);

if length(badex)
    exstr = sprintf('%s (Loaded %s)',exstr,sprintf('%d ',exid));
end

gaps(1) = Expts{exid(1)}.Header.timeoffset;
    Expts{exid(1)}.Header.timeadjust = gaps(1);
    for j = 2:length(exid)
        durs(j) = (Expts{exid(j)}.Trials(end).TrialStart-Expts{exid(j)}.Trials(1).TrialStart)./10000;
        xdur(j) = Expts{exid(j)}.Header.timeoffset-Expts{exid(j-1)}.Header.timeoffset;
        gaps(j) = ((Expts{exid(j)}.Trials(1).TrialStart -Expts{exid(j-1)}.Trials(end).TrialStart)./10000)...
            +Expts{exid(j)}.Header.timeoffset-Expts{exid(j-1)}.Header.timeoffset -2;
%makes a mess of gtie expts (one .mat file). Its calculation above anyway -
%? need it for uProbe expts?
%            -Expts{exid(j-1)}.Trials(1).TrialStart/10000;
        Expts{exid(j)}.Header.timeadjust = sum(gaps);
        starts(j) = Expts{exid(j)}.Header.timeoffset-Expts{exid(j)}.Header.timeadjust ...
            +Expts{exid(j)}.Trials(1).TrialStart/10000;
    end
    if length(exid) == 1
        X.timeadjust = 0;
        X.expstarts = 0;
    else
        X.expstarts = starts-2;  %time of first trial
        X.timeadjust = cumsum(gaps);
    end
    X.currentexpt = exid(1);
    np = size(DATA.CellList,2); %this is the important np for calc below
    for j = 1:length(cellids)
        aExpts = {};
        [a,b] = find(DATA.CellList == cellids(j));
        id = find(ismember(a,exid));
        a = a(id);
        b = b(id);
        c = ceil(b./np);
        p = 1+mod(b-1,np);
        for k = 1:length(a)
            [xcl, b] = FindExcludedTrials(DATA,a(k),p(k),c(k),Clusters{a(k)}{p(k)});
            aExpts{k} = CountSpikes(Expts{a(k),1},Clusters{a(k)}{p(k)},c(k),b.xid);
            aExpts{k}.Header.cellnumber  = cellids(j);
        end
        AllExpts = cat(1,AllExpts,a);
        AllProbes = cat(1,AllProbes,p);
        Allclid = cat(1,Allclid,c);
        if length(a)
            a = PlotRateSequence(aExpts,'color',colors{j},scaling,'offset',(j-1)*2,'bytime','callback',{@HitTrial, cellids(j)});
            text(max(a.times)*1.01,a.meanrates(end),num2str(cellids(j)),'color',colors{j},'horizontalalignment','left');
            hold on;
            X.RateRange(cellids(j),:) = minmax(a.rates);
            X.handles(cellids(j)) = a.handle;
        end
    end
    axis('tight');
    X.yrange = get(gca,'ylim');
    X.xrange = get(gca,'xlim');
    X.exptlist = exid;
    X.cellids = cellids;
    set(f,'UserData',X);
    fprintf('Took %.1f\n',mytoc(ts));
    set(gcf,'keypressfcn',@RateSeqKeyPressed);
    yl = get(gca,'ylim');
    extypes = {};
    AllExpts = unique(AllExpts);
    for j = 1:length(AllExpts)
        extypes{j} = Expts{AllExpts(j)}.Header.expname;
    end
    types = unique(extypes);
    first = zeros(size(types));
    for j = 1:length(AllExpts)
        x = find(X.exptlist == AllExpts(j));
        id = strmatch(Expts{AllExpts(j)}.Header.expname,types,'exact');
        t = Expts{AllExpts(j)}.Header.timeoffset-Expts{AllExpts(j)}.Header.timeadjust;
        t = X.expstarts(x);
        h = plot([t t],yl,'--','color',colors{id});
        if length(cellids) == 1
            dstr = sprintf('M%.1f D%.1f',DistanceMeasure(Clusters{AllExpts(j)}{AllProbes(j)},Allclid(j),DATA.mahaltype),...
                DATA.dropi(AllExpts(j),AllProbes(j),Allclid(j)));
        else
        dstr = [];
        end
        if ~first(id)
            set(h,'linestyle','-');
            h = text(t,yl(2),sprintf('%d %s %s',AllExpts(j),types{id},dstr),'color',colors{id},'VerticalAlignment','top','HorizontalAlignment','right','rotation',90);
            first(id) =1;
        elseif length(types) == 1
            h = text(t,yl(2),sprintf('%d %s',AllExpts(j),dstr),'color','k','VerticalAlignment','top','HorizontalAlignment','right','rotation',90);
        end
    end
    if isempty(AllExpts)
        title(sprintf('Cell %d No Data %s',cellids,exstr));
    elseif length(cellids) == 1
        title(sprintf('Cell %d %s',cellids,exstr));
        check = CheckExptRates(aExpts);
        exid = [check.errs.exptno];
        [a,b] = max(abs([check.errs.ff]));
        [c,d] = max(abs(log(check.diffs)));
        c = log(check.diffs(d));
        [e,f] = max(abs([check.errs.slope]));
        e = check.errs(f).slope;
        fprintf('ExptV/M %.2f. Max diff %.2f (%d), max FF %.2f(%d), slope %.2f (%d)\n',...
            check.blkff,c,exid(d),a,exid(b),e,exid(f));
    end
    cmenu = uicontextmenu;
    cellid = unique(DATA.CellList);
    cellid = cellid(cellid>0);
    for j = 1:length(cellid)
        uimenu(cmenu,'label',sprintf('Cell %d',cellid(j)),'Callback',{@SetRateseqPlot, cellid(j)});
    end
    set(gca,'UIContextMenu',cmenu);

    
function DATA = CalcCellMeans(DATA)

    cellid = unique(DATA.CellList(:));
    cellid = cellid(cellid > 0);
    Clusters = getappdata(DATA.toplevel,'Clusters');
    for j = 1:length(cellid)
        c = cellid(j);
        id = find(DATA.CellList==c);
        [eid,pid,cl] = ind2sub(size(DATA.CellList),id);
        C = Clusters{eid(1)}{pid(1)};
        ms = zeros(length(C.chspk),size(C.MeanSpike.ms,2));
        chspk = [];
        for e = 1:length(eid)
            ichspk = Clusters{eid(e)}{pid(e)}.chspk;
            ms = ms + Clusters{eid(e)}{pid(e)}.MeanSpike.ms(ichspk,:);
            chspk(e,:) = ichspk;
        end
        DATA.CellDetails.MeanSpike{c}.ms = ms./length(e);
        DATA.CellDetails.MeanSpike{c}.chspk = mean(chspk);
    end
    SaveCellList(DATA);
    
function E = CountSpikes(Expt, C, clnum,xcl)
%E = CountSpikes(Expt, C, clnum,xcl)
%finds spikes in C that match trials in Expt, and
%adds them to the Expt.Trials strucutre.
%clnum 1 = cluster 1, i.e. C.clst == 2
        latency = 500;
        id = find(C.clst == clnum+1);
        t = (C.times(id) .*10000);
        for j = 1:length(Expt.Trials)
            starts(j) =  Expt.Trials(j).Start(1);
            ends(j) =  Expt.Trials(j).End(end);
        end
        if isempty(id)
            npre = 0;
            npost = 0;
        else
        npre = find(starts < t(1));
        if npre > 10
            id = find(starts < t(1)+10000);
            Expt.Trials = Expt.Trials(id);
        end
        npost = find(starts > t(end));
        if npost > 1
            id = find(starts < t(end));
            Expt.Trials = Expt.Trials(id);
        end
        end
        for j = 1:length(Expt.Trials)
            T = Expt.Trials(j);
            id = find(t > T.Start(1)+latency & t <= T.End(end)+latency);
            Expt.Trials(j).count = length(t(id));
            Expt.Trials(j).Spikes = round([t(id)-T.Start(1)]');
            if size(t,2) == 1
                Expt.Trials(j).Spikes = Expt.Trials(j).Spikes';
            end
        end
        E = Expt;
        [a,b] = setdiff([Expt.Trials.id],xcl);
        E.Trials = Expt.Trials(b);
        
    function PlotCellRates(DATA,type)
    
        Expts = getappdata(DATA.toplevel,'Expts');
        Clusters = getappdata(DATA.toplevel,'Clusters');
        [e, p, c] = FindCell(DATA, DATA.currentcell);
        id = sub2ind(size(Expts),e,p);
        SetFigure(DATA,DATA.tag.xyseq);
        for j = 1:length(id)
            xcl = FindExcludedTrials(DATA,e(j),p(j),c(j),Clusters{e(j)}{p(j)});
            xid = [Expts{id(j)}.Trials(xcl).id];
            Expts{id(j)} = CountSpikes(Expts{id(j)},Clusters{e(j)}{p(j)},c(j),xid);
        end
        plottype = strmatch(type, {'both' 'rates' 'xy'});
        if ismember(plottype,[1 2])
            if plottype == 1
                subplot(2,1,1);
                hold off;
            else
                subplot(1,1,1);
                hold off;
            end
            PlotRateSequence(Expts(id));
        end
        if ismember(plottype,[1 3])
            if plottype == 1
                subplot(2,1,2);
            else
                subplot(1,1,1);
            end
            hold off;
            for j = 1:length(id)
                PlotXYSequence(DATA, Clusters{e(j)}{p(j)}, 'expt', Expts{id(j)});
                hold on;
            end
        end
        
function MarkExpts(DATA,type)
    
    if strcmp(type,'none')
        return;
    end
    id = find(strcmp(type,DATA.expnames));
    for j = 1:length(id)
        DrawBox(id(j),[1.3 DATA.nprobes-0.7],3,'color','g','linewidth',1);
    end

function dup = isduplicate(DATA, row, p, cl)
if size(DATA.CellList,3) < cl || DATA.CellList(row,p,cl) >= 0
    dup = 0;
else
    dup = -DATA.CellList(row,p,cl);
end
    
function [true, cellid, clid] = isacell(DATA, row, p, clid)
        cellid = 0;
        if nargin > 3
            if clid <= size(DATA.CellList,3)
                cellid = DATA.CellList(row,p,clid);
            else
                cellid = 0;
            end
            true = cellid > 0;
            return;
        else
        clid = 0;
        end
        true = sum(DATA.CellList(row,p,:) > 0);
        dup = DATA.CellList(row,p,:) < 0;
        if true
            cellid = DATA.CellList(row,p,:);
            clid = find(cellid > 0);
            cellid = cellid(cellid>0);
        end
        
        
function SavePlotClustersConfig(DATA,file, varargin);
    savefields = {'showspkxy'  'plotallxy' 'plotspks' 'plothist'...
        'plotcells' 'plotxyseq' 'plot' 'plottrighist', 'comparerates' 'show' 'renderer'};
    SaveConfig(DATA, file, savefields, varargin{:});
    
      
function DATA = ApplyConfig(DATA, varargin)
    %
    configfile = [];
    if isfield(DATA,'configfile')
        configfile = DATA.configfile;
    end
    j = 1;
    while j <= length(varargin)
        if j == 1 && ischar(varargin{j})
            configfile = varargin{j};
        end
        j = j+1;
    end

    if isempty(configfile)
        return;
    end
    DATA.configfile = configfile;

    DATA = ReadConfig(DATA, DATA.configfile,'print','nochoose');
    if isfield(DATA,'layoutfile') && exist(DATA.layoutfile)
        ApplyLayout(DATA);
    end

    set(DATA.toplevel,'UserData',DATA);
    

function DATA = ApplyLayout(DATA, varargin)
layoutfile = [];
if isfield(DATA,'layoutfile')
    layoutfile = DATA.layoutfile;
end
j = 1;
while j <= length(varargin)
    if j == 1 && ischar(varargin{j})
        layoutfile = varargin{j};
    end
    j = j+1;
end

if ~exist(layoutfile)
        fprintf('Cant read %s\n',layoutfile);
        return;
end
DATA.layoutfile = layoutfile;
load(DATA.layoutfile);
    setappdata(DATA.toplevel,'Figpos',Figpos);
    f = fields(Figpos);
    for j = 1:length(f)
        it = findobj('type','figure','tag',f{j});
        if length(it) == 1
                set(it,'Position',Figpos.(f{j}));
        end
    end
    if exist('Showvals','var')
        f = fields(Showvals);
        for j = 1:length(f)
            DATA.show.(f{j}) = Showvals.(f{j});
        end
    end
    if exist('Datavals','var')
        f = fields(Datavals);
        for j = 1:length(f)
            DATA.(f{j}) = Datavals.(f{j});
        end
    end

    function OptionMenu(a, b, fcn)
        
        DATA = GetDataFromFig(a);
        onoff = {'off' 'on'};
        if strcmp(fcn,'allvpcsmode')
            tag = get(a,'Tag');
            DATA.allvpcsmode = tag;
            SetMenuCheck(a,'exclusive');
        elseif strcmp(fcn,'usesavedcodes')
            DATA.usesavedcodes = ~DATA.usesavedcodes;
            set(a,'Checked', onoff{1+DATA.usesavedcodes});
        elseif strcmp(fcn,'useautoclusters')
            DATA.useautoclusters = ~DATA.useautoclusters;
            set(a,'Checked', onoff{1+DATA.useautoclusters});
            ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2),'oneprobe');
        elseif strcmp(fcn,'popplot')
                DATA = PlotAllClusters(DATA,[]);
        elseif strcmp(fcn,'tofront')
            FiguresToFront(DATA.tag);
        elseif strcmp(fcn,'calccellmeans')
            DATA = CalcCellMeans(DATA);
        elseif strcmp(fcn,'checkratesequence')
            CheckAllRateSequences(DATA);
        elseif strcmp(fcn,'fillcellsfrommark')
            DATA = FillCellList(DATA, 'frommark');
        elseif strcmp(fcn,'loadconfig')
            if ~isfield(DATA,'configfile') || isempty(DATA.configfile)
                DATA.configfile =  '/bgc/group/matlab/preferences/PlotClusters/default.config';
            end
            [afile, pathname] = uigetfile(DATA.configfile);
            DATA.configfile = [pathname afile];
            DATA = ApplyConfig(DATA);
        elseif strcmp(fcn,'optimizeall')
            Clusters = getappdata(DATA.toplevel,'Clusters');
            ax  = findobj(gcf,'type','axes');
            oldname = get(gcf,'name');
            for j = 1:length(ax)
                axdata = get(ax(j),'UserData');
                axes(ax(j));
                set(gcf,'Name',sprintf('Optimizing E%dP%d',axdata.eid,axdata.probe));
                FixBoundary(DATA, Clusters, axdata.eid,axdata.probe);
            end
            set(gcf,'Name',oldname);
        elseif sum(strcmp(fcn,{'saveconfig' 'savedefaultconfig'}))
            if strcmp(fcn,'savedefaultconfig') || ~isfield(DATA,'configfile')
                DATA.configfile = DATA.defaultconfigfile;
            end
            [outname, pathname] = uiputfile(DATA.configfile);
            if outname
                DATA.configfile = [pathname outname];
                SavePlotClustersConfig(DATA, DATA.configfile);
                fprintf('Setttings saved to %s\n',DATA.configfile);
            end
        elseif strcmp(fcn,'scaledensity')
            DATA.plot.scaledensity = ~DATA.plot.scaledensity;
            set(a,'checked',onoff{1+DATA.plot.scaledensity});
        elseif strcmp(fcn,'loadlayout')
            if ~isfield(DATA,'layoutfile') || isempty(DATA.layoutfile)
                DATA.layoutfile =  '/bgc/group/matlab/preferences/PlotClusters/Bruce.layout.mat';
            end
            [afile, pathname] = uigetfile(DATA.layoutfile);
            DATA.layoutfile = [pathname afile];
            ApplyLayout(DATA);
        elseif sum(strcmp(fcn,{'savelayout' 'savedefaultlayout'}))
            Figpos = getappdata(DATA.toplevel,'Figpos');
            if strcmp(fcn, 'savedefaultlayout')
                DATA.layoutfile = DATA.defaultlayoutfile;
            elseif ~isfield(DATA,'layoutfile') || isempty(DATA.layoutfile)
                DATA.layoutfile =  '/bgc/group/matlab/preferences/PlotClusters/Bruce.layout.mat';
            end
            f = fields(DATA.tag);
            for j = 1:length(f);
                it = findobj('type','figure','Tag', DATA.tag.(f{j}));
                if length(it) == 1
                    Figpos.(DATA.tag.(f{j})) = get(it,'Position');
                end
            end
            Showvals = DATA.show;
            f = {'showspkxy' 'showspkmean'};
            for j = 1:length(f)
            Datavals.(f{j}) = DATA.(f{j});
            end
            [outname, pathname] = uiputfile(DATA.layoutfile);
            if outname
                DATA.layoutfile = [pathname outname];
                save(DATA.layoutfile,'Figpos','Showvals','Datavals');
                fprintf('Layout saved to %s\n',DATA.layoutfile);
            end

            SetCheckExclusive(a);
        
        elseif fcn == 1
            DATA.plot.density = ~DATA.plot.density;
            set(a,'Checked',onoff{DATA.plot.density+1});
            D = getappdata(gcf,'allplot');
            if isfield(D,'type') %called from AllXY plot
                if strcmp(D.type,'AllprobeXY')
                    PlotExptsProbe(DATA, D.type);
                end
            elseif isfield(DATA,'currentpoint')
                ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
            end
        elseif fcn == 2
            reload = 1;
            DATA = PlotCellList(DATA,'showfig','reload');
        elseif fcn == 3
            PlotExpts(DATA);
        elseif fcn == 4
           DATA =  PlotAllProbeXY(DATA);
        elseif fcn == 5
            PlotAllProbeMean(DATA,'lineonly');
        elseif fcn == 6
            PlotAllProbeMean(DATA,'imageonly');
        elseif fcn == 7
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),1);
        elseif strcmp(fcn,'excludecurrent')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),DATA.currentcellcluster);
        elseif strcmp(fcn,'excludecl4')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),4);
        elseif strcmp(fcn,'excludecl3')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),3);
        elseif strcmp(fcn,'excludecl2')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),2);
        elseif strcmp(fcn,'usealltrialscl1')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),1, 'reset');
        elseif strcmp(fcn,'setallprobeplot')
            
        elseif strcmp(fcn,'findconnect')
            for j =1:length(DATA.xcorrs)
                xc = DATA.xcorrs(j).xc;
                [a,b] = max(xc);
                p = prctile(xc,[50 90]);
                t = DATA.xcorrval.times(b);
                if  a > p(2) * 1.1 && a > p(1) +std(xc)*4 && abs(t) < 0.01 && abs(DATA.xcorrs(j).separation) > 2
                    cscore(j) = 1;
                else
                    cscore(j) = 0;
                end
                delays(j) = t;
            end
            np = ceil(sqrt(sum(cscore > 0)));
            id = find(cscore > 0);
            for j = 1:length(id)
                mysubplot(np,np,j)
                plot(DATA.xcorrval.times,DATA.xcorrs(id(j)).xc,'k-');
                set(gca,'xtick',[],'ytick',[]);
                axis('tight');
                set(gca,'buttondownfcn',{@HitXcorr,'zoom',0,0})
                yl = get(gca,'ylim');
                xl = get(gca,'xlim');
                line([0 0],yl,'color','r');
                text(xl(2),yl(2),sprintf('E%dP%d,%d %.1fms',DATA.xcorrs(id(j)).eid,DATA.xcorrs(id(j)).probes(1),DATA.xcorrs(id(j)).probes(2),delays(id(j)).*1000),...
                    'verticalalignment','top','horizontalalignment','right');
            end
        elseif strcmp(fcn,'oldxcorrcells')
            recalc = 0;
            C= getappdata(DATA.toplevel,'Clusters');
            Clusters = C{DATA.currentpoint(1)};
            e = DATA.currentpoint(1);
            cid = [];
            cnum = [];
            for j = 1:size(DATA.CellList,3)
                ids{j} = find(DATA.CellList(e,:,j) > 0);
                idn = ones(size(ids{j})) .*j;
                cid = [cid ids{j}];
                cnum = [cnum idn+1];
            end
            SetFigure(DATA,DATA.tag.xcorr,'front');
            np = length(cid);
            if isfield(DATA,'xcorrs')
                probes = [DATA.xcorrs.probes];
                cellid = [DATA.xcorrs.cells];
            else
                probes = [];
                cellid = [];
            end
            tic;
            DATA.xcorrval.times = [-0.2:0.001:0.2];
            for j = 1:np
                for k = 1:j
                    P = Clusters{cid(j)};
                    Q = Clusters{cid(k)};
                    if isempty(cellid)
                        id = [];
                    else
                    id = find(cellid(:,1) == cid(j) && cellid(:,2) == cid(k) || cellid(:,2) == cid(j) && cellid(:,1) == cid(k));
                    end
                    if length(id) == 1 && DATA.xcorrs(id).calctime > P.savetime(1)
                        xc = DATA.xcorrs(id).xc;
                    else
                        pid = find(P.clst == cnum(j));
                        qid = find(Q.clst == cnum(j));
                        xc = xcorrtimes(P.times(pid),Q.times(qid));
                        DATA.xcorrs(id).xc = xc;
                    end
                    if k == j
                        xc(201) = 0;
                    end
                    mysubplot(np,np,k+(j-1)*np,'leftmargin',0.02);
                    plot(-200:200,xc,'k');
                    axis('tight');
                    if k == j
                        title(sprintf('P%d',cid(k)));
                    end
                    if k == 1
                        ylabel(sprintf('P%d',cid(j)));
                    end
                    set(gca,'xtick',[],'ytick',[]);
                end
            end
            PlotXcorrs(DATA, DATA.xcorrs, DATA.currentpoint(1), bycell);
            toc
        elseif find(strcmp(fcn,{'allclusterxcorr'}))
            C= getappdata(DATA.toplevel,'Clusters');
            Clusters = C{DATA.currentpoint(1)};
            nc = 1;
            for k = 1:length(Clusters)
                if DATA.plot.dprimemin == 0 || (Clusters{k}.fitdprime(1) < DATA.plot.dprimemin)
                    cells(nc).p = k;
                    cells(nc).cl = 1;
                    if isfield(Clusters{k},'clst')
                        t = find(Clusters{k}.clst == 2);
                        Clusters{k}.times = Clusters{k}.times(t);
                    end
                    nc = nc+1;
                end
                for j = 1:length(Clusters{k}.next)
                    if isfield(Clusters{k}.next{j},'times') && ...
                            (DATA.plot.dprimemin == 0 || Clusters{k}.next{j}.fitdprime(1) < DATA.plot.dprimemin)
                        cells(nc).p = k;
                        cells(nc).cl = 1+j;
                        nc = nc+1;
                    end
                end
            end
            SetFigure(DATA,DATA.tag.xcorr);
            PlotAllXcorr(DATA, Clusters,cells,'callback',@PlotXcorr);
            ReplotXcorrs(DATA,[],DATA.plot.xcorrpoptype);
        elseif find(strcmp(fcn,{'xcorrcells' 'xcorrallcells' 'xcorrallprobes' 'recalcxcorrcell' 'recalcxcorrall'}))
            if  strcmp(fcn, 'recalcxcorrcell')
                recalc = 1;
                DATA = rmfields(DATA,'xcorrs');
                bycell = 1;
            elseif  strcmp(fcn, 'recalcxcorrall')
                recalc = 1;
                DATA = rmfield(DATA,'xcorrs');
                bycell = 0;
            else
                recalc = 0;
            end
            if  strcmp(fcn, 'xcorrallprobes')
                bycell = 0
            else
                bycell = 1;
            end
            C= getappdata(DATA.toplevel,'Clusters');
            Clusters = C{DATA.currentpoint(1)};
            if strcmp(fcn,'xcorrcells')
                expts = DATA.currentpoint(1);
            else
                expts = 1:length(DATA.exptid);
            end
            e = 1;
            cid = [];
            cnum = [];
            for j = 1:size(DATA.CellList,3)
                ids{j} = find(DATA.CellList(e,:,j) > 0);
                idn = ones(size(ids{j})) .*j;
                cid = [cid ids{j}];
                cnum = [cnum idn+1];
            end
            SetFigure(DATA,DATA.tag.xcorr,'front');
            ClearPlot;
            nc = 0;
            if isfield(DATA,'xcorrs') && recalc == 0
                probes = cat(1,DATA.xcorrs.probes);
                cellid = cat(1,DATA.xcorrs.cells);
                [DATA.xcorrs.valid] = deal(0);
                for j = 1:length(DATA.xcorrs)
                    e = DATA.xcorrs(j).eid;
                    p = DATA.xcorrs(j).probes;
                    c = DATA.xcorrs(j).cells;
                    cl = DATA.xcorrs(j).clnum;
                    cc(1) = DATA.CellList(e,p(1),cl(1));
                    if DATA.CellList(e,p(1),cl(1)) ~= c(1)
                        fprintf('E%dP%d Setting old xcorr C%d to C%d\n',e,p(1),c(1),DATA.CellList(e,p(1),cl(1)));
                        DATA.xcorrs(j).cells(1) = DATA.CellList(e,p(1),cl(1));
                        nc = nc+1;
                    end
                    if DATA.CellList(e,p(2),cl(2)) ~= c(2)
                        fprintf('E%dP%d Setting old xcorr  C%d to C%d\n',e,p(2),c(2),DATA.CellList(e,p(2),cl(2)));
                        DATA.xcorrs(j).cells(2) = DATA.CellList(e,p(2),cl(2));
                        nc = nc+1;
                    end
                end

            else
                probes = [0 0];
                cellid = [0 0];
            end
            np = length(cid);
            tic;
            nxc = 1;
            DATA.xcorrval.times = [-0.2:0.001:0.2];
            for eid = 1:length(expts)
                e = expts(eid);
                cid = [];
                cnum = [];
                Clusters = C{e};
                cells = [];
                for j = 1:size(DATA.CellList,3)
                    ids{j} = find(DATA.CellList(e,:,j) > 0);
                    idn = ones(size(ids{j})) .*j;
                    cid = [cid ids{j}];
                    cnum = [cnum idn+1];
                    cells = [cells DATA.CellList(e,ids{j},j)];
                end
                [cells, id] = sort(cells);
%                id = id(ismember(cells,[4 14]));
                cid = cid(id);
                cnum = cnum(id);
                if bycell
                    np = length(cid);
                else
                    np = DATA.nprobes;
                    cnum = ones(1,np).*2;
                end
                
                
                fprintf('E%d %d Cells:%s\n',e,np,sprintf(' %d',cells));
                for j = 1:np
                    for k = 1:j
                        if bycell
                            P = Clusters{cid(j)};
                            Q = Clusters{cid(k)};
                            id = find(cellid(:,1) ==cells(j) & cellid(:,2) == cells(k) | cellid(:,2) == cells(j) & cellid(:,1) == cells(k));
                        else
                            P = Clusters{j};
                            Q = Clusters{k};
                            id = find(probes(:,1) == j & probes(:,2) == k);
                        end
                        if ~isfield(DATA,'xcorrs')
                            id = [];
                            DATA.xcorrs = [];
                        else
                            id = intersect(id,find([DATA.xcorrs.eid] == e));
                        end
                        if length(id) == 1
                            p = DATA.xcorrs(id).probes;
                            c = DATA.xcorrs(id).cells;
                            cl = DATA.xcorrs(id).clnum;
                            cc(1) = DATA.CellList(e, p(1), cl(1));
                            cc(2) = DATA.CellList(e, p(2), cl(2));
                            if sum(cc == c) == 2
                                ok = 1;
                            else
                                ok = 0;
                            end
                        else
                            ok = 0;
                        end
                        if recalc == 0 && ok && DATA.xcorrs(id).calctime > max([P.savetime(1) Q.savetime(1)])
                            xc = DATA.xcorrs(id).xc;
                            DATA.xcorrs(id).valid = 1;
                        else
                            if length(id) >= 1
                                nxc = id(1);
                            else
                                nxc = length(DATA.xcorrs)+1;
                            end
                            if recalc == 0
                                fprintf('Recalculating %d->%d',cells(j),cells(k));
                            end
                            if cells(k) > cells(j)
                                fprintf('Order error\n');
                            end
                            if DATA.CellList(e,cid(j),cnum(j)-1) ~= cells(j)
                                fprintf('Cell %d ->cell %d in list\n',DATA.CellList(e,cid(j),cnum(j)-1),cells(j));
                            end
                            if DATA.CellList(e,cid(k),cnum(k)-1) ~= cells(k)
                                fprintf('Cell %d ->cell %d in list\n',DATA.CellList(e,cid(k),cnum(k)-1),cells(k));
                            end
                            pid = find(P.clst == cnum(j));
                            qid = find(Q.clst == cnum(k));
                            [xc, details] = xcorrtimes(P.times(pid),Q.times(qid));
                            if k == j
                                xc(details.xpts == 0) = 0;
                            end
                            DATA.xcorrval.times = details.xpts;
                            DATA.xcorrs(nxc).xc = xc;
                            DATA.xcorrs(nxc).n = [length(pid) length(qid)];
                            DATA.xcorrs(nxc).valid = 1;
                            if bycell
                                DATA.xcorrs(nxc).cells = [DATA.CellList(e,cid(j),cnum(j)-1) DATA.CellList(e,cid(k),cnum(k)-1)];
                                DATA.xcorrs(nxc).separation = cid(j)-cid(k);
                                DATA.xcorrs(nxc).probes = [cid(j) cid(k)];
                                DATA.xcorrs(nxc).clnum = [cnum(j)-1 cnum(k)-1];
                                if diff(DATA.xcorrs(nxc).cells) > 0
                                    fprintf('Order error\n');
                                end
                            else
                                DATA.xcorrs(nxc).cells = [DATA.CellList(e,j,1) DATA.CellList(e,k,1)];
                                DATA.xcorrs(nxc).separation = j-k;
                                DATA.xcorrs(nxc).probes = [j k];
                                DATA.xcorrs(nxc).clnum = [1 1];
                            end
                            DATA.xcorrs(nxc).eid = e;
                            DATA.xcorrs(nxc).shapexc = ShapeCorr(P,Q);
                            DATA.xcorrs(nxc).calctime = now;
                            DATA.xcorrs(nxc).effic = max(details.efficacy)
                        end
                    end
                end
            end
            toc
            SaveExtras(DATA);
            set(DATA.toplevel,'UserData',DATA);
            PlotXcorrs(DATA, DATA.xcorrs, expts, bycell);
        elseif strcmp(fcn,'usealltrialscl2')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),2, 'reset');
        elseif fcn == 8
            DATA.plot.gmcid = ~DATA.plot.gmcid;
            set(a,'Checked',onoff{DATA.plot.gmcid+1});
        elseif fcn == 9 %'xcorrselected
            C= getappdata(DATA.toplevel,'Clusters');
            Clusters = C{DATA.currentpoint(1)};
            e = DATA.currentpoint(1);
            SetFigure(DATA,DATA.tag.xcorr,'front');
            ps = find(DATA.selectprobe(e,:));
            np = length(ps);
            if np == 1
                np = length(Clusters{ps}.next)+1;
                oneprobe = 1;
            else
                oneprobe = 0;
            end
            for j = 1:np
                for k = 1:j
                    if oneprobe
                        P = Clusters{ps};
                        probe = ps;
                        probeb = k;
                        pid = find(P.clst == 1+j);
                        qid = find(P.clst == 1+k);
                        xc = xcorrtimes(P.times(pid),P.times(qid));
                    else
                        P = Clusters{ps(j)};
                        probe = ps(j);
                        probeb = ps(k);
                        pid = find(P.clst == 2);
                        Q = Clusters{ps(k)};
                        qid = find(Q.clst == 2);
                        xc = xcorrtimes(P.times(pid),Q.times(qid));
                    end
                    if k == j
                        xc(201) = 0;
                    end
                    mysubplot(np,np,k+(j-1)*np);
                    plot(-200:200,xc);
                    axis('tight');
                    if k == j
                        title(sprintf('P%d',probeb));
                    end
                    if k == 1
                        ylabel(sprintf('P%d',probe));
                    end
                end
            end
        elseif fcn == 10 %Xcorr for all adjacent probes
            C= getappdata(DATA.toplevel,'Clusters');
            Clusters = C{DATA.currentpoint(1)};
            SetFigure(DATA,DATA.tag.allxy,'front');
            [nr,nc] = Nsubplots(length(Clusters));
            SetFigure(DATA,DATA.tag.allxy);
            tic;
            for j = 1:length(Clusters)-1
                mysubplot(nr,nc,j);
                xr = get(gca,'xlim');
                yl = get(gca,'ylim');
                hold off;
                P = Clusters{j};
                Q = Clusters{j+1};
                pid = find(P.clst == 2);
                qid = find(Q.clst == 2);
                xc = xcorrtimes(P.times(pid),Q.times(qid),'method',5);
                xpts = linspace(xr(1),xr(2),length(xc));
                plot(xpts,yl(1)+(xc.*diff(yl)./max(xc)),'k','linewidth',2);
                drawnow;
            end
            toc;
        elseif fcn == 11
            DATA =  PlotAllProbeXY(DATA,'oneprobe');
        elseif fcn == 12
            Get2DMaxima(DATA);
        elseif fcn == 13
            DATA = FindMissing(DATA);
%            CellFinder(DATA,1);
        elseif strmatch(fcn,{'showexcluded' 'showcellmeans'})
            DATA.plotspk.(fcn) = ~DATA.plotspk.(fcn);
            set(a,'Checked', onoff{1+DATA.plotspk.(fcn)});
        end
        set(DATA.toplevel,'UserData',DATA);
       
function ReplotXcorrs(a,b, type)
       
    DATA = GetDataFromFig(a);
    xcorrs = getappdata(DATA.toplevel,'xcorrs');
    if sum(strcmpi(type, {'meanim' 'xcorr' 'syncspikes' 'histograms'}))
        DATA.plot.xcorrtype = type;
        SetMenuCheck(a,'exclusive');
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmpi(type, 'Shape/Efficacy')
        DATA.plot.xcorrpoptype = type;
        GetFigure(DATA.tag.xcorrpop,DATA.toplevel);
        ClearPlot;
        hold off;
        for j = 1:length(xcorrs)
            plot(xcorrs(j).shapexc, max(xcorrs(j).efficacy),'o',...
                'buttondownfcn',{@PlotXcorr, xcorrs(j).p(1),xcorrs(j).p(2)});
            hold on;
        end
        set(gca,'yscale','log','ylim',[min([xcorrs.shapexc]) 1]);
        ylabel('Efficacy');
        xlabel('Shape Xcorr');
    end

    

function  SetCheckExclusive(a)
%Set current menu item checked, all others off
    m = get(a,'parent');
    c = get(m,'Children');
    set(c,'Checked','off');
    set(a,'checked','on');

function xc = xShapeCorr(P,Q)
    xc = corrcoef(P.MeanSpike.ms(:),Q.MeanSpike.ms(:));
    xc = xc(1,2);

    
function PlotExpts(DATA)
    SetFigure(DATA, DATA.tag.expt);
    Expts = Expts{DATA.currentpoint(1)};
   ex =      PlotExpt(Expts{1},'rcnmin',10);
   ts = now;
    for j = 1:length(Expts)
        PlotRates(Expts{j},ex);
%        PlotExpt(Expts{j},'hold');
        hold on;
    end
     mytoc(ts);
     
function FindDuplicates(DATA, cell)
    eid = cell(1);
    p = cell(2);
    cl = cell(3);
    C= getappdata(DATA.toplevel,'Clusters');
    Clusters = C{eid};
    nc = 1;
    
    D = getappdata(DATA.toplevel,'DuplicatePlot');
    if isfield(D,'epc') && sum(D.epc == cell) == 3
        cells = getappdata(DATA.toplevel,'xcCellList');
        p =  D.p;
%        cid = find([cells.p] == p & [cells.cl] == cl);
        cellid = DATA.CellList(eid, cell(2),cl);
        xcorrs = getappdata(DATA.toplevel,'xcorrs');
        cid = find([cells.p] == cell(2) & [cells.cl] == cl);
    else
    for k = 1:length(Clusters)
        if DATA.plot.dprimemin == 0 || (Clusters{k}.fitdprime(1) < DATA.plot.dprimemin)
            cells(nc).p = k;
            cells(nc).cl = 1;
            cells(nc).dropi = Clusters{k}.dropi(3);
            cells(nc).dprime = Clusters{k}.fitdprime(1);
            if isfield(Clusters{k},'clst')
                t = find(Clusters{k}.clst == 2);
                Clusters{k}.times = Clusters{k}.times(t);
            end
            nc = nc+1;
        end
        for j = 1:length(Clusters{k}.next)
            if isfield(Clusters{k}.next{j},'times') && ...
                    (DATA.plot.dprimemin == 0 || Clusters{k}.next{j}.fitdprime(1) < DATA.plot.dprimemin)
                cells(nc).p = k;
                cells(nc).cl = 1+j;
                cells(nc).dropi = Clusters{k}.next{j}.dropi(3);
                cells(nc).dprime = Clusters{k}.next{j}.fitdprime(1);
                nc = nc+1;
            end
        end
    end
    cid = find([cells.p] == p & [cells.cl] == cl);
    cellid = DATA.CellList(eid, p,cl);
    SetFigure(DATA,DATA.tag.xcorr);
    [cells, xcorrs] = PlotAllXcorr(DATA, Clusters,cells,'sublist',cid,'callback',@PlotXcorr);
    id = find(max(cat(1,xcorrs.efficacy)') > 0.2);
    p = unique([xcorrs(id).p]);
    end
    otherp = setdiff(p,cid);
    D.p = p;
    D.eid = eid;
    D.cell = cellid;
    D.epc = cell;
    setappdata(DATA.toplevel,'DuplicatePlot',D);
    axdata.eid = eid;
    axdata.p = p;
    axdata.cell = cellid;

    GetFigure(DATA.tag.xcorrpop,DATA.toplevel);
    ClearPlot;
    subplot(2,1,1);
    hold off;
    for j = 1:length(xcorrs)
        if sum(xcorrs(j).efficacy) > 0
            ydata(j) = max(xcorrs(j).efficacy);
        h = plot(xcorrs(j).shapexc, ydata(j),'o',...
            'buttondownfcn',{@PlotXcorr, xcorrs(j).p(1),xcorrs(j).p(2)});
        [isdup, dupi] = ismember(xcorrs(j).p,otherp);
        if sum(isdup)
            a = find(p == otherp(max(dupi)));
            set(h,'color',DATA.colors{a},...
                'markerfacecolor',DATA.colors{a});
        end
        hold on;
        else
            ydata(j) = NaN;
        end
    end
    set(gca,'yscale','log','ylim',[min(ydata) 1]);
    ylabel('Efficacy');
    xlabel('Shape Xcorr');
    subplot(2,1,2);
    dx = diff(minmax([cells(p).dprime]))./50;
    for j = 1:length(p)
        cstr = '';
        [a,b] = isacell(DATA,eid,cells(p(j)).p,cells(p(j)).cl);
        if a
            cstr = sprintf('C%d',b);
        elseif b < 0
            cstr = sprintf('Dup%d',-b);            
        end
        h = plot(cells(p(j)).dropi,cells(p(j)).dprime,'o','color',DATA.colors{j});
        text(cells(p(j)).dropi+dx,cells(p(j)).dprime,...
            sprintf(' %d/%d%s',cells(p(j)).p,cells(p(j)).cl,cstr),...
            'horizontalalignment','left','color',DATA.colors{j});
        AddLineContextMenu(DATA, h, eid,cells(p(j)).p,'cellnum',p(j),'duplicate',cellid);
        set(h,'buttondownfcn',{@PlotXcorr, cid,p(j)},...
                'markerfacecolor',DATA.colors{j});
        hold on;
    end
    ylabel('Fit Dprime');
    xlabel('DropIndex');
    set(gca,'UserData',axdata);
    title(sprintf('%d/%d',cells(cid).p,cells(cid).cl));
    PlotXcorr(DATA, 'histograms',p);
        
function PlotXcorrs(DATA, xcorrs, expts, bycell)
    ClearPlot;
    if bycell
        xcorrs = xcorrs([xcorrs.valid] == 1);
        cellids = cat(1,xcorrs.cells);
        cellids(isnan(cellids)) = 0;
    else
        cellids = cat(1,xcorrs.probes);
    end
    exids = cat(1,xcorrs.eid);
    weights = prod(cat(1,xcorrs.n)');
    cells = unique(cellids);
    cells = setdiff(cells, find(DATA.plot.xcorrexclude));
    cells = cells(cells > 0);
    probes = cat(1,xcorrs.probes);
    np = length(cells);
    ns = 0;
    for j = 1:length(cells)
        ida = find(cellids(:,1) == cells(j) & ismember(exids,expts));
        idb = find(cellids(:,2) == cells(j) & ismember(exids,expts));
        cellpos(j) = (sum(probes(ida,1))+sum(probes(idb,2)))./(length(ida)+length(idb));
        for k = 1:j
            ida = find(cellids(:,1) == cells(j) & cellids(:,2) == cells(k) & ismember(exids,expts));
            idb = find(cellids(:,2) == cells(j) & cellids(:,1) == cells(k) & ismember(exids,expts));
            if length(ida)+length(idb) > 0
                separation(j,k) =( sum([xcorrs(ida).separation]) - sum([xcorrs(idb).separation]))./(length(ida)+length(idb));
            else
                separation(j,k) = cellpos(j)-cellpos(k);
            end
            separation(k,j) = -separation(j,k);
        end
    end
    order = sum(separation < 0);
    [a,b] = sort(order);
    icells = cells; %unsorted
    cells = cells(b);
    for j = 1:length(cells)
        for k = 1:j
            id = find((cellids(:,1) == cells(j) & cellids(:,2) == cells(k)) | (cellids(:,2) == cells(j) & cellids(:,1) == cells(k)));
            id = id(ismember(exids(id),expts));
            if length(id)
                if length(id) > 1
                    xc = WeightedSum(cat(1,xcorrs(id).xc),weights(id));
                else
                    xc = xcorrs(id).xc;
                end
                mysubplot(np,np,k+(j-1)*np,'leftmargin',0.02);
                h = plot(-200:200,xc,'k-','linewidth',2);
                synci = SyncIndices(xc);
                if synci(2) < DATA.crit.synci
                    set(h,'color','r');
                end
                if ~isnan(synci(2))
                    ns = ns+1;
                    syncis(ns,1:length(synci)) = synci;
                    syncis(ns,3) = cells(j);
                    syncis(ns,4) = cells(k);
                end
                axis('tight');
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k)], expts});
                set(h,'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k)],expts});
                if k == j
                    [a,b] = find(DATA.CellList == cells(k));
                    p = 1+mod(b-1,DATA.nprobes);
                    if j == 1 || DATA.plot.xcorrlabeltop == 1
                        title(sprintf('Cell%d',cells(k)));
                    else
                        ii = find(icells == cells(k-1));
                        ij  = find(icells == cells(j));
                        if bycell
                            h = text(xl(1),yl(2),sprintf('C%d at %s (%.1f)',cells(k),ProbeLabel(p, DATA),separation(ii,ij)));
                            set(h,'HorizontalAlignment','left','verticalalignment','bottom');
                        else
                            title(sprintf('P%s',ProbeLabel(k, DATA)));
                        end
                    end
                end
                if k == 1
                    ylabel(sprintf('Cell%d',cells(j)));
                end
            end
        end
    end
    if DATA.plot.synci
    mysubplot(2,2,2);
    myscatter(syncis(:,1),syncis(:,2),'o','ids',syncis(:,3:4));
    end
    
function str = ProbeLabel(p, DATA)
    if ~isempty(DATA) && ~isempty(DATA.ArrayConfig) && length(unique(p)) ==1
        p = unique(p);
        str = sprintf('%d(%d,%d)',p,DATA.ArrayConfig.X(p),DATA.ArrayConfig.Y(p));
    elseif length(unique(p)) == 1
        str = sprintf('%d',mean(p));
    else
        str = sprintf('%.1f',mean(p));
    end
        
    
function args = PlotArgs(DATA)
args = {};
if DATA.exclude.offpeak
    args = {args{:} 'exclude' 'offpeak'};
end

if DATA.exclude.onpeak
    args = {args{:} 'exclude' 'onpeak'};
end
if DATA.exclude.onpeak
    args = {args{:} 'exclude' 'cellsonly'};
end

if ~strcmpi(DATA.colorscheme,'Plain')
    args = {args{:} 'colorscheme' DATA.colorscheme};
end
    
function Update(a,b)
DATA = GetDataFromFig(a);

DATA.spoolspikes = GetCheck('SpoolSpikes',DATA.toplevel);
DATA.showspkxy = GetCheck('SpikeXY',DATA.toplevel);
DATA.showspkmean = GetCheck('SpikeMean',DATA.toplevel);
%plotexpt = GetCheck('PlotExpt',DATA.toplevel);
%if plotexpt && strcmp(DATA.plotexpttype,'none')
  %  DATA.plotexpttype = 'means';
%elseif plotexpt == 0
%    DATA.plotexpttype = 'none';
%end
DATA.refitgm = GetCheck('RefitGM',DATA.toplevel);
DATA.plot.showgm = GetCheck('ShowGM',DATA.toplevel);
DATA.plotallxy = GetCheck('PlotAllXY',DATA.toplevel);
DATA.plotxyseq = GetCheck('PlotXYseq',DATA.toplevel);
DATA.plottrighist = GetCheck('PlotTrigHist',DATA.toplevel);
DATA.plotspks = GetCheck('PlotSpks',DATA.toplevel);
DATA.plothist = GetCheck('PlotHistogram',DATA.toplevel);
DATA.exclude.offpeak = GetCheck('PeakMismatch',DATA.toplevel);
DATA.exclude.onpeak = GetCheck('PeakMatch',DATA.toplevel);
DATA.exclude.noncell = GetCheck('NonCell',DATA.toplevel);
DATA.colorscheme = GetPopString('ColorScheme',DATA.toplevel);
DATA.show.exptno = GetCheck('LabelExptno',DATA.toplevel);
DATA.show.exptname = GetCheck('LabelExptName',DATA.toplevel);
DATA.show.ed = GetCheck('LabelEd',DATA.toplevel);
DATA.show.allvpcs = GetCheck('AllVPcs',DATA.toplevel);

set(DATA.toplevel,'UserData',DATA);

function [str, value, it] = GetPopString(tag, varargin)

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it) 
    value = get(it(1),'value');
    strs = get(it(1),'string');
    str = deblank(strs(value,:));
else
    value = 0;
    str = [];
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



function id = FindSpikes(DATA,C, xcl)
    xid = [];
    e = DATA.currentpoint(1);
    p = C.probe(1);
    t = C.times.*10000;
    Expts = getappdata(DATA.toplevel,'Expts');
    e = floor(C.exptno); %exlusion list is by Exp # not row
    e = C.exptid; %but Expts list is as long as rows
    if max(xcl) > length(Expts{e, 1}.Trials)
        fprintf('Error In Trial Exclusion list\n');
    else
        for j = 1:length(xcl)
            id = find(t > Expts{e, 1}.Trials(xcl(j)).Start(1) & t< Expts{e, 1}.Trials(xcl(j)).End(end));
            xid = [xid id];
        end
    end
    id = setdiff(1:length(C.times),xid);

function FixBoundary(DATA, Clusters, e, p)
    C = Clusters{e}{p};
    xy = xyrotate(C.xy(:,1),C.xy(:,2),C.angle);
    fits{1}.params = C.fitdpparams(1,:);
    fits{2}.params = C.fitdpparams(2,:);
    fits{1}.state.meanlimit = [min(xy(:,1)) C.crit];
    fits{2}.state.meanlimit = [C.crit max(xy(:,1))];
    fits{1}.amp = C.fitdpparams(1,3);
    fits{2}.amp = C.fitdpparams(1,3);
%    B = mydip(fits);
    B = MyDip(xy(:,1),'idlist',C.clst);
    crit = B.x(B.dip(1));
    if DATA.plothist || DATA.refitgm
        oldf = gcf;
        SetFigure(DATA,DATA.tag.hist,'front');
        hst = PlotClusterHistogram(DATA, C, DATA.refitgm,'cluster', DATA.currentcluster);
        plot([crit crit],get(gca,'ylim'));
        id = find(hst.x > min(B.x) & hst.x < max(B.x));
        if ~isempty(id)
        cmax = max(hst.nsp(id));
        plot(B.x,B.y .* cmax./max(B.y),'g-');
        end
        figure(oldf);
    end
    hold on;
    if B.dip(2)./min(B.dip(3:4)) < 0.9 && B.d > 0.1 & B.d < 0.97
        xy = xyrotate([crit crit],[minmax(xy(:,2))'],-C.angle);
        plot(xy(:,1),xy(:,2),'k-');
    end
    
    function plots = PlotClusterPoints(C, uid, cid, varargin)
        plotgmcid = 0;
        finishplot = 0;
        clnum = 1;
        plots = [];
        j = 1;
        colors = mycolors('spkcolors');
        while j <= length(varargin)
            if strncmpi(varargin{j},'colors',6)
                j = j+1;
                colors = varargin{j};
            elseif strncmpi(varargin{j},'finish',5)
                finishplot = 1;
            elseif strncmpi(varargin{j},'plotgmcid',8)
                plotgmcid = 1;
            elseif strncmpi(varargin{j},'xydata',6)
                j = j+1;
                C.xy = varargin{j};
            end
            j = j+1;
        end
        if length(cid) > 1
            xy = XYSpace(C);
            C.xy(:,1) = xy(:,1);
            xy = XYSpace(C.next{1});
            C.xy(:,2) = xy(:,1);
        elseif cid > 1 && length(C.next) >= cid-1 && isfield(C.next{cid-1},'xy')
            C.xy = C.next{cid-1}.xy;
        end
        if isempty(uid)
            uid = 1:size(C.xy,1);
        end
        axdata = get(gca,'UserData');
        plot(C.xy(uid,1),C.xy(uid,2),'.','markersize',1);
        hold on;
        if plotgmcid && isfield(C,'gmfit2d')
            cid = cluster(C.gmfit2d,C.xy);
            id = find(cid == clnum+1);
        elseif isfield(C,'clst')
            id = find(C.clst(uid) == clnum+1);
        elseif C.sign < 0
            id = find(C.xy(uid,1) <  C.crit(1));
        else
            id = find(C.xy(uid,1) >  C.crit(1));
        end
        cells = unique(C.clst(uid));
        for j = 1:length(cells)
            id = find(C.clst(uid) == cells(j));
            if length(id)
                plots(j) = plot(C.xy(uid(id),1),C.xy(uid(id),2),'.','markersize',1,'color',colors{cells(j)});
                if j > 1
                    set(plots(j),'ButtonDownfcn',{@ReplotXY, C.exptid, C.probe(1), cells(j)-1});
                end
            end
        end
        set(gca,'UserData',axdata);
        
function  h = AddTitle(DATA, C, titlemode)
    spstr = [];
    p = C.probe(1);
    if titlemode == 2
        if C.space(1) == 6
            c = DATA.NDSpaceChars(C.space(2));
        else
            c = DATA.NDSpaceChars(C.space(1));
        end
        h=title(sprintf('E%.0fP%s %.2f',C.exptno,ProbeLabel(p,DATA),DistanceMeasure(C, 1,DATA.mahaltype)));
    elseif titlemode == 1
        h= title(sprintf('E%.0fP%s %.2f  %s %.2f',C.exptno,ProbeLabel(p,DATA),C.mahal(1),spstr,C.mahal(4)));
    else
        h= title(sprintf('P%s Ex %.0f Gm %.2f  %s (%.2f,%.2f for space %.0f)%.0f',ProbeLabel(p, DATA),C.exptno,C.mahal(1),spstr,C.mahal(4),C.bestspace(1),C.bestspace(2),C.sign));
    end

function FastAxes(ax)
    set(ax,'xlimmode','manual','ylimmode','manual');
    set(ax,'Xtick',[],'ytick',[]);
    set(ax,'xticklabel',[],'yticklabel',[]);
    set(gca,'zlimmode','manual','climmode','manual');
    set(gca,'xtickmode','manual','ytickmode','manual','ztickmode','manual');
    set(gca,'xticklabelmode','manual','yticklabelmode','manual','zticklabelmode','manual');
    
function plots = PlotClusterXY(DATA, C, varargin)
    titlemode = 0;
    tightplot = 0;
    twoplot = 0;
    clnum = 1;
    axdata.allxy = 0;
    xyargs = {};
    if isfield(C,'rawxy')
        C.xy = C.rawxy;
    end

    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'allxy',5)
            axdata.allxy = 1;
        elseif strncmpi(varargin{j},'plotgmcid',8)
            xyargs = {xyargs{:} varargin{j}};
        elseif strncmpi(varargin{j},'shorttitle',8)
            titlemode = 1;
        elseif strncmpi(varargin{j},'vshorttitle',8)
            titlemode = 2;
        elseif strncmpi(varargin{j},'twoplot',5) %Combine spaces > 1 cluster
            twoplot = 1;
            clnum = [1 2];
        elseif strncmpi(varargin{j},'tight',5)
            tightplot = 1;
        elseif strncmpi(varargin{j},'cellid',6)
            j = j+1;
            clnum = varargin{j};
        elseif strncmpi(varargin{j},'xydata',6)
            j = j+1;
            C.xy = varargin{j};
        else
            xyargs = {xyargs{:} varargin{j}};
        end
        j = j+1;
    end
    lincolor = 'r';
    
    if ~isfield(C,'xy') %empty probesin grid data
        plots = [];
        return;
    end
    p = C.probe(1);
    e = C.exptid;
    if length(clnum) > 1
        CC = C;
    elseif clnum > 1 && length(C.next) > clnum-2 && isfield(C.next{clnum-1},'xy');
        CC = C.next{clnum-1};
    else
        CC = C;
    end
    axdata.toplevel = DATA.toplevel;
    axdata.eid = e;
    axdata.probe = p;
    xcl = FindExcludedTrials(DATA, e, p,1,C);
    if length(xcl)
        uid = FindSpikes(DATA, C, xcl);
    else
        uid = 1:size(C.xy,1);
    end
    if DATA.plot.density
        lincolor = 'w';
        [plots, D] = DensityPlot(CC.xy(uid,1),CC.xy(uid,2),'ynormal');
        if DATA.plot.scaledensity
            r = CalcRadius(CC,[D.x(:) D.y(:)]);
            id = find(r < 1);
            cmax = max(D.z(id));
            caxis([0 cmax]);
        end
        hold on;
    else
        plots = PlotClusterPoints(C,uid,clnum, xyargs{:});
        if DATA.show.linecontextmenus
        for j = 1:length(plots)
            AddLineContextMenu(DATA, plots, C.exptid, C.probe(1));
        end
        end
        set(gca,'UserData',axdata,'uicontextmenu',AddContextMenu(gca,'subplot'));
    end
    if length(clnum) > 1
        spstr = 'D1/D2';
    elseif clnum > 1 && length(C.next) >= clnum-1 && isfield(C.next{clnum-1},'space')
        spstr = [DATA.SpaceTypeLabels{C.next{clnum-1}.space(1)} sprintf(' %d',C.next{clnum-1}.space(2:end))];
    elseif isfield(C,'space')
        spstr = [DATA.SpaceTypeLabels{C.space(1)} sprintf(' %d',C.space(2:end))];
    else
        spstr = '';
    end
    if isfield(DATA,'GaussFitdp') && size(DATA.GaussFitdp,1) >= e
%        spstr = [sprintf('(Fit %.1f)',DATA.GaussFitdp(e,p,2)) spstr];
    end
    if isfield(C,'gmfit2d') && DATA.plot.showboundary
            xy = GMBoundary(C); %,'plot to show image
            plot(xy(:,1),xy(:,2),'k');
    end
    axis('tight');
    gmd = DistanceMeasure(C, 1,DATA.mahaltype);
    if isfield(C,'bestspace') && C.space(1) == 6
        if titlemode == 2
            c = DATA.NDSpaceChars(C.space(2));
            h=title(sprintf('E%.0fP%d %.2f ND%c',C.exptno,p,DistanceMeasure(C, 1, DATA.mahaltype),c));
        elseif titlemode == 1
           h= title(sprintf('E%.0fP%d %.2f  %s %.2f',C.exptno,p,gmd,spstr,C.mahal(4)));
        else
           h= title(sprintf('P%d Ex %.0f Gm %.2f  %s (%.2f,%.2f for space %.0f)%.0f',p,C.exptno,gmd,spstr,C.mahal(4),C.bestspace(1),C.bestspace(2),C.sign));
        end
    else
        if titlemode == 2
            c = DATA.SpaceChars(C.space(1));
            h=title(sprintf('E%.0fP%s %.2f%c',C.exptno,ProbeLabel(p, DATA),gmd,c));
        elseif titlemode == 1
            h=title(sprintf('E%.0fP%s %.2f %.2f %.1f',C.exptno,ProbeLabel(p,DATA),gmd,C.dropi(3)));
        else
            h =title(sprintf('P%s Ex %.0f Gm %.2f  %s (%.2f)%.0f Drop%.1f',ProbeLabel(p,DATA),C.exptno,gmd,spstr,C.mahal(4),C.sign,C.dropi(3)));
        end
    end
    if C.shape == 1
        DrawEllipse(C,'color',lincolor);
    elseif C.shape == 2
        DrawEllipse(C,'color',lincolor);
    elseif C.shape == 0
        C.down = 0;
        DrawEllipse(C,'color',lincolor);
    end
    for j = 1:length(C.next) 
        if isfield(C.next{j},'shape')
        C.next{j}.down = 0;
        if C.next{j}.shape == 0
            DrawEllipse(C.next{j},'color',DATA.colors{j+2});
        end
        end
    end
    setappdata(0,'control_is_down',0);
    if isacell(DATA,e,p)
        sz = get(h,'fontsize');
        sz = sz .* 1.4;
        set(h,'color','b','fontweight','bold','fontsize',sz);
    elseif C.auto == 1
        set(h,'color','r');
    elseif DATA.plot.density
        set(h,'color','w');
    end
    
    if DATA.plottrighist
        AddTrigHist(DATA,C, clnum);
    end
    
    if tightplot
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
    end

    function [r,y] = CalcRadius(E,xy)
                
        rx = E.xyr(3);
        ry = E.xyr(4);
        if isfield(E,'aspectratio') & E.aspectratio > 0
            xys = xyrotate(xy(:,1)-E.xyr(1),(xy(:,2)-E.xyr(2)) ./E.aspectratio,E.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./E.aspectratio)).^2;
        else
            xys = xyrotate(xy(:,1)-E.xyr(1),xy(:,2)-E.xyr(2),E.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
        end
        if nargout > 1
            y = xys(:,2);
        end

function x = Rprime(r)
    x = sqrt(r);
    
function xy = XYSpace(C)
    if C.shape == 0
        [x, y] = CalcRadius(C, C.xy);
        xy(:,1) = Rprime(x);
        xy(:,2) = y;
    else
        xy = xyrotate(xy(:,1),xy(:,2),E.angle);
    end

            
function [x, nsp] = PlotHist(xy, varargin)
    E = [];
    gmfit = [];
    noplot = 0;
    j = 1;
    while j <= length(varargin)
        if isfield(varargin{j},'pos') || isfield(varargin{j},'crit')
            E = varargin{j};
            p = E.probe(1);
            eid = E.exptid;
        elseif isfield(varargin{j},'gmfit1d') || isfield(varargin{j},'gmfit2d')
            gmfit = varargin{j};
        elseif iscell(varargin{j}) && isfield(varargin{j}{1},'gmfit2d')
            gmfit = varargin{j}{eid,p};
        elseif strcmp(varargin{j},'noplot')
            noplot = 1;
        end
        j = j+1;
    end
    

    if E.shape == 0
        r = CalcRadius(E, E.xy);
        [nsp,x] = hist(sqrt(r),200);
        scale = trapz(x,nsp);
    else
        xy = xyrotate(xy(:,1),xy(:,2),E.angle);
        [nsp,x] = hist(xy(:,1),200);
        scale = trapz(x,nsp);
    end
    if noplot
        return;
    end
    hold off;
        bar(x,nsp,1);

        
    if E.shape == 0
        hold on;
        plot([1 1],[0 max(nsp)],'r');
    elseif isfield(gmfit,'gmfit1d')
        [a,b] = GMDip(xy,gmfit.gmfit1d);
        hold on;
        plot(b.gxy(:,1),b.gxy(:,2).*scale,'color','g');
        plot(b.gxy(:,1),b.gxy(:,3).*scale,'color','g');
        plot(b.gxy(:,1),sum(b.gxy(:,[2 3]),2).*scale,'r');
        plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
        plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
    else
        hold on;
        b.gxy(:,1) = linspace(min(xy(:,1)),max(xy(:,1)));
        gmsd = NaN;
        plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
    end
    
    
function Expt = CountExptSpikes(DATA,Expt,C,clnum)
%Expt = CountExptSpikes(DATA,Expt,C,clnum)
%finds spikes in C that match trials in Expt, and
%adds them to the Expt.Trials strucutre.
%clnum 1 = cluster 1, i.e. C.clst == 2
    if isempty(Expt) || isempty(C)
        Expt.Header.cellnumber = 0;
        Expt.Header.probe = 0;
        return;
    end
    cellid = DATA.CellList(C.exptid,C.probe(1),clnum);
    if cellid > 0
        xcl = FindExcludedTrials(DATA,C.exptid,C.probe(1),clnum,C);
        xid = [Expt.Trials(xcl).id];
        Expt.Header.cellnumber = cellid;
    else
        xid = [];
        Expt.Header.cellnumber = 0;
    end
    Expt.Header.probe = C.probe(1);
    Expt = CountSpikes(Expt, C, clnum, xid);
    return;
    spkt = C.times(C.clst == clnum) .*10000;

     for j = 1:length(Expt.Trials)
         id = find(spkt > Expt.Trials(j).Start(1) & spkt < Expt.Trials(j).End(end));
         Expt.Trials(j).Spikes = round(spkt(id)'-Expt.Trials(j).Start(1));
         Expt.Trials(j).count = sum(spkt(id) > 500 & spkt(id) < Expt.Trials(j).dur+500);
     end

    
 function Expt = PlotExptCounts(DATA, e, p, cl, varargin)
     Expt = [];
     if strncmpi(DATA.plotexpttype,'none',4)  
         return;
     end
     tag = DATA.tag.expt; %set to empty to not set figure
       plotargs = {};

       j = 1;
       while j <= length(varargin)
           if strcmpi(varargin{j},'hold')
               plotargs = {plotargs{:} 'hold'};
           elseif strcmpi(varargin{j},'tag')
               j = j+1;
               tag = varargin{j};
           end
           j = j+1;
       end
       
       
       if ~isempty(tag)
           SetFigure(DATA,tag);
       end
     [Clusters, DATA] = CheckClusterLoaded(DATA, e);
     Expts = getappdata(DATA.toplevel,'Expts');
     clnum = cl;
     Expt = CountExptSpikes(DATA, Expts{e,1},Clusters{e}{p},clnum);
       if strncmp(Expt.Header.expname,'image.orXob',11)
           plotargs = {plotargs{:} 'reverse'};
       end
     if strncmpi(DATA.plotexpttype,'trialcounts',10)
         PlotExpt(Expt,'seqt','rcnmin',10,plotargs{:});
     elseif strncmpi(DATA.plotexpttype,'none',4)
     else
         PlotExpt(Expt,'shown','rcnmin',10,'fbox',plotargs{:});
     end


function  C = GetCurrentCluster(DATA)
    Clusters = getappdata(DATA.toplevel,'Clusters')
    C = Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)};

     
function  DATA = ConditionalPlotXY(DATA, C, force)     
    
    if isempty(C)
        C = GetCurrentCluster(DATA);
    end
    if DATA.showspkxy || force
    SetFigure(DATA,DATA.tag.xyplot,'front');
    hold off;
    PlotClusterXY(DATA,C,'cellid',DATA.currentcluster);
    AddCellLabels(DATA, DATA.currentpoint(1), DATA.currentpoint(2),'NW');
%If this is a new Expt/Probe, clear newcut
    if DATA.currentpoint(1) ~= DATA.NewCut.exptid || DATA.currentpoint(2) ~= DATA.NewCut.probe
        DATA.NewCut.exptid = 0;
        DATA.NewCut.probe = 0;
    end
    if DATA.plottrighist
        AddTrigHist(DATA,C,DATA.currentcluster);
    end
    end

  function DATA  = ShowData(DATA, ex,p, varargin)
  
  showall = 1;
  j = 1;
  while j <= length(varargin)
      if strncmpi(varargin{j},'oneprobe',7)
          showall = 0;
      end
      j = j+1;
  end
      
[Clusters, DATA] = CheckClusterLoaded(DATA, ex);
Expts = getappdata(DATA.toplevel,'Expts');
    oldf = gcf;
Expt = [];
if size(Expts,1) >= ex
    if size(Expts,2) >= p
        Expt = Expts{ex,p};
    else
        Expt = Expts{ex,1};
    end
    DATA.Expt  = Expt;
end
if showall == 0
    it = findobj(DATA.toplevel,'Tag','ProbeList');
    set(it,'value',p)
    it = findobj(DATA.toplevel,'Tag','ExptList');
    set(it,'value',ex)
end

if DATA.datatype == 2
C = Clusters{ex}{p}.cluster{DATA.templatesrc};
elseif DATA.useautoclusters
    AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
    C = AutoClusters{ex}{p};
else
C = Clusters{ex}{p};
end
if ~C.probe(1)
    C.probe(1)=p;
end
if C.exptid <= 0
    fprintf('!!!Ex %d P %d has 0 exptid\n',ex,p);
    C.exptid = ex;
end
if isfield(C,'Evec')
    eveci = C.Evec.Eval(1)./sum(C.Evec.Eval);
else
    eveci = NaN;
end
if isfield(C.MeanSpike,'muxc')
    muxc = C.MeanSpike.muxc;
else
    muxc = 0;
end
if isfield(C,'bestspace')
    bestspace = C.bestspace(2);
else
    bestspace = NaN;
    C.bestspace = 0;
end
DATA.cellplotmenu = AddContextMenu(DATA,'cellplot');
DATA = MarkCurrentCluster(DATA);
PrintComments(DATA,ex,p);

if isfield(C,'User')
    user = C.user;
else
    user = '??';
end
fprintf('Mahal ND %.2f(%d), 2D %.2f, 1D %.2f. Dropi %.2f (T%.2f). Made %s by %s\n',C.bestspace(1),bestspace,C.mahal(1),C.mahal(4),C.dropi(3),C.Trigger,datestr(C.ctime),user)
spkr = C.ncut./C.nspks;
if isfield(DATA,'GaussFitdp') && size(DATA.GaussFitdp,1) >= ex
    fprintf('Fit %.1f %.1f spkratio %.3f\n',DATA.GaussFitdp(ex,p,1),DATA.GaussFitdp(ex,p,2),spkr);
end
str = [];
if strmatch(DATA.plot.alltype,{'followcorr' 'templatescore'})
    if DATA.datatype == 1
        xc = CalcTemplateXcorr(Clusters{DATA.templateid(2)}{DATA.templateid(1)},Clusters{ex}{p});
    else
        xc = CalcTemplateXcorr(DATA.Templates{DATA.templatesrc},Clusters{ex}{p}.cluster{DATA.templatesrc});
    end
        str = sprintf(' xc %.2f',xc);
elseif strmatch(DATA.plot.alltype,{'BuildTimes'})
        str = sprintf(' took %.2f',(DATA.ctimes(ex,p,3)-DATA.ctimes(ex,p,1))*24);
elseif strmatch(DATA.plot.alltype,{'PcGms'})
    [a,b] = max(C.pcgms);
        str = ['GM' sprintf(' %.2f',C.pcgms)];
elseif strmatch(DATA.plot.alltype,{'Tagged'}) & DATA.tagged(ex,p)
    DATA.tagstrings = {'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear'};
    fprintf('Tagged %s\n',DATA.tagstrings{DATA.tagged(ex,p)});
end
if isfield(C,'first') & isfield(C.first,'needmore')
    str = [ str sprintf('Nmore %d',C.first.needmore)];
end
fprintf('Evi %.2f, muxc %.2f %s\n',eveci,muxc,str);


str = sprintf('Ex %.1f: %s',DATA.exptid(DATA.currentpoint(1)),DATA.expnames{DATA.currentpoint(1)});
fprintf('%s\n',str);
SetFigure(DATA, DATA.tag.all);
title(str);

if DATA.plotallxy && showall
    PlotAllProbeXY(DATA);
end
if DATA.plotxyseq
    PlotXYSequence(DATA, [ex p]);
end
if DATA.plotspks
    SetFigure(DATA,DATA.tag.spikes);
    h = QuickSpikes(DATA, [ex p]);
    AddLineContextMenu(DATA, h, ex, p);
    drawnow;
end
if DATA.show.allvpcs && showall == 0
    CallAllVPcs(DATA,ex,p);
end
DATA = ConditionalPlotXY(DATA, C, 0);
if ~strcmp(DATA.plotexpttype,'none')
   DATA.Expt = PlotExptCounts(DATA, ex,p, DATA.currentcluster);
end
if DATA.plothist || DATA.refitgm
    SetFigure(DATA,DATA.tag.hist,'front');
    PlotClusterHistogram(DATA, C, DATA.refitgm,'cluster', DATA.currentcluster);
end
if DATA.refitgm
    if C.shape == 0
        C.r = CalcRadius(C, C.xy);
        [a,b] = GMDip(C.r,0);
    else
    [a,b] = GMDip(C.xy,0);
    end
    if isfield(DATA,'GaussFitdp') && ~isnan(dp) && 0
        DATA.GaussFitdp(ex,p,1) = b.gmdprime;
        DATA.GaussFitdp(ex,p,2) = dp;
        DATA.GaussFitdp(ex,p,3) = now;
        DATA.gmfitpos(ex,p,:)  = c.fitpos;
    end
    plot(b.gxy(:,1),b.gxy(:,2).*scale,'g');
    plot(b.gxy(:,1),b.gxy(:,3).*scale,'g');
    plot(b.gxy(:,1),sum(b.gxy(:,2:3),2).*scale,'r');
    title(sprintf('%d/%d events M%.2f->%.2f Fit %.2f%s',C.ncut,size(C.xy,1),C.mahal(4),b.gmdprime,dp,qstr));
    [a,b,c] = GMfit(C.xy,2+length(C.next),1);
    id = cluster(a,C.xy);
    sgn = (mean(C.xy(id==1))-mean(C.xy(id==2))) .* C.sign;
    if  sgn > 0
        mu =2; su = 1;
    else
        mu=1; su=2;
    end
    SetFigure(DATA,DATA.tag.xyplot,'front');
    hold off;
    if sum(id==mu) < 200
        sz =  5;
    else
        sz = 1;
    end
    plot(C.xy(id==mu,1),C.xy(id==mu,2),'.','markersize',sz);
    hold on;
    if sum(id==su) < 200
        sz =  5;
    else
        sz = 1;
    end
    plot(C.xy(id==su,1),C.xy(id==su,2),'r.','markersize',sz);
    for j = 1:length(C.next)
        plot(C.xy(id==2+j,1),C.xy(id==2+j,2),'.','markersize',sz,'color','g');
    end
    title(sprintf('P%d Ex %.0f Gm %.2f(2D%.2f) (%.2f,%.2f for space %.0f)%.0f',p,C.exptno,C.mahal(4),C.mahal(1),C.mahal(2),C.bestspace(1),C.bestspace(2),C.sign));
end

if DATA.showspkmean
    SetFigure(DATA, DATA.tag.spkmean,'front');
    hold off;
    PlotMeanSpike(C,p,0,'addtitle',str,DATA.plotmeantype, DATA);
end
if DATA.spoolspikes && showall == 1
    xs = '';
    if rem(C.exptno,1) > 0.001
        xs = 'a';
    end
    a = SpoolSpikeFile(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
end


    
if DATA.steptype == 2
    if ishandle(DATA.markcc)
        delete(DATA.markcc);
    end
    SetFigure(DATA, DATA.tag.clusters);
    DATA.markcc = DrawBox(ex, p,1);
end
DATA.currentpoint = [ex, p];
set(DATA.toplevel,'UserData',DATA);
DATA = PlotCellList(DATA,'showfig');
drawnow;
figure(oldf);


function result = PlotClusterHistogram(DATA, C, refit, varargin)
    j = 1;
    cl = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cluster',5)
            j =j+1;
            cl = varargin{j};
        end
        j = j+1;
    end
    ex = C.exptid;
    p = C.probe;
    if cl > 1 && cl <= length(C.next)+1 && isfield(C.next{cl-1},'xy')
        X = C;
        C = C.next{cl-1};
        C.clst = X.clst;
        C.exptid = X.exptid;
    end
    if DATA.plot.showgm
    GMfits = getappdata(DATA.toplevel,'GMfits');
    [x, nsp] = PlotHist(C.xy,C,GMfits);
    else
        [x, nsp] = PlotHist(C.xy,C);
    end
    result.x = x;
    result.nsp = nsp;
    hold on;
    scale = trapz(x,nsp);
    if isfield(C,'fitdpparams') && refit ==0
        fits{1}.params = C.fitdpparams(1,:);
        fits{2}.params = C.fitdpparams(2,:);
        dp = abs(C.fitdprime(1));
        c.fitpos = C.fitdprime(2:3);
        c.diff = abs(fits{1}.params(1)-fits{2}.params(1));
        c.sd = sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
    else
        [dp, fits,c] = Fit2Gauss(C,200);
    end
    if sum(c.fitpos) == 2
        qstr = [];
    else
        qstr = '*';
    end
    if C.shape == 0
        b.gxy(:,1) = linspace(min(x),max(x),100);
    elseif DATA.plot.showgm && isfield(GMfits{ex}{p},'gmfit1d') 
        [a,b] = GMDip(C.xy,GMfits{ex}{p}.gmfit1d);
        plot(b.gxy(:,1),b.gxy(:,2).*scale,'g');
        plot(b.gxy(:,1),b.gxy(:,3).*scale,'g');
        plot(b.gxy(:,1),sum(b.gxy(:,[2 3]),2).*scale,'r');
        gmsd = sqrt(mean(b.G.Sigma));
        fprintf('Diffs %.2f,%.2f, Sds %.2f %.2f\n',abs(diff(b.G.mu)),c.diff,gmsd,c.sd)
        plot([C.crit(1) C.crit(1)],[0 max(nsp)],'r');
    else
        b.gxy(:,1) = linspace(min(C.xy(:,1)),max(C.xy(:,1)));
        gmsd = NaN;
        plot([C.crit(1) C.crit(1)],[0 max(nsp)],'r');
    end

    
    if ~isnan(dp)
%     fx = linspace(min(C.xy(:,1)),max(C.xy(:,1)),200);
    fya = FitGauss(b.gxy(:,1), fits{1}.params, 'eval');
    fyb = FitGauss(b.gxy(:,1), fits{2}.params, 'eval');
    yl = get(gca,'ylim');
    plot(b.gxy(:,1),fya,'m-','linewidth',2);
    plot(b.gxy(:,1),fyb,'c-','linewidth',2);
    set(gca,'ylim',yl);
    end
    title(sprintf('%d/%d events M%.2f From Fit %.2f%s',C.ncut,size(C.xy,1),C.mahal(4),dp,qstr));


function AddTrigHist(DATA, C, cl, varargin)
    %superimposes a histogram of trigger values on current axes.
if ~isfield(C,'vhist');
    return;
end
j = 1;
while j <= length(varargin)
    j = j+1;
end

drawall = 1;
cl = 1;  %%Forces drawing of all
if cl > 1 && length(C.next) > cl-2 && isfield(C.next{cl-1},'vhist')
    C = C.next{cl-1};
end
    
    h = ishold;
hold on;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
x = linspace(xl(1),xl(2),length(C.vhist));
plot(x,yl(1)+C.vhist.*diff(yl)./max(C.vhist),'color',DATA.colors{cl+1});
if drawall
    for j = 1:length(C.next)
        if isfield(C.next{j},'vhist')
            x = linspace(xl(1),xl(2),length(C.next{j}.vhist));
            plot(x,yl(1)+C.next{j}.vhist.*diff(yl)./max(C.next{j}.vhist),'color',DATA.colors{j+2});
        end
    end
end
if ~ishold
    hold off;
end

function PlotCorrelogram(C, varargin)
    ta = C.times(C.clst ==1);
    tb = C.times(C.clst ==2);
    xc(:,1) = xcorrtimes(ta,tb); 
    xc(:,2) = xcorrtimes(ta,ta); 
    xc(:,3) = xcorrtimes(tb,tb); 
    xc(201,2:3) = 0;
    SetFigure(DATA,DATA.tag.xcorr);
    plot(xc);

function h = PlotMeanSpike(C, p, cluster, varargin)
    addstr = [];
    plots = [ 1 1];
    DATA = [];
    h = [];
    colors = mycolors('spkcolors');
    j = 1; 
    while j <= length(varargin)
        if isstruct(varargin{j}) && isfield(varargin{j},'ArrayConfig')
            DATA = varargin{j};
        elseif strncmpi(varargin{j},'addtitile',5)
            j = j+1;
            addstr = varargin{j};
        elseif strncmpi(varargin{j},'lineonly',5)
            plots = [0 1];
        elseif strncmpi(varargin{j},'imageonly',7)
            plots = [1 0];
        elseif strncmpi(varargin{j},'smallimage',7)
            plots = [4 0];
        elseif strncmpi(varargin{j},'exptim',6)
            plots = [2 0];
        elseif strncmpi(varargin{j},'meanlines',6)
            plots = [3 0];
        elseif strncmpi(varargin{j},'twoim',5)
            plots = [1 2];
        end
        j = j+1;
    end
    
    if isfield(DATA,'probes')
        nprobes = DATA.nprobes;
    else
        nprobes = size(C.MeanSpike.ms,1);
    end

    if size(C.MeanSpike.ms,1) == 1
        chspk = 1;
        voff = 0;
    else
        chspk = C.probe(1) + [-1:1];
        chspk = chspk (chspk > 0 & chspk <= nprobes);
        if isfield(DATA,'voffset') && length(DATA.voffset) > max(chspk)
            voff = DATA.voffset(chspk)-DATA.voffset(p);
        else
        voff = [-1:1] .*2;
        end
    end
    
    if cluster == 0 %plot all
        nclusters = 1;
        for j = 1:length(C.next)
            if isfield(C.next{j},'MeanSpike')
                nclusters = nclusters+1;
            end
        end
        if size(C.MeanSpike.ms,1) == 1
            PlotMeanSpike(C,1,1,'meanlines');
            title(sprintf('P%s/%d Ex %.1f Gm %.2f (%.2f) %s',ProbeLabel(p, DATA),cluster,C.exptno,C.mahal(1),C.mahal(2),addstr));

            return;
        elseif plots(1) == 3
            PlotMeanSpike(C,p,1,'meanlines',DATA);
        return;
        end
        nc = 1;
        if sum(plots > 0) > 1
            subplot(nclusters,2,1);
            np = 2;
        else
            subplot(nclusters,1,1);
            np = 1;
        end
        if plots(1) == 1
            PlotMeanSpike(C, p, 1, 'imageonly');
        end
        if np == 2
        subplot(nclusters,2,2);
        end
        if plots(2) == 2
            PlotMeanSpike(C, p, -1, 'imageonly');
        elseif plots(2) > 0
            PlotMeanSpike(C, p, 1, 'lineonly');
        end
        nc = 1;
        for j = 1:length(C.next)
            if isfield(C.next{j},'MeanSpike')
                subplot(nclusters,np,1+nc*np);
                if plots(1) == 1
                PlotMeanSpike(C, p, j+1, 'imageonly');
                end
                if np == 2
                subplot(nclusters,np,2+nc*np);
                end
                if plots(2)
                PlotMeanSpike(C, p, j+1, 'lineonly');
                end
                nc = nc+1;
            end
        end
        return;
    end
    if sum(plots > 0) > 1
        subplot(1,2,1);
    end
    if cluster > 1
        if length(C.next) > cluster-2
            C.next{cluster-1}.exptno = C.exptno;
            C = C.next{cluster-1};
        end
    end
    if ~isfield(C,'MeanSpike') %can happen with new cluster
        return;
    end
    if plots(1) == 1 || plots(1) == 4
        hold off;
        if cluster < 0
            h(1) = imagesc(C.MeanSpike.mu);
        elseif plots(1) == 4 && isfield(C,'chspk')
            h(1) = imagesc(C.MeanSpike.ms(C.chspk,:));
        else
            h(1) = imagesc(C.MeanSpike.ms);
        end
        if size(C.MeanSpike.ms,1) == 1
            p = 1;
        elseif p <= 0 && isfield(C,'probe');
            p = C.probe(1);
        end
        line([0 5],[p p],'color','r');
        title(sprintf('P%d/%d Ex %.1f Gm %.2f (%.2f) %s',p,cluster,C.exptno,C.mahal(1),C.mahal(2),addstr));
        if sum(plots > 0) > 1
            subplot(1,2,2);
        end
    elseif plots(1) == 3 %mean lines
        if cluster == 0
            subplot(1,1,1);
        hold off;
        voff(1) = 0;
        for j = 1:length(chspk)
            mm(j,:,1) = minmax(C.MeanSpike.ms(chspk(j),:));
            if cluster == 0
            for k = 1:length(C.next)
                mm(j,:,k+1) = minmax(C.next{k}.MeanSpike.ms(chspk(j),:));
            end
            end
            mx(j,1) = min(mm(j,1,:));
            mx(j,2) = max(mm(j,2,:));
            if j > 1
                voff(j) = mx(j-1,2)-mx(j,1);
            end
        end
        voff = cumsum(voff).*0.8;

        
        for j = 1:length(chspk)
             plot(C.MeanSpike.ms(chspk(j),:)+voff(j),'r','linewidth',2);
             hold on;
             plot(C.MeanSpike.mu(chspk(j),:)+voff(j),'color',colors{1},'linewidth',2);
             if cluster == 0
             for k = 1:length(C.next)
                 if isfield(C.next{k},'MeanSpike')
                     plot(C.next{k}.MeanSpike.ms(chspk(j),:)+voff(j),'color',colors{k+2},'linewidth',2);
                 end
             end
             end
        end
        else
            voff = zeros(size(chspk));
            for j = 1:length(chspk)
                plot(C.MeanSpike.ms(chspk(j),:)+voff(j),'color',colors{j+1},'linewidth',2);
                hold on;
            end
            if length(chspk) == 1
                plot(C.MeanSpike.mu(chspk,:)+voff(j),'color',colors{1},'linewidth',2);
            end
        end
    end
    if plots(2)
        hold off;
        v = std(C.MeanSpike.ms');
        id = find(v > max(v)/2);
        
        if length(v) >= p && v(p) < 0.1 %low v range so rescale dp;
            x = max(abs(minmax(C.MeanSpike.ms(:))));
            dpscale = x/3;
        else
            dpscale = 1;
        end
        for j = id
            h(2) = plot(C.MeanSpike.ms(j,:),'r');
            hold on;
            if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j
                plot(C.MeanSpike.dp(j,:).*dpscale,'g');
            end
            plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);
        end
    end

function HitImage(src,b, type)
    ts = now;
DATA = GetDataFromFig(src);
if DATA.profiling
    fprintf('Get took %.3f\n',mytoc(ts));
    ts = now;
    profile on;
end
ax = get(src,'Parent');
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
ex = round(xy(1,2));
p = round(xy(1,1));
if xy(1,1) > p
    setcl = 2;
else
    setcl = 1;
end
[Clusters, DATA] = CheckClusterLoaded(DATA, ex);

bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
zval = NaN;
for j = 1:length(l)
    a = get(l(j));
    if isfield(a,'CData')
        Z = get(l(j),'CData');
        zval = Z(ex,p);
    end
end

    cmenu = get(src,'UIContextMenu');
    axdata = get(gca,'UserData');
offset = DATA.clusteroffset;
cntrl_is_down = getappdata(0,'control_is_down');

if cntrl_is_down
    DATA.selectprobe(ex,p) = ~DATA.selectprobe(ex,p);
    DATA.selecth = DrawBoxes(DATA, 3);
    bt =4;
elseif bt == 2
    c = get(cmenu,'children');
    tid = find( strcmp('Title',get(c,'Tag')));
    [a,b,cl] = isacell(DATA,ex,p);
    if a
        cid = find(cl == setcl+offset);
        if isempty(cid) && setcl+offset > max(cl)
            cid = max(cl); %if hit on 2 but only 1 is defined as a cell
            thiscell = 0;
            if DATA.nclusters(ex,p) >= setcl+offset
                cid = setcl+offset;
            end
        else
            cid = cl(cid);
            thiscell = b(cid);
        end
        if length(b) == 1
            set(c(tid),'Label',sprintf('E%dP%d Cell %d',ex,p,b));
            if length(cid) == 1
                set(c(tid),'foregroundcolor',DATA.colors{cid+1});
            end
        elseif length(cid) == 1 && cid <= length(b)
            othercell = setdiff(b,thiscell);
            set(c(tid),'Label',sprintf('E%dP%d Cell %d (%s)',ex,p,b(cid),sprintf('%d ',othercell)));
            set(c(tid),'foregroundcolor',DATA.colors{cid+1});
        else
            set(c(tid),'Label',sprintf('E%dP%d Cells%s',ex,p,sprintf('%d ',b)));
            if length(cid) == 1
                set(c(tid),'foregroundcolor',DATA.colors{cid+1});
            end
        end
    else
        if setcl+offset > DATA.nclusters(ex,p)
            cid = DATA.nclusters(ex,p);
        else
            cid = setcl+offset;
        end
        if cid <= 0
            cid = 1;
        end
        d = isduplicate(DATA,ex,p,cid);
        if sum(DATA.selectprobe(:)) > 1
            str = [' +' num2str(sum(DATA.selectprobe(:)))];
        else
            str = '';
        end
        if d > 0
            set(c(tid),'Label',sprintf('E%dP%d Dup%d%s',ex,p,d,str));
        else
            set(c(tid),'Label',sprintf('E%dP%d%s',ex,p,str));
        end
    end
    axdata.eid = ex;
    axdata.probe = p;
    axdata.toplevel = DATA.toplevel;
    axdata.clnum = cid;
    set(gca,'userdata',axdata);
    return;
end
if strcmp(type,'CellRates')
    DATA.markexpts = DATA.expnames{ex};
    DATA.currentcell = p;
    fprintf('E%dC%d',ex,p);
    fprintf('Hit %.0f,%.0f %.3f %d\n',ex,p,zval,bt);
    PlotMenu(DATA, [], 'cells', 'rateseqone');
    return;
end
fprintf('Hit %.0f,%.0f %.3f type %d,%d\n',ex,p,zval,type,bt);
%If not cutting clusters, then set currentcluster to match hit
%if ht
if DATA.nclusters(ex,p) >= setcl
    DATA.currentcluster = setcl;
end
if isempty(Clusters) || DATA.show.cellsummary
    DrawBox(ex, p,3,'color','w');
    PlotCellSummary(DATA, ex, p);
    return;
end
if ~isfield(Clusters{ex}{p},'clst')
    cprintf('red','No clst For E%dP%d\n',ex,p);
    return;
end
if ismember(type, [1 3]) % 3 = hit cell image - set cell#
    if DATA.datatype == 2
        C = Clusters{ex}{p}.cluster{DATA.templatesrc};
    else
        C = Clusters{ex}{p};
    end
    nc = Counts(C.clst);
    nstr = sprintf(' %d',nc);
    if isfield(C,'manual')
        xstr  = sprintf('Man:%d',C.manual);
    else
        xstr = [];
    end
    fprintf('P%d E%.0f cut on %s:%s spks %s\n',p,C.exptno,datestr(C.ctime),nstr,xstr);
    if type == 3
       it = findobj('Tag','CellNumberId');
       cells = sum(DATA.CellList(ex,p,:) > 0);
       if cells > 0
           lastcell = DATA.currentcell;
           a = find(DATA.CellList(ex,p,:) > 0);
           id = a(a > offset & a <= offset+2);
           if setcl == 2 & length(id) == 2
               id = id(2);
           elseif isempty(id)
               id = offset+1;
           else
               id = id(1);
           end
           DATA.cellcluster = id;         
           DATA.currentcell = DATA.CellList(ex,p,id);
           if DATA.currentcell == 0
               DATA.currentcell = lastcell;
           end
           set(it,'value',DATA.currentcell);
           cit = findobj('Tag','CellCluster');
           set(cit,'value',id+offset);
           cit = findobj('Tag','ClusterModifier');
           set(cit,'value',id+offset);
           if lastcell ~= DATA.currentcell
               CellChanged(DATA)
           end
       end
    end
elseif type == 2
    id = find(DATA.AllPairs(:,1) == ex & DATA.AllPairs(:,2) == p);
    if length(id)
    xc = meanccf(DATA, id, ex, p);
    dist(1) = meanmahal(DATA, id, ex);
    dist(2) = meanmahal(DATA, id, p);
    SetFigure(DATA,DATA.tag.onecell);
    plot(xc);
    title(sprintf('Expts %s mahal P%d %.2f, P%d %.2f',num2str(DATA.allexpt(id)),ex,dist(1),p,dist(2)));
    end
    set(gca,'ylim',[0 max(xc)]);
    SetFigure(DATA,DATA.tag.spkmean);
    subplot(1,2,1);
    ms = MeanSpike(DATA, DATA.allexpt(id), p);
    imagesc(ms);
    subplot(1,2,2);
    ms = MeanSpike(DATA, DATA.allexpt(id), ex);
    imagesc(ms);
    C = Clusters{DATA.allexpt(id(1))}{p};
    ex= DATA.allexpt(id(1));
end
if bt == 3
    if DATA.selecth > 0 && ishandle(DATA.selecth)
        delete(DATA.selecth);
    end
    DATA.selectprobe = zeros(size(DATA.selectprobe));
    if p > DATA.currentpoint(2)
        DATA.proberange = DATA.currentpoint(2):p;
    else
        DATA.proberange = p:DATA.currentpoint(2);
    end
    if ex > DATA.currentpoint(1)
        DATA.selectexpts = DATA.currentpoint(1):ex;
    else
        DATA.selectexpts = ex:DATA.currentpoint(1);
    end
    DATA.selectprobe(DATA.selectexpts,DATA.proberange) = 1;
    DATA.selecth = DrawBox(DATA.selectexpts,DATA.proberange, 3);
    if strfind(get(gcf,'Tag'),'CellList')
        set(DATA.selecth,'color','w');
    end
elseif bt == 4
    DATA.selectprobe(ex,p) = 1;
elseif bt == 2
    DATA.selectprobe(ex,p) = ~DATA.selectprobe(ex,p);
elseif type == 3
    DATA.selectprobe = zeros(length(Clusters),DATA.nprobes);
    DATA.selectprobe(ex,p) = 1;
end
DATA.currentpoint(1) = ex;
DATA.currentpoint(2) = p;
DATA.currentprobe = DATA.currentpoint(2);
DATA.exptno = DATA.exptid(DATA.currentpoint(1));

DATA = ShowData(DATA,ex,p,'oneprobe');
if ismember(bt, [2 3]) && ~strcmp(DATA.plotexpttype,'none')
    SetFigure(DATA, DATA.tag.expt);
    DATA.Expt = PlotSelectedExpts(DATA);
    set(DATA.toplevel,'UserData',DATA);
end
if DATA.profiling
    fprintf('Took %.2f\n',mytoc(ts));
    profile viewer;
end

function Expt = PlotCombinedExpt(DATA)

[exid, pid] = find(DATA.selectprobe);
Clusters = CheckClusterLoaded(DATA,exid);
Expts = getappdata(DATA.toplevel,'Expts');
clusterid = 1;
SetFigure(DATA,DATA.tag.expt);
args = {'shown'};
    if length(clusterid) == 1
        clusterid = ones(size(exid)).* clusterid;
    end
    h = [];
    for j = 1:length(exid)
        E{j} = CountExptSpikes(DATA, Expts{exid(j),1},Clusters{exid(j)}{pid(j)},clusterid(j));
    end
    Expt = CombineExpts(E);
    PlotExpt(Expt,'fbox',args{:});

function [Expt, AllExpt] = PlotSelectedExpts(DATA, varargin)
%PlotSelectedExpts(DATA, ...)
%PlotSelectedExpts(DATA, ...)
    Expt = [];
    AllExpt = [];
    if strncmpi(DATA.plotexpttype,'none',4)
        DATA.plotexpttype = 'rcmeans';
    end
    args = {};
    exid = [];
    pid = [];
    cellid = 0;
    colors = mycolors;
    clusterid = 1;
    linestyles = {'-' '--' ':' '.-' '-' '--' ':' '.-'};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cell',6)
            j = j+1;
            cellid = varargin{j};
        elseif strncmpi(varargin{j},'cluster',6)
            j = j+1;
            clusterid = varargin{j};
        elseif strncmpi(varargin{j},'colors',6)
            j = j+1;
            colors = varargin{j};
        elseif strncmpi(varargin{j},'expts',6)
            j = j+1;
            exid = varargin{j};
            j = j+1;
            pid = varargin{j};
        else
            args = {args{:} varargin{j}};
        end
        j = j+1;
    end
    if isempty(exid)
    [exid, pid] = find(DATA.selectprobe);
    end
    Clusters = CheckClusterLoaded(DATA,exid);
    Expts = getappdata(DATA.toplevel,'Expts');

    lastf = gcf;

    if length(clusterid) == 1
        clusterid = ones(size(exid)).* clusterid;
    end
    h = [];
    for j = 1:length(exid)
        args = {};
        Expt = CountExptSpikes(DATA, Expts{exid(j),1},Clusters{exid(j)}{pid(j)},clusterid(j));
        if strncmp(Expt.Header.expname,'image.orXob',11)
            args = {args{:} 'reverse'};
        end
        if Expt.Header.rc
            if strncmpi(DATA.plotexpttype,'rcmeans',4)
                args = {args{:} 'condense'};
            else
                args = {args{:} 'rcnmin' 10 'sdfw' 100};
            end
        end
        if j == 1
            res{j} = PlotExpt(Expt,'forcecolor',colors{j},args{:});
        else
            res{j} = PlotExpt(Expt,'forcecolor',colors{j},args{:},'holdon');
        end
        if isfield(res{j},'delaysamples') %an RC expt
            [a,b,c] = PlotRC(res{j},'bestdelay');
            figure(lastf);
            for k = 1:size(b,2)
                eh(k) = plot(a(:,k),b(:,k),'o-','color',colors{j},'linestyle',linestyles{k});
                hold on;
            end
            Expt.handles = eh;
        else
            Expt.handles = res{j}(1).allhandles;
        end
        Expt.handles = Expt.handles(ishandle(Expt.handles));
        set(Expt.handles,'buttondownfcn',{@HitPopPoint, exid(j), pid(j), cellid});
        h(j) = Expt.handles(1);
        labels{j} = sprintf('E%dP%d',exid(j),pid(j));
        AllExpt{j} = Expt;
        hold on;
    end
    if length(h)
        mylegend(h,labels);
    end

    function CellChanged(DATA)
        if DATA.show.watchallcellxy   %hit cell image plot
            PlotAllCellXY(DATA);
        end
        if DATA.show.watchallcellspks
            PlotAllCell(DATA,'allspkswithmean');
        end
        if DATA.show.watchallcellmean
            PlotAllCellMean(DATA,'meanlines');
        end
        it = findobj(DATA.fig.celllist,'Tag','CellNumberId');
        if ~isempty(it)
            cl = caxis;
            nc = size(colormap,1);
            cmap = colormap;
            ci = round(nc * (get(it,'value')-cl(1))./diff(cl));
            if ci <= nc & ci > 0
            set(it,'backgroundcolor',cmap(ci,:),'foregroundcolor',1-cmap(ci,:));
            end
        end

function PlotXcorr(a,b, pa, pb)
DATA = GetDataFromFig(a);
F = GetFigure(a);
cells = getappdata(DATA.toplevel,'xcCellList');
if isempty(b) %from GUI
    plottype = DATA.plot.xcorrtype;
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})
else
    plottype = b;
    bt = 1;
end
if strcmp(pa,'zoom')
    xl = get(gca,'xlim');
    if bt == 2 
        xl(1) = xl(1)*2;
        xl(2) = xl(2)*2;
    elseif bt ==1
        xl(1) = xl(1)/2;
        xl(2) = xl(2)/2;
    end
    set(gca,'xlim',xl);
    return;
end
if bt == 2
    return;
end
GetFigure(DATA.tag.xcorr,DATA.toplevel,'front');
control_is_down = getappdata(0,'control_is_down');
DataClusters = getappdata(DATA.toplevel,'Clusters');
selected = getappdata(DATA.toplevel,'xcSelected');
if control_is_down
    ax = findobj(gcf,'type','axes','userdata',[pa pb]);
    set(ax,'color','c');
    selected(pa,pb) = 1;
elseif nargin > 3
    ax = findobj(gcf,'type','axes','userdata',[pa pb]);
    set(ax,'color','w');
    selected(pa,pb) = 0;
end
setappdata(DATA.toplevel,'xcSelected',selected);
eid = DATA.currentpoint(1);
if strcmp(plottype,'meanim')
    SetFigure(DATA, DATA.tag.xcorr);
    mysubplot(2,4,3);
    P = Cell2Cluster(cells(pa),DataClusters{eid});
    imagesc(P.MeanSpike.ms);
    set(gca','xtick',[],'ytick',[]);
    h = title(sprintf('%d/%d',cells(pa).p,cells(pa).cl));
    set(h,'verticalalignment','middle')
    mysubplot(2,4,4);
    Q = Cell2Cluster(cells(pb),DataClusters{eid});
    cl = minmax(cat(1,Q.MeanSpike.ms(:),P.MeanSpike.ms(:)));
    caxis(cl);
    imagesc(Q.MeanSpike.ms);
    caxis(cl);
    set(gca','xtick',[],'ytick',[]);
    h = title(sprintf('%d/%d',cells(pb).p,cells(pb).cl));
    set(h,'verticalalignment','middle')
    xc = ShapeCorr(P,Q);
    return;
elseif strcmp(plottype,'histograms')
    mysubplot(2,2,2);
    if nargin < 4
        cid = pa;
    elseif control_is_down
        [a,b] = find(selected);
        cid = unique([a b]);
    else
        cid = [pa pb];
    end
    colors = DATA.colors;
    hold off;
 %   ClusterDetails = getappdata(DATA.toplevel,'ClusterDetails');
    for j = 1:length(cid)
        c = cid(j);
        P = Cell2Cluster(cells(cid(j)),DataClusters{eid});
        [x, nsp] = PlotHist(P.xy,P,'noplot');
        h(j) = plot(1:length(x),nsp ./ max(nsp),'color',colors{j},...
            'linewidth',2);
        hold on;
    end
%    legend(h,labels);
ylabel('Cluster');
    if pa > length(cells)/2
        mysubplot(2,2,1);
    else
        mysubplot(2,2,4);
    end
    xl = get(gca,'xlim');
    for j = 1:length(cid)
        c = cid(j);
        P = Cell2Cluster(cells(cid(j)),DataClusters{eid});
        if isfield(P,'triggerV')
            [y,hx] = hist(P.triggerV,100);
            x = xl(1) + (hx-min(hx)).*diff(minmax(xl))./diff(minmax(hx));
            h(j) = plot(1:length(x),y./max(y),'color',colors{j},'linewidth',2);
        else
            h(j) = plot([1 2],[0 0],'color',colors{j});
        end
        hold on;
        labels{j} = sprintf('%d/%d: %.2f %.1f',cells(c).p,cells(c).cl,P.fitdprime(1),P.dropi(3));
    end
    mylegend(h,labels);
    ylabel('Trigger');
    return;
elseif strcmp(DATA.plot.xcorrtype,'syncspikes')
mysubplot(2,2,2);
PlotSyncSpikes(DATA, eid, [cells(pa).p cells(pb).p], [cells(pa).cl cells(pb).cl]);
return;
end
mysubplot(2,2,2);
P = Cell2Cluster(cells(pa),DataClusters{eid});
Q = Cell2Cluster(cells(pb),DataClusters{eid});

[xc, details]  = xcorrtimes(P.times,Q.times);
if pa == pb
    xc(details.midpt) = 0;
end
hold off;
plot(details.xpts,xc,'k','linewidth',2);
axis('tight');
yl = get(gca,'ylim');
xl = get(gca,'xlim');
[a,b] = max(xc);
text(xl(2),yl(2),sprintf('%.0fms %.3f %d/%d %d/%d',...
    details.xpts(b).*1000,max(details.efficacy),...
    cells(pa).p,cells(pa).cl,cells(pb).p,cells(pb).cl),'horizontalalignment','right','verticalalignment','top');
line([0 0],yl,'linestyle','--');
set(gca,'buttondownfcn', {@PlotXcorr, 'zoom', 2});

function Q = Cell2Cluster(cell,Clusters)
 Q = Clusters{cell.p};
triggerV = [];
if isfield(Q,'triggerV')
    triggerV = Q.triggerV(Q.clst == cell.cl+1);
end
if cell.cl > 1
    Q = Q.next{cell.cl-1};
elseif isfield(Q,'clst') && length(Q.clst) == length(Q.times)
    Q.times = Q.times(Q.clst == 2);
end
if ~isempty(triggerV)
    Q.triggerV = triggerV;
end

function HitXcorrAll(a,b, type, cells, expts)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');
exids = cat(1,DATA.xcorrs.eid);
if type == 2
%    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    e = cells(1);
    xc = DATA.xcorrs(cells(6));
    mysubplot(2,2,2);
    xpts = DATA.xcorrval.times;
    plot(xpts,xc.xc,'k');
    if DATA.plot.xcmax < max(xpts) *1000
        set(gca,'xlim',[-DATA.plot.xcmax DATA.plot.xcmax]./1000);
    end
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    [a,b] = max(xc.xc);
    text(xl(2),yl(2),sprintf('Cell %d->%d P%d->%d %.0fms Exp%s',xc.cells(1),...
        xc.cells(2),cells(2),cells(3),xpts(b).*1000,DATA.expnames{e}),...
        'horizontalalignment','right','verticalalignment','top');
    line([0 0],yl,'linestyle','--');
    SetFigure(DATA,DATA.tag.spkmean);
    subplot(1,2,1);
    h = PlotMeanSpike(Clusters{e}{cells(2)},cells(2),cells(4),'imageonly');
    set(h,'ButtonDownFcn',{@HitXYPlot, e , cells(2)},'UserData',DATA.toplevel);
    subplot(1,2,2);
    h = PlotMeanSpike(Clusters{e}{cells(3)},cells(3),cells(5),'imageonly');
    set(h,'ButtonDownFcn',{@HitXYPlot, e , cells(3)},'UserData',DATA.toplevel);

    SetFigure(DATA,DATA.tag.spikes);
    if DATA.plot.synctmax > 0
        PlotSyncSpikes(DATA, e, cells([2 3]), cells([4 5]));
    end
    DATA.selectprobe(e,cells([2 3])) = 1;
    DATA.xcid = cells(6);
    
    if strcmp(DATA.plotexpttype,'means')
        Expts = getappdata(DATA.toplevel,'Expts');
        SetFigure(DATA,DATA.tag.expt);
        subplot(1,2,1);
        PlotExpt(Expts{e,cells(2)});
        subplot(1,2,2);
        PlotExpt(Expts{e,cells(3)});
    end
    set(DATA.toplevel,'UserData',DATA);
%    SpoolAllProbes(DATA, e, AllSpikes(e,cells([2 3])), Clusters{e});
    return;
end
cellids = cat(1,DATA.xcorrs.cells);
probes = cat(1,DATA.xcorrs.probes);
if type == 0 %all probes, not just cells
 ids = find(probes(:,1) == cells(1) & probes(:,2) == cells(2));
 pa = cells(1);
 pb = cells(2);
 np = length(unique(probes));
else
 cells([1 2]) = sort(cells([1 2]),'descend');
 ida = find(cellids(:,1) == cells(1) & cellids(:,2) == cells(2));
 idb = find(cellids(:,2) == cells(1) & cellids(:,1) == cells(2));
 ids = cat(1,ida, idb);
 ids = ids(ismember(exids(ids),expts));
 eid = ismember(exids,expts);
 icells = unique([DATA.xcorrs(eid).cells]);
 icells = icells(icells > 0);
 np = length(icells);
 pa = mean(cat(1,probes(ida,1), probes(idb,2))); 
 pb = mean(cat(1,probes(idb,1), probes(ida,2))); 
end

col = 1;
row = 1;
for j = 1:length(ids)
    X = DATA.xcorrs(ids(j));
    col = col+1;
    if row > floor(np/2)
        if col > np
            if row > floor(np/2)
                row = row+1;
                col =row+1;
            end
        end            
    else
        if col > floor(np/2)
            row = row+1;
            if row == floor(np/2)
                row = row+1;
            end
            col = row+1;
        end
    end
    k = (row-1)*np + col;
    mysubplot(np, np, k,'leftmargin',0.02);
    h = plot(DATA.xcorrval.times,X.xc,'r-');
    si = SyncIndices(X.xc);
    if DATA.plot.synci
        fprintf('E%dP%d,%d synci %s\n',X.eid,X.probes(1),X.probes(2),sprintf('%.2f ',si));
    end
    set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll, 2, [X.eid X.probes X.clnum ids(j)],expts});
    set(h,'buttondownfcn',{@HitXcorrAll, 2, [X.eid X.probes X.clnum ids(j)], expts});
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    text(xl(1),yl(2),sprintf('E%dP%d,%d',X.eid,X.probes(1),X.probes(2)),'horizontalalignment','left','verticalalignment','top');
end
if length(ids)
    while row < np
    col = col+1;
    if row > floor(np/2)
        if col > np
            if row > floor(np/2)
                row = row+1;
                col =row+1;
            end
        end            
    else
        if col > floor(np/2)
            row = row+1;
            if row == floor(np/2)
                row = row+1;
            end
            col = row+1;
        end
    end
    k = (row-1)*np + col;
    mysubplot(np, np, k, 'leftmargin',0.02);
    delete(gca);
    end
end
 mysubplot(2,2,2,'leftmargin',0.02);
 if length(ids) > 1
     weights = prod(cat(1,DATA.xcorrs(ida).n)');
     na = sum(weights);
     if length(ida) > 1
         xca = WeightedSum(cat(1,DATA.xcorrs(ida).xc),weights);
     else
         ida = DATA.xcorrs(ida).xc;
     end
     weights = prod(cat(1,DATA.xcorrs(idb).n)');
     nb = sum(weights);
     if length(idb) == 0
         xcb = [];
         nb = 0;
         xc = xca;
     else
         if length(idb) > 1
             xcb = WeightedSum(cat(1,DATA.xcorrs(idb).xc),weights);
         else
             xcb = DATA.xcorrs(idb).xc;
         end
         xc = WeightedSum(cat(1,xca,fliplr(xcb)),[na nb]);
     end
 else
     xc = DATA.xcorrs(ids).xc;
 end
 xpts = DATA.xcorrval.times;
 plot(xpts,xc,'k-');
 line([0 0],get(gca,'ylim'),'linestyle','--','color','r');
yl = get(gca,'ylim');
xl = get(gca,'xlim');
[a,b] = max(xc);
text(xl(2),yl(2),sprintf('C%d->%d %.0fms (%.0f+-%.1f) P%.1f,%.1f',cells(1),cells(2),xpts(b).*1000,prctile(xc,50),std(xc),pa,pb),'horizontalalignment','right','verticalalignment','top');




 function HitXcorr(a,b, id, ex, cells)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');

if strcmp(id,'zoom')
    xl = get(gca,'xlim');
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    if bt == 2
        xl(1) = xl(1)*2;
        xl(2) = xl(2)*2;
    else
        xl(1) = xl(1)/2;
        xl(2) = xl(2)/2;
    end
    set(gca,'xlim',xl);
    return;
end


a = cells(1);
b = cells(2);
if ex > 0
C = Clusters{ex};
sync = DATA.synci{ex};
xc = sync.ccf{cells(1),cells(2)};
else
   [xc, synci, allxc] = meanccf(DATA, id,a,b);
   si = SyncIndices(xc);
   fprintf('synci %s: %.2f->%.2f\n',sprintf(' %.2f',synci(:,1)),mean(synci(:,1)),si(1));
end
SetFigure(DATA,DATA.tag.onecell);
midpt = ceil(length(xc))./2;
hold off;
if length(id) > 1
    plot(allxc');
    hold on;
    plot(xc,'k','linewidth',2);
else
plot(xc);
end
set(gca,'ylim',[0 max(xc)]);
if ex > 0
    title(sprintf('Ex%d P%d (%.1f),%d (%.1f)',ex,cells(1),C{a}.mahal(1),cells(2),C{b}.mahal(1)));
else
    title(sprintf('%d Expts, P %d,%d',length(id),cells(1),cells(2)));
end
%xc = xcorrtimes(C{cells(1)}.times,C{cells(2)}.times);
%plot(xc);

function synci = SyncIndices(xc)
midpt = ceil(length(xc)./2);
synci(1) = mean([xc(midpt-1) xc(midpt+1)])./xc(midpt);
synci(2) = mean([xc(1:midpt-10) xc(midpt+10:end)])./xc(midpt);

    function NewCellSummary(a,b,e,p)
        DATA= GetDataFromFig(a);
        PlotCellSummary(DATA,e,p);

function PlotCellSummary(DATA, e, p)
    Clusters = LoadCluster(DATA.name, e,'getxy');
    Clusters{p}.exptid = e;
    DATA.exptid = DATA.CellDetails.exptids;
    chspk = [p p-1 p+1];
    chspk = chspk(chspk > 0 & chspk <= DATA.nprobes);
    
    Spks = AddSpikes(DATA, {},  e, p);
    maxv(p,:) = Spks{e,p}.VRange;
    for j = 2:length(chspk)
    Spks = AddSpikes(DATA, Spks,  e, chspk(j));
    maxv(chspk(j),:) = Spks{e,chspk(j)}.VRange;
    end
    DATA.voffset =  cumsum(cat(1,maxv(2:end,2), maxv(end,2))) - cumsum(maxv(:,1));
    DATA.voffset = DATA.voffset .* 0.7;

    F = SetFigure(DATA, 'CellSummary');
    m = findobj(F,'Tag','CellList');
    c = get(m,'Children');
    delete(c);
    cellid = unique(DATA.CellList(e,:,:));
    cellid = cellid(cellid>0);
    for j = 1:length(cellid)
        [b,a] = find(squeeze(DATA.CellList(e,:,:)) == cellid(j));
        uimenu(m,'Label',sprintf('Cell%d P%dC%d',cellid(j),b,a),'callback',{@NewCellSummary, e, b});
    end
    subplot('Position',[0.01 0.01 0.5 0.98]);
    QuickSpikes(DATA, cat(1,Spks{e,chspk}),Clusters{p});
    subplot('Position',[0.52 0.52 0.46 0.46]);
    hold off;
    PlotClusterPoints(Clusters{p},[],1);
    Clusters{p}.down = 0;
    DrawEllipse(Clusters{p},'color','r');

    set(gca,'xtick',[],'ytick',[]);
    AddCellLabels(DATA, e, p);
    subplot('Position',[0.52 0.01 0.46 0.46]);
    PlotMeanSpike(Clusters{p},p,1,'imageonly');

    
    
function [xc, synci, allxc] = meanccf(DATA, id, a,b)
    expts = unique(DATA.allexpt(id));
    if b > a
        c = a;
        a = b;
        b = c;
    end
    for j = 1:length(expts)
        xc(j,:) = DATA.synci{expts(j)}.ccf{a,b};
        synci(j,:) = SyncIndices(xc(j,:));
    end
    allxc = xc;
    xc = mean(xc,1);

function d = meanmahal(DATA, id, a)
Clusters = getappdata(DATA.toplevel,'Clusters');
    expts = unique(DATA.allexpt(id));
    for j = 1:length(expts)
        d(j) = Clusters{expts(j)}{a}.mahal(1);
    end
    d = mean(d);

function ms = MeanSpike(DATA, id, a)
    for j = 1:length(id)
        ms(:,:,j) = Clusters{id(j)}{a}.MeanSpike.ms;
    end
    ms = squeeze(mean(ms,3));

function HitPoint(a,b, C, mode, pt)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');

C = Clusters{DATA.id};
fprintf('Probe %d\n',pt);
sz = 1;
SetFigure(DATA,DATA.tag.onecluster);
hold off; 
plot(C{pt}.xy(:,1),C{pt}.xy(:,2),'.', 'markersize', sz);
hold on;
yl = get(gca,'ylim');
plot([C{pt}.crit(1) C{pt}.crit(1)],yl,'r-');
[a,b] = smhist(C{pt}.xy(:,1),100);
plot(b,yl(1)+ a.*range(yl)./max(a),'r');
title(sprintf('Probe %d',pt));


function HitPopPoint(a,b, ex, p, cell)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');

C = Clusters{ex}{p};
fprintf('Expt %d Probe %d\n',C.exptno,p);
DATA.selectprobe = zeros(length(Clusters),DATA.nprobes);
DATA.selectprobe(ex,p) = 1;
if nargin > 4 && cell > 0
    DATA.currentcell = cell;
end
ShowData(DATA, ex, p);


 function Clusters = ReloadClusters(DATA, eid)
     
 cfile = [DATA.name '/' DATA.strings{eid}];
 if strfind(cfile,'AutoClusterTimes')
     Clusters = LoadClusters(cfile);
 else
     afile = strrep(cfile,'ClusterTimes','AutoClusterTimes')
     AutoClusters = LoadClusters(afile);
     Clusters = LoadClusters(cfile);
     for j = 1:length(AutoClusters)
         if j > length(Clusters) || ~isfield(Clusters{j},'mahal') || ...
                 (Clusters{j}.auto == 1 && AutoClusters{j}.savetime(1) > Clusters{j}.savetime(1))
             Clusters{j} = AutoClusters{j};
         end
     end
 end
 for j = 1:length(Clusters)
     Clusters{j}.exptid = eid;
 end
 Clusters = FixClusters(Clusters);
     
function rawxy =RawXY(C, xy)
%Convert from XY saved in clusters (rotated for lines)
%Back into the original space
    if C.shape == 0
        rawxy = xy;
    else
        rawxy = xyrotate(xy(:,1),xy(:,2),-C.angle);
    end

function [Clusters, details] = LoadClusters(name)
%loadClusters  loads a cluster set by filename
exptno = -1;

if ~exist(name,'file')
    Clusters = [];
    details = [];
    return;
end
load(name);
details.loadname = name;
details.loadtime = now;

Clusters = Clusters;
for j = 1:length(Clusters)
    if isfield(Clusters{j},'exptno')
        exptno = Clusters{j}.exptno;
    end
    if isfield(Clusters{j},'clst')
        fprintf('%s P%d has clst\n',name,j);
        if length(Clusters{j}.clst) ~= length(Clusters{j}.times)
            Clusters{j} = rmfield(Clusters{j},'clst');
        end
    end
end
if exptno < 0
    id = strfind(name,'Expt');
    if length(id)
        [exptno, a,b,c] = sscanf(name(id(1)+4:end),'%d');
        if name(id(1)+3+c) == 'a'
            exptno = exptno+0.1;
        end
    end
end
details.exptno = exptno;
dname = strrep(name,'.mat','Details.mat');
if ~exist(dname,'file')
dname = strrep(name,'ClusterTimes','ClusterTimesDetails.mat');
end    
if exist(dname)
load(dname);
else
    mycprintf('errors','Can''t find %s\n',dname);
    ClusterDetails = [];
    return;
end
if exist('FullVData','var')
    details.FullVData = FullVData;
end
for j = 1:length(ClusterDetails)
    if isfield(ClusterDetails{j},'xy')
%        Clusters{j}.xy = ClusterDetails{j}.xy;
        Clusters{j}.xy = RawXY(Clusters{j}, ClusterDetails{j}.xy);
        if isfield(ClusterDetails,'next')
            nc = min([length(ClusterDetails{j}.next) length(Clusters{j}.next)]);
        for k = 1:nc
            if isfield(ClusterDetails{j}.next{k},'xy') && ~isempty(Clusters{j}.next{k});
                Clusters{j}.next{k}.xy = ClusterDetails{j}.next{k}.xy;
                Clusters{j}.next{k}.xy = RawXY(Clusters{j}.next{k},ClusterDetails{j}.next{k}.xy);
            end
        end
        end
        if isfield(ClusterDetails{j},'clst')
            Clusters{j}.clst = ClusterDetails{j}.clst;
        else
            Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
            id = find(ismember(ClusterDetails{j}.t,Clusters{j}.times));
            Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
        Clusters{j}.clst(id) = 2;
        end
    Clusters{j}.times = ClusterDetails{j}.t;
    if isfield(Clusters{j},'t')
        Clusters{j} = rmfield(Clusters{j},'t');
    end
    if exptno >= 0
    Clusters{j}.exptno = exptno;
    end
    end
    if isfield(ClusterDetails{j},'Evec')
        Clusters{j}.Evec = ClusterDetails{j}.Evec;
    end
    if isfield(ClusterDetails{j},'triggerV')
        Clusters{j}.triggerV = ClusterDetails{j}.triggerV;
    end
    if isfield(ClusterDetails{j},'next')
        for k = 1:length(ClusterDetails{j}.next)
            if isfield(ClusterDetails{j}.next{k},'xy')
                Clusters{j}.next{k}.xy = ClusterDetails{j}.next{k}.xy;
            end
        end
    end
end

function DATA = LoadTrialLists(DATA)
    Expts = getappdata(DATA.toplevel,'Expts');
    cluster = 1;
    for j = 1:size(Expts,1)
        Tn = [Expts{j,1}.Trials.id];
        tdiff = setxor(DATA.CellDetails.trialids{j},Tn);
        if length(tdiff) > 0 && j <= size(DATA.CellDetails.excludetrials,1)
            fprintf('Expt Trial List Changed (%d) for Expt %d\n',length(tdiff),Expts{j,1}.Header.exptno);
            for k = 1:size(DATA.CellDetails.excludetrials,2)
                nx = length(DATA.CellDetails.excludetrials{j,k,cluster});
                if nx > 0
                    fprintf('Probe %d Cluster %d has %d excluded Trials\n',k,cluster,nx);
                    [a, id] = ismember(Tn,DATA.CellDetails.excludetrials{j,k,cluster});
                end
            end
        end
        DATA.CellDetails.trialids{j} = Tn;
    end
    

function DATA = LoadCellFile(DATA)
    cellfile = [DATA.name '/CellList.mat'];
    Expts = getappdata(DATA.toplevel,'Expts');

    if exist(cellfile,'file')
        load(cellfile);
        d = dir(cellfile);
        dstr = datestr(d.date);
        dstr = dstr([1:2 4:6]);
        nr = min([length(Expts) size(CellList,1)]);
        
%Dangerous to  change size of celllist. It will be written out again and so
%might overwrite datat that is wanted, just becuase and Expt is temporarily
%missing
%        DATA.CellList = CellList(1:nr,:);
        DATA.CellList = CellList;
        DATA.CellDetails = CellDetails;
        if ~isfield(CellDetails,'excludetrials')
            DATA.CellDetails.excludetrials = {};
        end
        if isfield(CellDetails,'Quality')
            CellDetails.Quality = zeros(size(CellList));
            id = find(CellList > 0);
            CellDetails.Quality(id) = 4;
        end
% If all expts not loaded yet, can't check th exptid list.
% exptid = Nan inidicates this has happened
        for j = 1:size(Expts,1)
            if isfield(Expts{j},'Header')
                exptids(j) = Expts{j,1}.Header.exptno;
            else
                exptids(j) = NaN;
            end
        end
        if ~isfield(CellDetails,'exptids')
                    DATA.CellDetails.exptids = exptids;
        elseif exist('exptids','var') && sum(isnan(exptids)) == 0  % check expt list matched
            [id, ia, ib] = setxor(CellDetails.exptids,exptids);
            if length(ib) % new expts
                id = find(ismember(exptids,CellDetails.exptids));
                DATA.CellList(id,:,:) = CellList;
                DATA.CellList(ib,:,:) = 0;
                DATA.CellDetails.exptids = exptids;
                CellDetails = rmfields(CellDetails,'trialids'); %rebuild below
            end
        end
        if ~isfield(CellDetails,'trialids')
            for j = 1:size(Expts,1)
                DATA.CellDetails.trialids{j} = [Expts{j,1}.Trials.id];
            end
        else
            for j = 1:size(Expts,1)
                Tn = [Expts{j,1}.Trials.id];
                [a,b,c] = setxor(DATA.CellDetails.trialids{j},Tn);
                if length(a) > 0 && j <= size(DATA.CellDetails.excludetrials,1)
                    fprintf('Expt Trial List Changed for Expt %d (%d gone, %d new)\n',Expts{j,1}.Header.exptno,length(b),length(c));
                    DATA.CellDetails.trialids{j} = Tn;
                    for k = 1:size(DATA.CellDetails.excludetrials,2)
                        if ~isempty(DATA.CellDetails.excludetrials{j,k,1})
                            [a, id] = ismember(Tn,DATA.CellDetails.excludetrials{j,k,1});
                        end
                    end
                end
            end
        end
        DATA = ConvertExclusion(DATA);
        if exist('CellListB','var')
            DATA.CellList(1:size(CellListB,1),1:size(CellListB,2),2) = CellListB;
        end
        if exist('CellChanges','var')
            DATA.CellChanges = CellChanges;
            if size(DATA.CellChanges,2) < 6
                DATA.CellChanges(:,6) = NaN;
            end
        else
            CellChanges = [];
        end
        if DATA.cellbackup == 0 %not backed up, so do this now
            bakfile = strrep(cellfile,'.mat',[dstr '.mat']);
            try
                save(bakfile,'CellList','CellDetails','CellChanges');
            catch
                cprintf('red','COuld not back up Cell list to %s\n',bakfile);
            end
            DATA.cellbackup = 1;
            set(DATA.toplevel,'UserData',DATA);
        end
    elseif isfield(DATA,'exptid')
        DATA.CellList = zeros([length(DATA.exptid) DATA.nprobes 2]);
        DATA.CellDetails.exptids = DATA.exptid;
    end


function nt = PlotTrialSpikes(DATA, nt, varargin)
if length(nt) > 1
    for j = 1:length(nt)
        k= PlotTrialSpikes(DATA, nt(j), varargin{:});
        drawnow;
    end
    nt = k;
    return;
end


if nt < 1 
    nt = 1;
    return;
end

if  nt > length(DATA.Expt.Trials)
    nt = length(DATA.Expt.Trials);
    return;
end
Clusters = getappdata(DATA.toplevel,'Clusters');
e = DATA.currentpoint(1);
AllSpikes = CheckAllSpikes(DATA,e, 1:DATA.nprobes);
SpoolAllProbes(DATA, e, AllSpikes(e,:), Clusters{e}, 'Trials' , nt)

    
function SelectTrial(src, b, type)
DATA = GetDataFromFig(src);
id = get(src,'value');
if strcmp(type,'one')
    SpoolSpikes(DATA,DATA.currentpoint,'ids',id)
else
    DATA.currenttrial = PlotTrialSpikes(DATA,id,'showall');
end
set(DATA.toplevel,'UserData',DATA);


function PlotTimeRanges(DATA)

    X = getappdata(DATA.toplevel,'Clusters');
    Expts = getappdata(DATA.toplevel,'Expts');
    colors = mycolors;
    Clusters = X{DATA.currentpoint(1)};
    e = DATA.currentpoint(1);
    SetFigure(DATA,DATA.tag.checktimes);
    hold off;
    for j = 1:length(Clusters)
        C = Clusters{j};
        plot(Clusters{j}.times,'color',colors{j});
        text(length(C.times),C.times(end),sprintf('%d',j),'color',colors{j});
        hold on;
        xcl = FindExcludedTrials(DATA, e, j,1, C);
        if length(xcl)
            t = [Expts{e,j}.Trials(xcl).TrialStart]./10000;
            plot(length(C.times),t,'.','color',colors{j});
        end
    end

function ChangeTag(a,b,fcn)
    DATA = GetDataFromFig(a);

    x = get(a,'parent');
    while ~isfigure(x)
        x = get(x,'parent');
    end
    set(x,'tag',[fcn 'expt' num2str(DATA.currentpoint(1))],'name',[fcn 'expt' num2str(DATA.currentpoint(1))]);
    
    
function [f, isnew] = SetFigure(DATA, tag, varargin)
    onoff = {'off' 'on'};
    [f,isnew] = GetFigure(tag,varargin{:});
    if isnew
        if strfind(tag,'AllXY');
            set(f,'UserData', DATA.toplevel);
            hm = uimenu(f,'Label','&Options','Tag','Options');
            uimenu(hm,'Label','&Keep This Graph','callback', {@ChangeTag, 'allxy'});
            uimenu(hm,'Label','Show &Means','callback', {@PlotMenu, 'cells' 'allmean'});
            uimenu(hm,'Label','Show &XY','callback', {@PlotMenu, 'cells' 'allxy'});
            uimenu(hm,'Label','&Density','callback', {@OptionMenu, 1});
            uimenu(hm,'Label','&Optimize all','callback', {@OptionMenu, 'optimizeall'});
            uimenu(hm,'Label','Add Context Menus','callback', {@PlotMenu, 'cells' 'addcontext'},'checked',onoff{1+DATA.show.linecontextmenus});
        elseif strfind(tag,'Comments')
            x = get(DATA.toplevel,'Position');
            set(f,'position',[x(1) x(2) 300 50],'menubar','none');
            uicontrol(f, 'string',[],'style','edit','units','normalized','position',[0.05 0.55 0.9 0.45],...
                'Callback',{@AddComment, []});
            uicontrol(f, 'string','Done','style','pushbutton','units','normalized','position',[0.05 0.05 0.2 0.45],...
                'Callback',{@AddComment, []});
            uicontrol(f, 'string','Done','style','text','units','normalized','position',[0.25 0.05 0.7 0.45],...
                'Tag','ProbeLabel');
            
            set(f,'UserData',DATA.toplevel);
        elseif strfind(tag,'CellSummary');
            hm = uimenu(f,'Label','Cell','Tag','CellList');
            set(f,'UserData', DATA.toplevel);
        elseif strfind(tag,'CellList');
            hm = uimenu(f,'Label','&Plots','Tag','PlotMenu','accelerator','p');
            sm = uimenu(hm,'Label','&Cells','Tag','ProbePlotMenu','accelerator','c');
            AddPlotMenu(sm,'cells');
            sm = uimenu(hm,'Label','&Probes','Tag','ProbePlotMenu','accelerator','p');
            AddPlotMenu(sm,'probes');
            uimenu(hm,'Label','AllXY','Callback',{@PlotMenu, 'cells' 'allxy'});
            uimenu(hm,'Label','All Means','Callback',{@PlotMenu, 'cells' 'allmean'});
            uimenu(hm,'Label','All Mean image','Callback',{@PlotMenu, 'cells' 'allmean'});
            uimenu(hm,'Label','All Trigger hist','Callback',{@PlotMenu, 'cells' 'trighist'});
            uimenu(hm,'Label','All Spks','Callback',{@PlotMenu, 'cells' 'allspks'});
            uimenu(hm,'Label','All Spks when change cell','Callback',{@PlotMenu, 'cells' 'watchallcellspks'},'tag','watchallcellspks');
            uimenu(hm,'Label','All Mean when change cell','Callback',{@PlotMenu, 'cells' 'watchallcellmean'},'tag','watchallcellmean');
            uimenu(hm,'Label','All XY when change cell','Callback',{@PlotMenu, 'cells' 'watchallcellxy'},'tag','watchallcellxy');
            uimenu(hm,'Label','SpkRate','Callback',{@PlotMenu, 'cells' 'spkrate'});
            uimenu(hm,'Label','SpkRate+xy sequence','Callback',{@PlotMenu, 'cells' 'spkrate+xy'});
            uimenu(hm,'Label','All','Callback',{@PlotMenu, 'cells' 'allplots'});
            uimenu(hm,'Label','All - next cell','Callback',{@PlotMenu, 'cells' 'allnext'});
            uimenu(hm,'Label','All Cell &Tuning This Expt','Callback',{@PlotMenu, 'cells' 'allexpttune'});
            uimenu(hm,'Label','All Tuning This Cell','Callback',{@PlotMenu, 'cells' 'allcelltune'});
            
            uimenu(hm,'Label','Estimate Drift','Callback',{@PlotMenu, 'cells' 'driftestimate'});
            uimenu(hm,'Label','Cell Summary','Callback',{@PlotMenu, 'cells' 'summary'});
            uimenu(hm,'Label','->AllVpcs','Callback',{@GuiMenu, 'CallAllVPcs'},'accelerator','a');
            
            hm = uimenu(f,'Label','&Cells','Tag','CellMenu');
            uimenu(hm,'Label','&Fill Using Marked','Callback',{@OptionMenu 'fillcellsfrommark'});
            uimenu(hm,'Label','Check &Rate sequences','Callback',{@OptionMenu 'checkratesequence'});
            uimenu(hm,'Label','Calculate Means','Callback',{@OptionMenu 'calccellmeans'});
            
            hm = uimenu(f,'Label','&Mark','Tag','MarkMenu');
            uimenu(hm,'Label',sprintf('Dropi < %.1f',DATA.crit.dropi),'Callback',{@PlotMenu, 'cells' 'markdropi'},'Tag','markdropi');
            uimenu(hm,'Label','Candidates','Callback',{@PlotMenu, 'cells' 'markcandidates'});
            uimenu(hm,'Label','mahal < 3','Callback',{@PlotMenu, 'cells' 'markmahal'});
            uimenu(hm,'Label','Ellipses','Callback',{@PlotMenu, 'cells' 'markellipse'});
            uimenu(hm,'Label','&Quick','Callback',{@PlotMenu, 'cells' 'markquick'});
            uimenu(hm,'Label','Old readmethod','Callback',{@PlotMenu, 'cells' 'markreadmethod'});
            uimenu(hm,'Label','Tagged','Callback',{@PlotMenu, 'cells' 'marktagged'});
            uimenu(hm,'Label','Good Cell','Callback',{@PlotMenu, 'cells' 'markgoodcell'});
            uimenu(hm,'Label','Good MU','Callback',{@PlotMenu, 'cells' 'markgoodmu'});
            uimenu(hm,'Label','Grid','Callback',{@PlotMenu, 'cells' 'markgrid'});
            uimenu(hm,'Label','Expt &Names','Callback',{@PlotMenu, 'cells' 'markexpnames'});
            sm = uimenu(hm,'Label','&Expts','Tag','ExptMarkMenu');
            AddExptList(sm,'markexpts',DATA)
            sm = uimenu(hm,'Label','&Other','Tag','OtherPlotMenu');
            uimenu(hm,'Label','Shape/Size','Callback',{@PlotMenu, 'cells' 'shapesize'});

            
            bp = [0.01 0.001 0.15 0.04];
            uicontrol(gcf,'style','pop','string',num2str([1:50]'), ...
        'Tag','CellNumberId','Callback', {@SetCellNumber, 'change'},...
        'units', 'norm', 'position',bp,'value',DATA.currentcell);
            bp = [0.1 0.001 0.08 0.04];
            it = uicontrol(gcf,'style','pushbutton','string','Set', ...
                'Callback', {@SetCellNumber, 'set'}, 'Tag','SetCell', ... 
                'units', 'norm', 'position',bp,'value',1);
            bp(1) = bp(1)+bp(3);
            bp(3) = 0.16;
            uicontrol(gcf,'style','pop','string','Just This Square|Cluster 2 This Square|Cluster 3|Cluster 4 This Square|Use Template Peaks|Better than mahal|Quality-suboptimal|Quality-poor|Spool All Expts|AllXY|AllMean|AllMeanIm|Exclude Selected Trials|Exclude Trials for C2', ...
        'Tag','ClusterModifier','callback',{@CellTrackMenu, 'setcluster'},...
                'units', 'norm', 'position',bp,'value',1);
            bp(1) = bp(1)+bp(3);
            bp(3) = 0.08;
            uicontrol(gcf,'style','pop','string','1+2|3+4|5+6', ...
        'Tag','CellPlane','callback',{@CellTrackMenu, 'setplane'},...
        'units', 'norm', 'position',bp,'value',1);
            bp(1) = bp(1)+bp(3);
            uicontrol(gcf,'style','pushbutton','string','Delete', ...
                'Callback', {@SetCellNumber, 'delete'}, 'Tag','DeleteCell',...
                'units', 'norm', 'position',bp,'value',1);
            bp(1) = bp(1)+bp(3);
            uicontrol(gcf,'style','pop','string','1|2|3|4|5|6', ...
        'Tag','CellCluster','callback',{@CellTrackMenu, 'setcluster'},...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3);
    bp(3) = 0.2;
    it = uicontrol(gcf,'style','pushbutton','string','Set Duplicate', ...
        'Callback', {@SetCellNumber, 'duplicates'}, 'Tag','SetDuplicate',...
        'Keypressfcn',{@KeyPressed, 3},...
        'units', 'norm', 'position',bp,'value',1);

    
    set(f,'UserData', DATA.toplevel);
        elseif strcmp(tag,DATA.tag.allspikes)
            Expts = getappdata(DATA.toplevel,'Expts');
             DATA.Expt = Expts{DATA.currentpoint(1),1};
            bp = [0.95 0.05 0.05 0.9];
            h = uicontrol(f,'Style', 'list',...
                'String', num2str([DATA.Expt.Trials.Trial]'), 'Tag', 'ChooseTrial','Units','norm', 'Position', bp,'value',1,...
                'Max',3,'Min',1,'Callback',{@SelectTrial, 'all'});
            bp = [0.95 0.95 0.05 0.05];
            uicontrol(f,'style','check','string','stop','Units','Normalized','Position',bp,'Tag','StopSpool',...
            'value',0);
        hm = uimenu(f,'Label','Spool','Tag','Cluster');
        uimenu(hm,'Label','Selected','Callback',{@PlotMenu, 'probes', 'spoolcurrent'});
        uimenu(hm,'Label','Synchronous','Callback',{@PlotMenu, 'probes', 'plotsync'});
        uimenu(hm,'Label','Synchronous Flip','Callback',{@PlotMenu, 'probes', 'plotsyncb'});
        hm = uimenu(f,'Label','Options','Tag','Cluster');
        uimenu(hm,'label','Keep','callback', {@ChangeTag, 'allspikes'});
        uimenu(hm,'label','Ues one scale','callback', {@PlotMenu, 'probes' 'onescale'});
        uimenu(hm,'Label','Add Context Menus','callback', {@PlotMenu, 'cells' 'addcontext'},'checked',onoff{1+DATA.show.linecontextmenus});

             set(f,'UserData',DATA.toplevel);
            
        elseif strfind(tag,'ClustersAll');
            hm = uimenu(f,'Label','Mark','Tag','MarkMenu');
            uimenu(hm,'Label','Manual Cut','Callback',{@PlotMenu, 'allcluster' 'markmanual'});
            uimenu(hm,'Label','Auto Cut','Callback',{@PlotMenu, 'allcluster' 'markauto'});
            uimenu(hm,'Label','Copied','Callback',{@PlotMenu, 'allcluster' 'markcopy'});
            uimenu(hm,'Label','Marked Good','Callback',{@PlotMenu, 'allcluster' 'markgood'});
            uimenu(hm,'Label','Bad Probes','Callback',{@PlotMenu, 'allcluster' 'markbadprobe'});
            uimenu(hm,'Label','Good Mu','Callback',{@PlotMenu, 'allcluster' 'markgoodmu'});
            uimenu(hm,'Label','Excluded Trials','Callback',{@PlotMenu, 'allcluster' 'markexclusions'});
            uimenu(hm,'Label','Quick','Callback',{@PlotMenu, 'cells' 'markquick'});
            set(f,'UserData', DATA.toplevel);
        elseif strcmp(tag,DATA.tag.xcorr)
            set(f,'UserData', DATA.toplevel);
            optm = uimenu(f,'Label','Options','Tag','XcorrOptions');
            hm = uimenu(optm,'Label','Zoom','Tag','xcZoom');
            sm = uimenu(hm,'Label','Main plot');
            uimenu(sm,'Label','10ms','Callback',{@PlotMenu, 'xcorr', 'zoom', 10});
            uimenu(sm,'Label','50ms','Callback',{@PlotMenu, 'xcorr', 'zoom', 50});
            uimenu(sm,'Label','100ms','Callback',{@PlotMenu, 'xcorr', 'zoom', 100});
            uimenu(sm,'Label','200ms','Callback',{@PlotMenu, 'xcorr', 'zoom', 200});
            sm = uimenu(hm,'Label','Sync Spikes');
            uimenu(sm,'Label','Off','Callback',{@PlotMenu, 'xcorr', 'spkzoom', 0});
            uimenu(sm,'Label','2ms','Callback',{@PlotMenu, 'xcorr', 'spkzoom', 2});
            uimenu(sm,'Label','5ms','Callback',{@PlotMenu, 'xcorr', 'spkzoom', 5});
            uimenu(sm,'Label','10ms','Callback',{@PlotMenu, 'xcorr', 'spkzoom', 10});
            uimenu(sm,'Label','20ms','Callback',{@PlotMenu, 'xcorr', 'spkzoom', 20});
            hm = uimenu(optm,'Label','&Plot');
            uimenu(hm, 'label','Shape/Efficacy','CallBack', {@ReplotXcorrs, 'Shape/efficacy'});
            uimenu(hm, 'label','Efficacy Image','CallBack', {@ReplotXcorrs, 'Efficacy'});
            sm = uimenu(hm,'Label','&One Pair Type');
            uimenu(sm, 'label','&CorssCorr','CallBack', {@ReplotXcorrs, 'xcorr'});
            uimenu(sm, 'label','&Mean Image','CallBack', {@ReplotXcorrs, 'meanim'});
            uimenu(sm, 'label','&Sync Spikes','CallBack', {@ReplotXcorrs, 'syncspikes'});
            uimenu(sm, 'label','&Histograms','CallBack', {@ReplotXcorrs, 'histograms'});

            
            hm = uimenu(optm,'Label','Labels','Tag','xcLabels');
            uimenu(hm,'Label','Cell# on left','Callback',{@PlotMenu, 'xcorr', 'labelcellleft', 0});
            uimenu(hm,'Label','Cell# on top','Callback',{@PlotMenu, 'xcorr', 'labelcelltop', 0});
            hm = uimenu(optm,'Label','replot','Tag','xcreplot','Callback',{@PlotMenu, 'xcorr', 'replot', 0});
            hm = uimenu(f,'Label','Exclude','Tag','Xcorrexclude');
            cellid = unique(DATA.CellList);
            cellid = cellid(cellid > 0);
            for j = 1:length(cellid)
                uimenu(hm,'Label',['Cell ' num2str(cellid(j))],'Callback',{@PlotMenu, 'xcorr', 'exclude', cellid(j)});
            end
            setappdata(f,'ParentFigure', DATA.toplevel);
            set(gcf, 'KeyPressFcn',{@KeyPressed,4},'Keyreleasefcn',{@KeyReleased, 4});
        elseif strfind(tag,'Spikes')
            Expts = getappdata(DATA.toplevel,'Expts');
             DATA.Expt = Expts{DATA.currentpoint(1),1};
             DATA.Expt = Expts{DATA.currentpoint(1),1};
            bp = [0.9 0.05 0.1 0.9];
            h = uicontrol(f,'Style', 'list',...
                'String', num2str([DATA.Expt.Trials.Trial]'), 'Tag', 'ChooseTrial','Units','norm', 'Position', bp,'value',1,...
                'Max',3,'Min',1,'Callback',{@SelectTrial, 'one'});
            bp = [0.9 0.95 0.05 0.05];
            uicontrol(f,'style','check','string','stop','Units','Normalized','Position',bp,'Tag','StopSpool',...
            'value',0);
            bp = [0.01 0.01 0.06 0.05];
            uicontrol(f,'style','pushbutton','string','AllEx','Units','Normalized','Position',bp,'Tag','AllExptButton',...
            'value',0,'callback',{@PlotMenu, 'probes', 'allprobespks'});
            bp(1) = bp(1)+bp(3);
            bp(3)=0.08;
            uicontrol(f,'style','pushbutton','string','AllProbes','Units','Normalized','Position',bp,'Tag','AllExptButton',...
            'value',0,'callback',{@PlotMenu, 'probes', 'allspks'});
            

            set(f,'UserData',DATA.toplevel);
            hm = uimenu(f,'Label','&Spool','Tag','Spool');
            uimenu(hm,'Label','This Probe','Callback',{@PlotMenu, 'spikes', 'Spool'});
            uimenu(hm,'Label','Sync Spikes','Callback',{@PlotMenu, 'probes', 'spoolsync'});
            uimenu(hm,'Label','&Exclude Selected Trials (current cluster)','Callback',{@OptionMenu, 'excludecurrent'});
        elseif strfind(tag,'XYplot');
                hm = uimenu(f,'Label','&Cluster','Tag','Cluster');
                uimenu(hm,'Label','&SaveChanges','Callback',{@XYCluster, 'save'});
                uimenu(hm,'Label','&QuickSave (q)','Callback',{@XYCluster, 'quicksave'});
%                uimenu(hm,'Label','Cut Ellipses','Callback',{@XYCluster, 'ellipses'});
                sm = uimenu(hm,'Label','Cut &Ellipse','Tag','EllipseMenu');
                uimenu(sm,'Label','Cluster &1','foregroundcolor',DATA.colors{2},'Callback',{@XYCluster, 'ellipse1'},'accelerator','1');
                uimenu(sm,'Label','Cluster &2','foregroundcolor',DATA.colors{3},'Callback',{@XYCluster, 'ellipse2'},'accelerator','2');
                uimenu(sm,'Label','Cluster &3','foregroundcolor',DATA.colors{4},'Callback',{@XYCluster, 'ellipse3'},'accelerator','3');
                uimenu(sm,'Label','Cluster &4','foregroundcolor',DATA.colors{5},'Callback',{@XYCluster, 'ellipse4'},'accelerator','4');
                sm = uimenu(hm,'Label','&View Cluster','Tag','EllipseMenu');
                uimenu(sm,'Label','Cluster &1','foregroundcolor',DATA.colors{2},'Callback',{@XYCluster, 'select1'});
                uimenu(sm,'Label','Cluster &2','foregroundcolor',DATA.colors{3},'Callback',{@XYCluster, 'select2'});
                uimenu(sm,'Label','Cluster &3','foregroundcolor',DATA.colors{4},'Callback',{@XYCluster, 'select3'});
                uimenu(sm,'Label','Cluster &4','foregroundcolor',DATA.colors{5},'Callback',{@XYCluster, 'select4'});
                sm = uimenu(hm,'Label','&Delete','Tag','EllipseMenu');
                uimenu(sm,'Label','Cluster &2','foregroundcolor',DATA.colors{3},'Callback',{@XYCluster, 'clear2'});
                uimenu(sm,'Label','Cluster &3','foregroundcolor',DATA.colors{4},'Callback',{@XYCluster, 'clear3'});
                uimenu(sm,'Label','Cluster &4','foregroundcolor',DATA.colors{5},'Callback',{@XYCluster, 'clear4'});
                uimenu(hm,'Label','Cut &Lines','Callback',{@XYCluster, 'lines'},'accelerator','L');
                uimenu(hm,'Label','Flip Line Criterion (f)','Callback',{@XYCluster, 'flip'},'accelerator','F');
 %               uimenu(hm,'Label','Ellipse for Cl2','Callback',{@XYCluster, 'ellipse2'});
                uimenu(hm,'Label','&No Cutting','Callback',{@XYCluster, 'nocut'});
                uimenu(hm,'label','TwoCluster Space','Callback',{@XYCluster, 'twoplot'});
                uimenu(hm,'Label','Optimize angle','Callback',{@XYCluster, 'bestangle'});
                uimenu(hm,'Label','Optimize angle (quick)','Callback',{@XYCluster, 'quickangle'});
                uimenu(hm,'Label','replot','Callback',{@XYCluster, 'replot'});
                uimenu(hm,'Label','revert','Callback',{@XYCluster, 'revert'});
                uimenu(hm,'Label','replot with GM ids','Callback',{@XYCluster, 'replotgm'});
                uimenu(hm,'Label','&Density','Callback',{@OptionMenu, 1});
                uimenu(hm,'Label','Scale Density','Callback',{@OptionMenu, 'scaledensity'});
                uimenu(hm,'Label','Use Auto','Callback',{@OptionMenu, 'useautoclusters'});
                bp = [0.01 0.01 0.06 0.05];
                uicontrol(f,'style','pushbutton','string','AllEx','Units','Normalized','Position',bp,'Tag','AllExptButton',...
            'value',0,'callback',{@PlotMenu, 'probes', 'AllprobeXY'});
            bp(1) = bp(1)+bp(3);
            bp(3) = 0.1;
            uicontrol(f,'style','pushbutton','string','AllProbes','Units','Normalized','Position',bp,'Tag','AllExptButton',...
            'value',0,'callback',{@PlotMenu, 'probes', 'AllXY'});
                tmp.parentfig = DATA.toplevel;
                set(f, 'KeyPressFcn',{@XYKeyPressed,3},'Keyreleasefcn',{@KeyReleased, 3});
                set(f,'UserData',DATA.toplevel);
        elseif strcmp(tag,DATA.tag.rateseq)
                hm = uimenu(f,'Label','&Options','Tag','Options');
                uimenu(hm,'Label','&Exclude Selected Trials (current cluster)','Callback',{@OptionMenu, 'excludecurrent'});
                uimenu(hm,'Label','&Scroll This Expt/Probe','Callback',{@RateMenu, 'spool'});
                uimenu(hm,'Label','&XY sequence plot','Callback',{@RateMenu, 'xyseq'});
                uimenu(hm,'Label','Zoom Out &Horizontal (-)','Callback',{@RateMenu, 'hzoomout'});
                uimenu(hm,'Label','Zoom Out &Vertical','Callback',{@RateMenu, 'vzoomout'});
                uimenu(hm,'Label','Zoom &Out','Callback',{@RateMenu, 'zoomout'});
                uimenu(hm,'Label','Zoom In Horizontal (&+)','Callback',{@RateMenu, 'hzoomin'});
                uimenu(hm,'Label','Zoom &In (Shift+)','Callback',{@RateMenu, 'zoomin'});
                sm = uimenu(hm,'Label','&Expts','Tag','ExptMarkMenu');
                AddExptList(sm, 'selectexpts', DATA)
                X.toplevel = DATA.toplevel;
                set(f,'UserData',X);
        elseif strcmp(tag,DATA.tag.xyseq)
                hm = uimenu(f,'Label','&Options','Tag','Options');
                uimenu(hm,'Label','&Density','Callback',{@PlotMenu, 'xyseq', 'density'});
                set(f,'UserData', DATA.toplevel);
        elseif strcmp(tag,DATA.tag.clusters)
            setappdata(f,'ParentFig', DATA.toplevel);
            hm = uimenu(f,'Label','&Plot','Tag','Plot');
            sm = uimenu(hm,'Label','RateCheck');
            uimenu(sm,'Label','FanoF-C.V.','Callback',{@RatePlotMenu, 'cv'});
            uimenu(sm,'Label','FanoF-Skew','Callback',{@RatePlotMenu, 'skew'});
            uimenu(sm,'Label','Image','Callback',{@RatePlotMenu, 'image'});
            sm = uimenu(hm,'Label','XY/Spks','Tag','TrialMenu');
            uimenu(sm,'Label','SelectedXY+Spks','Callback',{@PlotMenu, 'clusters', 'spksxy'});
            uimenu(sm,'Label','SelectedXY','Callback',{@PlotMenu, 'probes', 'SelectXY'});
            uimenu(sm,'Label','SelectedSpks','Callback',{@PlotMenu, 'probes', 'selectspks'});
       else
            set(f,'UserData', DATA.toplevel);
        end
        setappdata(f,'ParentFigure',DATA.toplevel);
        Figpos = getappdata(DATA.toplevel,'Figpos');
        if isfield(Figpos,tag)
            set(f,'position',Figpos.(tag));
        end
    end
    

function DATA = ReFitAll(DATA, fittype)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    Expts = getappdata(DATA.toplevel,'Expts');

    if strcmp(fittype,'sdindex')
        DATA.xysdindex = [];
    end
    ts = now;
    for j = 1:length(Clusters)
        for k = 1:length(Clusters{j})
            DATA.Expt = Expts{j,k};
            C = Clusters{j}{k};
            if strcmp(fittype,'3means') 
            a = GMfit(C.xy, 3,1);
            DATA.mahal3(j,k) = GMdprime(a);
            elseif strcmp(fittype,'sdindex')
                sdx = PlotXYSequence(DATA, [j k],'noplot');
                DATA.xysdindex(j,k) = sdx(end);
            end
        end
        fprintf('Ex%d %.2f\n',j,mytoc(ts));
    end
    SaveExtras(DATA);
    set(DATA.toplevel,'UserData',DATA);


function ReFit1D(DATA, ratio)
 
     Clusters = getappdata(DATA.toplevel,'Clusters');
     [row,col]  = find(DATA.mahal(:,:,1) > 2 & DATA.mahal(:,:,1)./DATA.mahal(:,:,3) > ratio);
     for j = 1:length(row)
         C = Clusters{row(j)}{col(j)};
         [a,b,c] = BestAngleGM(C.xy, C.gmfit2d, []);
         DATA.mahal(row(j),col(j),4) = b;
         if b < DATA.mahal(row(j),col(j),3);
             DATA.mahal(row(j),col(j),4) = b;
         end
     end
     SetFigure(DATA,DATA.tag.clusters);
     hold off;
     for j = 1:length(row)
         plot(DATA.mahal(row(j),col(j),3),DATA.mahal(row(j),col(j),4),'o','buttondownfcn',{@HitPopPoint, row(j),col(j)});
         hold on;
     end
     refline(1);

 function DATA = ReFit3means(DATA, ratio)
 
     Clusters = getappdata(DATA.toplevel,'Clusters');
     for j = 1:length(Clusters)
         for k = 1:length(Clusters{j})
         C = Clusters{j}{k};
         a = GMfit(C.xy, 3,1);
         DATA.mahal3(j,k) = GMdprime(a);
         end
     end
     SetFigure(DATA,DATA.tag.clusters);
     hold off;
     for j = 1:length(Clusters)
         for k = 1:length(Clusters{j})
         plot(DATA.mahal(j,k,1),DATA.mahal3(j,k),'o','buttondownfcn',{@HitPopPoint, row(j),col(j)});
         hold on;
         end
     end
     refline(1);
     SaveExtras(DATA);
     set(DATA.toplevel,'UserData',DATA);
   
 function DATA = ReFitGMDip(DATA, varargin)
 
     ts = now;
     Clusters = getappdata(DATA.toplevel,'Clusters');
     cid = 1:length(Clusters);
     j = 1;
     while j <= length(varargin)
         if strncmpi(varargin{j},'expts',5)
             j = j+1;
             cid = varargin{j};
         end
         j = j+1;
     end
     warning('off','stats:gmdistribution:FailedToConverge');
     warning('off','stats:gmdistribution:MaxIterations');

     for j = cid
         for k = 1:length(Clusters{j})
         C = Clusters{j}{k};
         [a,b] = GMDip(C.xy(:,1),[],'idlist',C.clst,'noplot');
         DATA.GaussFitdp(j,k,1) = b.gmdprime;
         [DATA.GaussFitdp(j,k,2), b, c] = Fit2Gauss(C);
         DATA.GaussFitdp(j,k,3) = now; %record fit time so can update
         DATA.gmfitpos(j,k,:) = c.fitpos;
         fprintf('E%dP%d %.2f %.2f %d %d\n',j,k,DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),DATA.gmfitpos(j,k,1),DATA.gmfitpos(j,k,2));
         end
         fprintf('Expt %d took %.0f sec\n',j,mytoc(ts));
     end
     mytoc(ts);
     SetFigure(DATA,DATA.tag.clusters);
     hold off;
     for j = 1:length(Clusters)
         for k = 1:length(Clusters{j})
             if sum(DATA.gmfitpos(j,k,:)) == 2
         plot(DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),'o','buttondownfcn',{@HitPopPoint, j,k});
             else
         plot(DATA.GaussFitdp(j,k,1),DATA.GaussFitdp(j,k,2),'ro','buttondownfcn',{@HitPopPoint, j,k});
             end
         hold on;
         end
     end
     refline(1);
     SaveExtras(DATA);
     set(DATA.toplevel,'UserData',DATA);

 function DATA = LoadExtra(DATA, force)
     Xbits = [];
     outname = [DATA.name '/PlotClusterExtra.mat'];
     if exist(outname,'file')
         load(outname);
         f = fields(Xbits);
         for j = 1:length(f)
             if ~isfield(DATA,f{j}) || force
             DATA.(f{j}) = Xbits.(f{j});
             end
         end
     end

 function [CC, GM] = CondenseClusters(C)
     for j = 1:length(C)
         [CC{j}, GM{j}] = CondenseCluster(C{j});
     end
     
 function [CC, GM] = CondenseCluster(C)
     CC = C;
     GM = {};
     for p = 1:length(C)
         f = fields(C{p});
         for k = 1:length(f)
             if isobject(C{p}.(f{k}))
                 GM{p}.(f{k}) = C{p}.(f{k});
                 CC{p} = rmfield(CC{p},f{k});
             end
         end
         if isfield(CC{p},'next')
             for c = 1:length(CC{p}.next)
                 if ~isempty(CC{p}.next{c})
                     f = fields(C{p}.next{c});
                     for k = 1:length(f)
                         if isobject(C{p}.next{c}.(f{k}))
                             GM{p}.next{c}.(f{k}) = C{p}.next{c}.(f{k});
                             CC{p}.next{c} = rmfield(CC{p}.next{c},f{k});
                         end
                     end
                 end
             end
         end
     end


     
 function SaveExtras(DATA)
     Xbits = [];
     outname = [DATA.name '/PlotClusterExtra.mat'];
     if isfield(DATA,'GaussFitdp') 
         Xbits.GaussFitdp = DATA.GaussFitdp;
         Xbits.gmfitpos = DATA.gmfitpos;
     end
     if isfield(DATA,'mahal3')
         Xbits.mahal3 = DATA.mahal3;
     end
     if isfield(DATA,'xcorrs')
         Xbits.xcorrs= DATA.xcorrs;
         Xbits.xcorrval = DATA.xcorrval;
     end
     if isfield(DATA,'xysdindex')
         Xbits.xysdindex = DATA.xysdindex;
     end
     names = {'shiftmatrix' 'fitjumps'};
     for j = 1:length(names)
     if isfield(DATA,names{j})
         Xbits.(names{j}) = DATA.(names{j});
     end
     end
     if ~isempty(Xbits)
         save(outname,'Xbits')
     end

function [dp, fits, details] = Fit2Gauss(C, r, varargin)

    plottype = 0;
    id = find(C.clst ==2);
    nid = find(C.clst ==1);
    if C.shape == 0 %ellipse
        rx = C.xyr(3);
        ry = C.xyr(4);
        xys = xyrotate(C.xy(:,1)-C.xyr(1),C.xy(:,2)-C.xyr(2),C.angle);
        r = sqrt(((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2);
        fx = linspace(min(r),max(r),200);
        [y,x] = hist(r,fx);
        [a,b] = min(abs(fx-1));
    elseif isfield(C,'crit')
        xys = xyrotate(C.xy(:,1),C.xy(:,2),C.angle);
        fx = linspace(min(xys(:,1)),max(xys(:,1)),200);
        [a,b] = min(abs(fx-C.crit));
    [y,x] = hist(xys(:,1),fx);
    else
        [a,b] = min(abs(fx-mean(fx))); %temporary
        [y,x] = hist(C.xy(:,1),fx);
    end
    fits{1} = FitGauss(x(1:b),y(1:b));
    fits{2} = FitGauss(x(b:end),y(b:end));
     if isfield(fits{1},'params') && isfield(fits{2},'params')
         details.fitpos = [fits{1}.params(1) < fx(b) fits{2}.params(1) > fx(b)];
         dp = abs((fits{1}.params(1)-fits{2}.params(1))./sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2));
         details.diff = abs((fits{1}.params(1)-fits{2}.params(1)));
         details.sd = sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
         if plottype == 1
             SetFigure(DATA,DATA.tag.fitgauss);
             hold off;
             bar(x,y);
             hold on;
             bar(x(1:b),y(1:b),'r')
             fya = fitGauss(fx, fits{1}.params, 'eval');
             fyb = fitGauss(fx, fits{2}.params, 'eval');
             plot(fx,fya,'r-','linewidth',2);
             plot(fx,fyb,'b-','linewidth',2);
             title(sprintf('Dprime from fits %.2f',dp));
         end
     else
         dp = NaN;
         details.fitpos = [0 0];
         details.diff = NaN;
         details.sd = NaN;
     end
    
 function [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)
%
% Find best angle to cut in a 2D space, using 1-D Gaussian Mixtures
% if quickmode ==1, uses a 2-D fit to the data, clusters, then finds the
% angle that maximize s dptime between the cluters. 
  quickmode = 0;
  a = 0:pi/36:pi;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'quick',5)
        quickmode = 1;
    end
    j = j+1;
end
  if isempty(G)
      G = GMFit(xy,2,1);
  end
  id = cluster(G, xy);
  cid = find(id > 1);
  nid = find(id == 1);
  if isempty(dipfit)
      [dip, dipfit] = GMDip(xy,[],'idlist',id);
  end
  if quickmode
      bestfit = [];
      if dipfit.sign >= 0
          acid = xy(:,1) > dipfit.dip(1);
          anid = xy(:,1) < dipfit.dip(1);
      else
          acid = xy(:,1) < dipfit.dip(1);
          anid = xy(:,1) > dipfit.dip(1);
      end
  end
  
  c = 0;
  for j = 1:length(a)
      XY = xyrotate(xy(:,1),xy(:,2),a(j));
      if quickmode
          details.d(j) = abs(CalcDprime(XY(cid,1),XY(nid,1)));
          details.oned(j) = abs(CalcDprime(XY(acid,1),XY(anid,1)));
      else
          [aa,bb] = GMDip(XY,[],'idlist',id);
          details.d(j) = bb.mahal(bb.best);
          details.dprime(j) = CalcDprime(XY(cid,1),XY(nid,1));
          if details.d(j) > c
              c = details.d(j);
              bestfit = bb.G{bb.best};
              crit = mean(bb.dip);
          end
      end
  end
  if quickmode
      [cs(1), js(1)] = max(details.oned);
      XY = xyrotate(xy(:,1),xy(:,2),a(js(1)));
      [aa,bs{1}] = GMDip(XY,[],'idlist',id);
      q(1) = bs{1}.mahal(bs{1}.best);
      [cs(2), js(2)] = max(details.d);
      XY = xyrotate(xy(:,1),xy(:,2),a(js(2)));
      [bb,bs{2}] = GMDip(XY,[],'idlist',id);
      q(2) = bs{2}.mahal(bs{2}.best);
      q(3) = dipfit.mahal(dipfit.best);
      js(3) = 1;
      bs{3} = dipfit;
      [aa,b] = max(q);
      c = q(b); %best mahal distance
      theta = a(js(b));
      bestfit = bs{b}.G{bs{b}.best};
      details.besti = js(b);
      details.crit = mean(bs{b}.dip)
  else
  [c, j] = max(details.d);
  details.besti = j;
  theta = a(j);
  details.crit = crit;
  end
  
  details.angles = a;
  details.gmfit = bestfit;

      details.xy = xyrotate(xy(:,1),xy(:,2),theta);
      details.dip = HartigansDipTest(sort(details.xy(:,1)));


function tag = GetFigureTag(src)
    a = src;
    while ~isfigure(a) && a ~= 0
        a = get(a,'Parent');
    end
    tag = get(a,'Tag');
    
        
         
function XYCluster(src,b, type, varargin)
DATA = GetDataFromFig(src);
EClusters = getappdata(DATA.toplevel,'Clusters');
useaxdata = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'useaxdata',6)
        useaxdata = 1;
    end
    j = j+1;
end

e = DATA.currentpoint(1);
p = DATA.currentpoint(2);
if useaxdata
    axdata = get(gca,'UserData');
    if isfield(axdata,'probe');
        p = axdata.probe;
    end
    if isfield(axdata,'eid');
        e = axdata.eid;
    end
end

C = EClusters{e}{p};
ex = C.exptno;
tag = GetFigureTag(src);

    onoff = {'on','off'};
        subtype = strmatch(type,{'bestangle' 'quickangle'});
    if sum(strcmp(type,{'save' 'quicksave'}))
        DATA = SaveCluster(DATA, [e p], strcmp(type,'quicksave'));
    elseif strcmp(type,'flip')
        DATA = FlipLineCrit(DATA);
    elseif strcmp(type,'optimize')
        axdata = get(gca,'userdata')
        if isfield(axdata,'eid')
            e = axdata.eid;
            p = axdata.probe;
        end
        FixBoundary(DATA, EClusters, e, p);
    elseif strcmp(type,'twoplot')
        DATA.ellmousept.shape = -1;
        DATA.allclustering = 0;
        set(DATA.toplevel,'UserData',DATA);
        hold off;
        PlotClusterXY(DATA,C,'cellid',[1 2]);
        FastAxes(gca);
    elseif strcmp(type,'nocut')
        DATA.elmousept.shape = -1;
        DATA.allclustering = 0;
        set(DATA.toplevel,'UserData',DATA);
    elseif strmatch(type,{'bestangle' 'quickangle'})
        if subtype == 2
            [a,b,c] = BestAngleGM(C.xy, [], [],'quick');
        else
            [a,b,c] = BestAngleGM(C.xy, [], []);
        end
        C.gmfit1d = c.gmfit;
        C.crit = c.crit;
        SetFigure(DATA,DATA.tag.xyplot);
        hold on;
        if abs(cos(a)) > 0
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        k = (yl + c.crit .* sin(a))./cos(a);
        cx =  c.crit *cos(a) + k * sin(a);
        dx = (diff(yl) .* tan(a))/2;
        plot(cx,yl);
        end
        fprintf('Best 1d %.2f\n',b);
        if DATA.plothist
            SetFigure(DATA,DATA.tag.hist)
            PlotHist(c.xy,C);
            title(sprintf('Best Angle %.1f deg mahal %.2f',a *180/pi,b));
            DATA.tag.misc = 'Misc';
            SetFigure(DATA,DATA.tag.misc);
            hold off;
            plot(c.angles,c.d);
            if subtype == 1
            hold on;
            plot(c.angles,abs(c.dprime),'r');
            end
        end
    elseif strncmp(type,'testline',5) 
        tic;
        line(get(gca,'xlim'),get(gca,'ylim'));
        drawnow('update');
        toc;
    elseif strncmp(type,'clear',5) 
        F = get(get(src,'parent'),'parent');
        if ~isfigure(F) %if called from submenu, go up one more
            F = get(F,'parent');
        end
        cnum = sscanf(type,'clear%d');
        id = find(EClusters{e}{p}.clst == cnum+1);
        EClusters{e}{p}.clst(id) = 1;
        EClusters{e}{p}.next{cnum-1} = [];
        DATA.NewCut.saved = -1;
        setappdata(DATA.toplevel,'Clusters',EClusters);
        set(DATA.toplevel,'UserData',DATA);
        PlotClusterXY(DATA,EClusters{e}{p});
    elseif strncmp(type,'select',6)
        DATA.currentcluster = sscanf(type(7:end),'%d');
        PlotClusterXY(DATA,C,'cellid',DATA.currentcluster);
    elseif strncmp(type,'ellipse',7) || strcmp(type,'lines')
        F = get(get(src,'parent'),'parent');
        if ~isfigure(F) %if called from submenu, go up one more
            F = get(F,'parent');
        end
        if strcmp(tag,DATA.tag.allxy)
            DATA = SetAllXYClustering(DATA,1);
        end
        DATA.elmousept.h= -1;
        DATA.elmousept.down = 0;
        DATA.elmousept.done = 0;
        DATA.elmousept.steps = 0;
        DATA.elmousept.angle = 0;
        DATA.elmousept.cluster = 1;
        shapes = [0 1 0];
        DATA.elmousept.shape = shapes(strmatch(type(1:5),{'ellipses' 'lines'}));
        if DATA.elmousept.shape == 0
            DATA.elmousept.cluster = sscanf(type(8:end),'%d');
            DATA.elmousept.plotargs = {'color' DATA.colors{DATA.elmousept.cluster+1}};
        else
            DATA.elmousept.plotargs = {'color' 'r'};
        end
        DATA.elmousept.color = [1 0 0];
        DATA.currentcluster = DATA.elmousept.cluster;
        DATA.currentcutcluster = DATA.elmousept.cluster;
        if DATA.currentcutcluster > 1 && length(C.next) >= DATA.currentcutcluster-1
            if ~isfield(C.next{DATA.currentcutcluster-1},'xy')
%  Count available spaces. NB Counting over C.next, so there is one more in
%  C
                for j = 1:length(C.next)
                    if isfield(C.next{j},'xy')
                        gotxy(j) = 1;
                    else
                        gotxy(j) = 0;
                    end
                end
                gid = find(gotxy)+1;
                switch sum(gotxy)
                    case 1, 
                        str = questdlg(sprintf('No XY for Cluster %d: Choose Space',DATA.currentcutcluster),'Select Space','Cluster1',sprintf('Cluster%d',gid(1)),'Cluster1');
                        usespace = sscanf(str,'Cluster%d')-1;
                    case 2, 
                        str = questdlg(sprintf('No XY for Cluster %d: Choose Space',DATA.currentcutcluster),'Select Space','Cluster1',sprintf('Cluster%d',gid(1)),sprintf('Cluster%d',gid(2)),'Cluster1');
                        usespace = sscanf(str,'Cluster%d')-1;
                    otherwise,
                    usespace = 0;
                    errordlg(sprintf('No XY for Cluster %d: Using 1',DATA.currentcutcluster));
                end
                figure(F);
                if usespace == 0
                    XY = C.xy;
                else
                    XY = C.next{usespace}.xy;
                end
                C.next{DATA.currentcutcluster-1}.xy = XY;
                C.next{DATA.currentcutcluster-1}.exptid = C.exptid;
                DATA.elmousept.xyspace = usespace;
            else
                XY = C.next{DATA.currentcutcluster-1}.xy;
            end
            if isfield(C.next{DATA.currentcutcluster-1},'pos');
                DATA.elmousept.pos = C.next{DATA.currentcutcluster-1}.pos;
            else
                DATA.elmousept.pos = C.pos;
            end
            if ~isfield(C.next{DATA.currentcutcluster-1},'xyr')
                DATA.elmousept.xyr = C.xyr;
                DATA.elmousept.angle = C.angle;
            else
                DATA.elmousept.xyr = C.next{DATA.currentcutcluster-1}.xyr;
                DATA.elmousept.angle = C.next{DATA.currentcutcluster-1}.angle;
            end
        else
            XY = C.xy;
            DATA.elmousept.angle = C.angle;
            if isfield(C,'pos')
                DATA.elmousept.pos = C.pos;
                DATA.elmousept.xyr = C.xyr;
            end
        end
        if length(DATA.elmousept.xyr) < 3
            DATA.elmousept.xyr(3) = 0;
        end
        hold off;
        PlotClusterXY(DATA,C,'cellid',DATA.currentcutcluster);
        AddCellLabels(DATA, DATA.currentpoint(1), DATA.currentpoint(2),'NW');
        FastAxes(gca);
%        DATA.elmousept.dragfcn = get(F,'WindowButtonMotionFcn');
        %should get old Fcns, then reset them after button release
        set(F, 'WindowButtonDownFcn',@XYButtonPressed);
        set(F, 'WindowButtonMotionFcn',@XYButtonDragged);
        set(F, 'WindowButtonUpFcn',@XYButtonReleased);
        set(F, 'WindowScrollWheelFcn',@ScrollWheel);
        set(F,'UserData',DATA.toplevel);
        set(DATA.toplevel,'UserData',DATA);
   
    elseif strncmp(type,'replot',6)
        hold off;
        if strcmp(type,'replotgm')
            PlotClusterXY(DATA,C,'plotgmcid');
        else
            PlotClusterXY(DATA,C);
        end
        DATA.NewCut.exptid = 0;
        DATA.NewCut.probe = 0;
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmp(type,'revert',6)
        e = DATA.currentpoint(1);
        p = DATA.currentpoint(2);
        X = ReloadClusters(DATA, e);
        EClusters{e}{p} = X{p};
        setappdata(DATA.toplevel,'Clusters',EClusters);
        DATA.NewCut.exptid = 0;
        DATA.NewCut.probe = 0;
        SetFigure(DATA,DATA.tag.xyplot); hold off;
        PlotClusterXY(DATA,X{p});
        set(DATA.toplevel,'UserData',DATA);
    end

function DATA = SaveCluster(DATA, pt, quick, varargin)
    e = pt(1);
    p = pt(2);
    EClusters = getappdata(DATA.toplevel,'Clusters');

    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    C = EClusters{e}{p};
    ex = C.exptno;
%if saved == -1, means a cluster was deleted
    if DATA.NewCut.saved >= 0 && DATA.NewCut.probe > 0
        C = DATA.NewCut;
    else
        mycprintf('red','Saving - but no new cluster\n');
    end
    if rem(ex,1) > 0.001
        xs = 'a';
    else
        xs = '';
    end
    cfile = [DATA.name '/Expt'  num2str(floor(ex)) xs 'ClusterTimes.mat'];
    afile = [DATA.name '/Expt'  num2str(floor(ex)) xs 'AutoClusterTimes.mat'];
    if exist(cfile,'file')
        load(cfile);
    elseif exist(afile,'file')
        load(afile);
        for j = 1:length(Clusters)
            Clusters{j}.auto = 1;
        end
    end
    dfile = [DATA.name '/Expt'  num2str(floor(ex)) xs 'ClusterTimesDetails.mat'];
    if quick == 0
        if exist(dfile)
            load(dfile);
        else
            ClusterDetails = {};
        end
        afile = [DATA.name '/Expt'  num2str(floor(ex)) xs 'AutoClusterTimesDetails.mat'];
        if exist(afile)
            a = load(afile);
            for j = 1:length(a.ClusterDetails)
                if j > length(ClusterDetails) || isempty(ClusterDetails{j}) || ~isfield(ClusterDetails{j},'t')
                    ClusterDetails{j} = a.ClusterDetails{j};
                end
            end
        end
        ClusterDetails{p}.clst = C.clst;
    end
    %C.xy should no longer be changed by plotclusters, so don't save it out
    %Need modify disk file (Clusters) and the version in memory now (EClusters)
    %Because PlotClusters adds fields that we don't want save out, easiest
    %thing is just to make changes separatley, not copy one to anaother
    if C.cluster > 1
        n = C.cluster-1;
        if ~isfield(C.next{n},'space') %new cluster - cut in same space as 1
            C.next{n}.space = C.space;
        end
        if ~isfield(C.next{n},'xy')
            C.next{n}.xy = C.xy;
        end

        C.next{n}.ctime = now;
        C.next{n}.savetime = now;
        if n > length(EClusters{e}{p}.next)
            EClusters{e}{p}.next{n} = {};
        end
        EClusters{e}{p}.next{n} = CopyFields(EClusters{e}{p}.next{n},C.next{n},...
                'shape','crit','xyr','sign','angle','savetime','times','space','ctime','savetime','clst');
        Clusters{p}.next{n} = CopyFields(EClusters{e}{p}.next{n},C.next{n},...
            'shape','crit','xyr','sign','angle','savetime','times','space','ctime','savetime');
        %C.times is from ClusterDetails.t
        Clusters{p}.next{n}.times = C.times(C.clst == C.cluster+1);
        Clusters{p}.next{n}.clusterprog = sprintf('PlotClusters V %.2f',DATA.version);
        Clusters{p}.next{n}.manual = 2;
        Clusters{p}.next{n}.pos = C.next{n}.pos;
        Clusters{p}.next{n}.aspectratio = C.next{n}.aspectratio;

        if quick == 0
            if n > length(ClusterDetails{p}.next)
                ClusterDetails{p}.next{n} = {};
            end
            ClusterDetails{p}.next{n}.xy = C.xy;
            if C.shape > 0
                ClusterDetails{p}.next{n}.xy = xyrotate(C.next{n}.xy(:,1),C.next{n}.xy(:,2),C.next{n}.angle);
            end
        end
    else
        ClusterDetails{p}.xy = xyrotate(C.xy(:,1),C.xy(:,2),C.angle);
        Clusters{p}.shape = C.shape;
        Clusters{p}.crit = C.crit;
        Clusters{p}.xyr = C.xyr;
        Clusters{p}.angle = C.angle;
        Clusters{p}.sign = C.sign;
        Clusters{p}.ctime = C.ctime;
        if isfield(C,'t')
            Clusters{p}.times = C.t(C.clst == 2);
            C.times = C.t;
        else
            Clusters{p}.times = C.times(C.clst == 2);
        end
        Clusters{p}.savetime(1) = now;
        EClusters{e}{p} = CopyFields(EClusters{e}{p},C,...
            'shape','crit','xyr','sign','angle','savetime','times','clst');
        Clusters{p}.clusterprog = sprintf('PlotClusters V %.2f',DATA.version);
        Clusters{p}.aspectratio = C.aspectratio;
        Clusters{p}.manual = 2;
    end

    setappdata(DATA.toplevel,'Clusters',EClusters);
    testing = 0;
    if testing
        ofile = strrep(cfile,'ClusterTimes','NewClusterTimes');
    else
        ofile = cfile;
    end
    if isempty(Clusters{p})
        Clusters{p} = rmfield(C,{'xy' 'clst' 'Evec'});
        for j = 1:length(Clusters{p}.next)
            Clusters{p}.next{j} = rmfields(Clusters{p}.next{j},'xy');
        end
    end
    fprintf('Saving %s\n',ofile);
    if exist('FullVData')
        save(ofile,'Clusters','FullVData');
    else
        fprintf('No FullVData for %s\n',ofile);
        save(ofile,'Clusters');
    end
    if testing
        ofile = strrep(dfile,'ClusterTimes','NewClusterTimes');
    else
        ofile = dfile;
    end
    if quick == 0
        if length(ClusterDetails{p}.clst) == length(ClusterDetails{p}.t)
            save(ofile,'ClusterDetails');
        else
            fprintf('E%d P%d Can''t Save ClusterDetails: size mismatch',C.exptid,C.probe(1));
            errordlg(sprintf('Can''t Save ClusterDetails: size mismatch','ClusterError','modal'));
        end
    end
    fprintf('Saved Cluster for probe %d Expt %d\n',p,ex);
    DATA.NewCut.saved = 1;


     function in = InGraph(pt, ax)
        xl = get(ax,'Xlim');
        yl = get(ax,'Ylim');
      in = pt(1,1) > xl(1) & pt(1,1) < xl(2) & pt(1,2) > yl(1) & pt(1,2) < yl(2);

function xy = AxPos(ax, pos)
    xl = get(ax,'Xlim');
    yl = get(ax,'Ylim');
    xy(1) = xl(1) + diff(xl) .* pos(1);         
    xy(2) = yl(1) + diff(yl) .* pos(2);         
    
function XYButtonPressed(src, data)
    DATA = GetDataFromFig(src);
    tag = GetFigureTag(src);
    testnew = 0;
    if strcmp(tag,DATA.tag.allxy) && DATA.allclustering == 0
 % need to set the currentpoint so that plots right thing.
        axdata = get(gca,'UserData');
        DATA.currentpoint(1) = axdata.eid;
        DATA.currentpoint(2) = axdata.probe;
        set(DATA.toplevel,'UserData',DATA);
        return;
    end
    if DATA.elmousept.shape < 0
        return;
    end
DATA.ts = now;
start = get(gca,'CurrentPoint');
if InGraph(start,gca)
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})
    if DATA.currentcluster ~= DATA.currentcutcluster
        Clusters = getappdata(DATA.toplevel,'Clusters');
        fprintf('Replotting XY for cluster %d\n',DATA.currentcutcluster);        
        PlotClusterXY(DATA,Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)},'cellid',DATA.currentcutcluster);
        AddCellLabels(DATA, DATA.currentpoint(1), DATA.currentpoint(2),'NW');
        xy = AxPos(gca,[1 1]);
        text(xy(1),xy(2),sprintf('Reset to cl%d',DATA.currentcutcluster),'HorizontalAlignment','right','verticalalignment','top');
        DATA.currentcluster = DATA.currentcutcluster;
    end
    DATA.elmousept.mode = mode;
    distance = DistanceToEllipse(DATA.elmousept,start(1,1:2));
    DATA.elmousept.steps = 0;
    if ishandle(DATA.elmousept.h)
        delete(DATA.elmousept.h);
    end
    set(gca,'xlimmode','manual','ylimmode','manual');
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    DATA.elmousept.aspectratio = diff(yl)./diff(xl);

    if  mode  == 2 %right button press
        DATA.elmousept.down = 3;
        DATA.elmousept.done = 0;
        if DATA.elmousept.shape == 0 %for line, don't move
            DATA.elmousept.start = start(1,1:2);
        end
    elseif distance < 1.05 %test fails for NaN
%pressed inside ellipse, so just tmove it. 
        DATA.elmousept.down = 2;
        DATA.elmousept.pos =[0  0 0 0 ];
        DATA.elmousept.start = start(1,1:2);
        DATA.elmousept.axis = gca;
%set up pos in case released with no drag        
    DATA.elmousept.pos(1) = DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    else
        if testnew 
        ax = gca;
        axdata = get(ax,'UserData');
        if isfield(axdata,'allxy') && axdata.allxy ==1
            Clusters = getappdata(DATA.toplevel,'Clusters');
            GetFigure(DATA.tag.xyplot);
            PlotClusterXY(DATA,Clusters{axdata.eid}{axdata.probe},'shorttitle','tight');
            DATA.xyax = gca;
            axes(ax);
            DATA.allxyax = ax;
        end
        end
        DATA.elmousept.down = 1;
        DATA.elmousept.done = 0;
        if ~isfield(DATA.elmousept,'angle') %if there, should have been set to match cluster
            DATA.elmousept.angle = 0;
        end
        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
        DATA.elmousept.axis = gca;
    end
set(DATA.toplevel,'UserData',DATA);
end


function distance = DistanceToEllipse(E, pos);
   
if isempty(E) | ~isfield(E,'pos');
    distance = NaN;
    return;
end


a(1) = (E.pos(3)+E.pos(1))/2; %x radius
a(2) = (E.pos(4)+E.pos(2))/2;


if E.shape == 1
angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
r(1) = abs(E.pos(3)-E.pos(1) + i * (E.pos(4)-E.pos(2)))/4;
r2 = r(1)/10;
else
r(1) = (E.pos(3)-E.pos(1))/2; %x radius
r(2) = (E.pos(4)-E.pos(2))/2;
angle = E.angle;
end
xy = pos - a;
xy = xy ./r;
cn = cos(-angle);
sn = sin(-angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = sum(p.^2);


function XYButtonReleased(src, data)
DATA = GetDataFromFig(src);
if DATA.elmousept.down == 0 
    return;
end
mode = DATA.elmousept.down;
DATA.elmousept.mode = mode;
start = get(gca,'CurrentPoint');
DATA.elmousept.done = 1;
p = DATA.elmousept.pos;
DATA.elmousept.down = 0;
axdata = get(gca,'UserData');
if isfield(axdata,'probe')
    pi = axdata.probe;
    ei = axdata.eid;
    DATA.currentpoint = [ei pi];
    fprintf('Clustering E%dP%d\n',ei,pi);
else
    pi = DATA.currentpoint(2);
    ei = DATA.currentpoint(1);
end

if mode == 1
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
elseif mode == 2  %button was pressed inside ellipse, just move it
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
end
xyr = DATA.elmousept.xyr;
set(DATA.toplevel,'UserData',DATA);
%touch inside ellispe to make cut. If drawing a line, make cut when
%released
if mode == 2  || (DATA.elmousept.shape  ==1 && mode == 1)
%PCCluster(DATA,DATA.elmousept,1);
       C =  ClassifySpikes(DATA,DATA.elmousept,1);
        if (DATA.plothist | DATA.refitgm) && DATA.elmousept.shape == 1
            [a,b]  = GMDip(C.xy(:,1),0);
            [c,d,e] = GMfit(C.xy,2,1,'idlist',C.clst);
            title(sprintf('E%dP%d mahal %.2f,%.2f(1/2)\n',ei,pi,b.mahal(b.best),d));
            if DATA.plothist
                SetFigure(DATA,DATA.tag.hist,'front');
                PlotHist(C.xy, C);
                title(sprintf('E%dP%d mahal %.2f\n',ei,pi,b.mahal(b.best)));
            end
        elseif DATA.refitgm && DATA.elmousept.shape == 0
            [a,b,c] = GMfit(C.xy,2,1,'idlist',C.clst);
            tic;
            [d,e, details] = BestAngleGM(C.xy,a,[],'quick');
            toc
            [aa,bb]  = GMDip(details.xy(:,1),0);
            title(sprintf('E%dP%d mahal %.2f (%.2f 1d)\n',ei,pi,b,max(bb.mahal)));
        end
end
figure(src);


function C = ClassifySpikes(DATA, E, mode)
    C = [];
    if E.shape < 0 && ~strcmp(mode,'flip')
        cprintf('blue','Can''t Classify Spikes until set cluster\n');
        return;
    end
    Clusters = getappdata(DATA.toplevel,'Clusters');
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if isfield(DATA,'NewCut') && DATA.NewCut.probe == p && DATA.NewCut.exptid == e
        C = DATA.NewCut;
    else
        C = Clusters{e}{p};
    end
    C.exptid = e;
    if E.shape < 0 || strcmp(mode,'flip') %no mouse input - just flipping
        E = C;
        E.mode = 3;
    end
    C.cluster = E.cluster;
    if E.cluster > 1
        nx = E.cluster-1;
%set something in next{nx} to be sure it exists
        C.next{nx}.aspectratio = E.aspectratio;
        if isfield(C.next{nx},'xy')
            xy = C.next{E.cluster-1}.xy;
        elseif isfield(E,'xyspace') && E.xyspace > 0
            C.next{nx} = C.next{E.xyspace};
            C.next{nx}.cluster = nx+1;
            xy = C.next{E.xyspace}.xy;
        else
            C.next{nx} = CopyFields(C.next{nx},C,{'space' 'TemplateUsed'});
            xy = C.xy;
        end
        C.next{nx}.aspectratio = E.aspectratio;
        C.next{nx}.shape = E.shape;
        C.next{nx}.ctime = now;
    else
        xy = C.xy;
        C.shape = E.shape;
        C.aspectratio = E.aspectratio;
        C.ctime = now;
    end
    
    if ~isfield(DATA,'Expt')
        Expts = getappdata(DATA.toplevel,'Expts');
        DATA.Expt = Expts{e};
    end
    if E.shape == 0 % ellipse
        cx = E.xyr(1);
        cy = E.xyr(2);
        rx = E.xyr(3);
        ry = E.xyr(4);
        if isfield(E,'aspectratio') & E.aspectratio > 0
            xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./E.aspectratio,E.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./E.aspectratio)).^2;
        else
            xys = xyrotate(xy(:,1)-cx,xy(:,2)-cy,E.angle);
        r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
        end
        id = find(r < 1);
        nspk = length(id);
        C.clst(id) = E.cluster+1;
        if diff(size(C.clst) >0)
            C.clst = C.lst';
        end
        id = find(r >=1 & C.clst == E.cluster+1);
        C.clst(id) = 1;
        PlotClusterPoints(C,[],E.cluster,'xydata',xy);
        FinishXYPlot(gca, DATA, C.exptid,C.probe(1));
        if E.cluster > 1
            n = E.cluster-1;
            C.next{n}.xyr = E.xyr;
            C.next{n}.angle = E.angle;
%            C.next{n}.crit = 0;
            C.next{n}.pos = E.pos;
            DrawEllipse(E,'color',DATA.colors{n+2});
        else
            C.xyr = E.xyr;
            C.angle = E.angle;
            DrawEllipse(E,'color','r');
        end
        if ~isfield(C,'crit')
            C.crit = 1;
        end
    elseif E.shape == 1
        if E.mode == 3 || strcmp(mode,'flip') %invert sign
            C.sign = C.sign * -1;
        end
        E.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        xy = xyrotate(C.xy(:,1),C.xy(:,2),E.angle);
        crit = xyrotate(E.pos([1 3]),E.pos([2 4]),E.angle);
        if C.sign < 0
        id = find(xy(:,1) < mean(crit(:,1)));
        nid = find(xy(:,1) > mean(crit(:,1)));
        else
        id = find(xy(:,1) > mean(crit(:,1)));
        nid = find(xy(:,1) < mean(crit(:,1)));
        end
        C.angle = E.angle;
        C.clst(id) = 2;
        oldi = find(C.clst(nid) ==2);
        C.clst(nid(oldi)) = 1;
        nspk = length(id);
        C.crit = mean(crit(:,1));
%        C.xy = xy;
        Clusters{e}{p}.crit = C.crit;
%        E.pos(1) = mean(crit(:,1));
 %       E.pos(3) = mean(crit(:,1));
 %       E.pos(4) = max(xy(:,2));
 %       E.pos(2) = min(xy(:,2));
        E.pos(5) = max(xy(:,1));
        E.pos(7) = min(xy(:,1));
        E.pos(6) = mean(E.pos([4 2]));
        E.pos(8) = E.pos(6);
%        C.angle = 0;
        if ~strcmp(mode,'flip')
            PlotClusterPoints(C,[],E.cluster,'xydata',C.xy);
            FinishXYPlot(gca, DATA, C.exptid,C.probe(1));
        fprintf('%.2f ',E.pos);
        fprintf('\n');
        DrawEllipse(E);
        C.pos = E.pos;
        end
        DATA.elmousept.pos = E.pos;
        if C.space(1) == 6
        end
    end
    C.id = find(C.clst ==E.cluster+1);
    muid = find(C.clst == 1);
    Clusters{e}{p}.clst = C.clst;
    if DATA.showspkmean
        Spks = CheckAllSpikes(DATA, e,p);
%can only recalc mean for this probe. Spks won't have matching time samples
%on other probes
        if length(C.clst) == size(Spks{e,p}.values,1)
            for j = p
                Vs = Spks{e,p}.Vscale;
                Clusters{e}{p}.MeanSpike.ms(j,:) = mean(Spks{e,j}.values(C.id,:)) .*Vs;
                Clusters{e}{p}.MeanSpike.mu(j,:) = mean(Spks{e,j}.values(muid,:)) .*Vs;
            end
            oldf = gcf;
            SetFigure(DATA, DATA.tag.spkmean);
            h = PlotMeanSpike(Clusters{e}{p},p,E.cluster,DATA.plotmeantype,DATA);
            set(0,'CurrentFigure',oldf);
        else
            fprintf('E%dP%d Spikes(%s)/Cluster(%d) length mismatch',e,p,size(Spks{e,p}.values,1),length(C.clst));
        end
    end
%    Clusters{e}{p}.xy = C.xy; Shouldn't have to change xy any more
if isfield(C,'crit')
    Clusters{e}{p}.crit = C.crit;
    Clusters{e}{p}.sign = C.sign;
else
    Clusters{e}{p}.crit = 0;
end
    setappdata(DATA.toplevel,'Clusters',Clusters);
    DATA.NewCut = C;
    DATA.NewCut.saved= 0;
    DATA.NewCut.probe = p; %temporay fix until grid clusters get correct numbers
    if isfield(C,'t')
        if legnth(C.t) > length(C.times)
            fprintf('!!!!E%dC%d  has t > times\n',C.exptno,C.probe(1));
            DATA.NewCut.times = C.t;
        end
    end
    dur = sum([DATA.Expt.Trials.dur])./10000;
    s = sprintf('%d/%d Spikes (%.1f/%.1f Hz)',nspk,length(C.clst),nspk./dur,length(C.clst)./dur);
    title(s);
    fprintf('%s\n',s);
%need to record new cluster, and new E.pos
    SetFigure(DATA, DATA.tag.expt);
%    PlotSelectedExpts(DATA);
    set(DATA.toplevel,'UserData',DATA);
    if DATA.plotspks
        SetFigure(DATA, DATA.tag.spikes);
        AllSpikes = CheckAllSpikes(DATA, e, p);
        h = QuickSpikes(DATA, AllSpikes{e,p},C);
        AddLineContextMenu(DATA, h, e, p);
    end
    if DATA.plothist
        SetFigure(DATA, DATA.tag.hist);
        PlotClusterHistogram(DATA, C, 1,'cluster',DATA.currentcutcluster);
    end
    
function ScrollWheel(src, evnt)
DATA = GetDataFromFig(src);
DATA.elmousept.angle = DATA.elmousept.angle+0.02*evnt.VerticalScrollCount;
fprintf('Angle %.2f\n',DATA.elmousept.angle);
DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
set(DATA.toplevel,'UserData',DATA);

function FinishXYPlot(ax, DATA, e,p)
    axis('tight');
    set(ax,'xtick',[],'ytick',[]); 
    tag = GetFigureTag(ax);
    if strcmp(tag,DATA.tag.allxy)
        set(ax,'ButtonDownFcn',{@HitXYPlot, e,p});
    end
    
function XYButtonDragged(src, data)
    persistent tslast;

    ts = now; 
DATA = GetDataFromFig(src);
if isfield(DATA,'elmousept') &&  DATA.elmousept.down >0
%    fprintf('D%d,%d %.2f %.2f\n',DATA.elmousept.down,DATA.elmousept.steps,DATA.elmousept.pos(1),DATA.elmousept.pos(3));
control_is_down = getappdata(0,'control_is_down');
cc = get(gcf,'currentcharacter');
 DATA.elmousept.steps = DATA.elmousept.steps +1;

 if  DATA.elmousept.down == 1
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
%    axes(DATA.xyax);
    if ~control_is_down %%cntrl stops the ellipse from drawing. Speeds things up in and AllXY plot
        DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
    end
%    axes(DATA.allxyax);
%    DATA.elmousept.lh = testDrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
% with drawnow here, takes a long time with AllXY. But without
% jst takes that time outside the subroutine
% try using Line (annotaiton) and deleting each time? 
%    drawnow;
%    mytoc(DATA.ts);
elseif  DATA.elmousept.down == 3 %set radius with R mouse button
    if DATA.elmousept.shape == 0 %set radius with R mouse button
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
    else %for line R press doesnt move line. Just inverts cluster.sign
    end
elseif  DATA.elmousept.down == 2 %moving ellipse
    start = get(gca,'CurrentPoint');
    if DATA.elmousept.steps > 5
        start = get(gca,'CurrentPoint');
    end
    DATA.elmousept.pos(1) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)-DATA.elmousept.xyr(3);
    DATA.elmousept.pos(2) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)-DATA.elmousept.xyr(4);
    DATA.elmousept.pos(3) = start(1,1)-DATA.elmousept.start(1)+DATA.elmousept.xyr(1)+DATA.elmousept.xyr(3);
    DATA.elmousept.pos(4) = start(1,2)-DATA.elmousept.start(2)+DATA.elmousept.xyr(2)+DATA.elmousept.xyr(4);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});

end
set(DATA.toplevel,'UserData',DATA);
if DATA.profiling
    if tslast == 0
        tslast = ts;
    end
    fprintf('Took %.3f (%.3f)\n',mytoc(ts),(ts-tslast)*60*60*24);
end
end
tslast = ts;


function CellTrackMenu(a,b, type)
    DATA = GetDataFromFig(a);
    
    if strcmp(type,'setmenu')
    fcn = get(a,'value') -3;    
        if fcn == 5
            PlotAllCellSpikes(DATA);
        elseif fcn == 6
            PlotAllCellXY(DATA);
        elseif fcn == 7
            PlotAllCellMean(DATA,'lineonly');
        elseif fcn == 8
            PlotAllCellMean(DATA,'imageonly');
        elseif fcn == 9
            ExcludeTrialsForCell(DATA, 0, 1);
        elseif fcn == 10
            ExcludeTrialsForCell(DATA, 0, 2);
        end
    elseif strcmp(type,'setplane')
        fcn = get(a,'value');
        DATA.clusteroffset = (fcn-1)*2;
        PlotCellList(DATA, 'showfig');
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmp(type,'setcluster')
        e = DATA.currentpoint(1);
        p = DATA.currentpoint(2);
        cl = get(a,'value');
        DATA.cellcluster = cl;
        if size(DATA.CellList,3) >= cl
        cell = DATA.CellList(e,p,cl);
        else
            cell = 0;
        end
        if cell > 0 && 0
            it = findobj('Tag','CellNumberId');
            set(it,'value',cell);
        end
        set(DATA.toplevel,'UserData',DATA);
    end

function DATA = ExcludeTrialsForCell(DATA, probe, cluster, varargin)
    
    useall = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'reset',5)
            useall = 1;
        end
       j = j+1;
    end

    Clusters = getappdata(DATA.toplevel,'Clusters');
    e = DATA.currentpoint(1);
    C = Clusters{e}{probe};
    if probe == 0
        it = findobj('Tag','CellNumberId');
        cellid = get(it,'value');
    else
        cellid = DATA.CellList(e, probe, cluster);
    end
%    if cellid == 0
%        return;
%    end
       
    F = findobj('type','figure','Tag',DATA.tag.spikes);
    if length(F) == 1
        it = findobj(F,'Tag','ChooseTrial');
    else
        it = findobj(DATA.spoolfig,'Tag','TrialList');
    end
    trials = get(it,'value');
    oldt = FindExcludedTrials(DATA, e, probe, cluster,C);
    trials = unique([trials oldt]);
    if useall
        DATA.CellDetails.excludetrials{ e, probe, cluster} = [];
    else
        DATA.CellDetails.excludetrials{e, probe, cluster} = DATA.trialids{e}(trials);
    end
    fprintf('Cell%d (E%dP%dCluster%d) Excluding: %s\n',cellid,e,probe,cluster,sprintf(' %d',trials));
set(DATA.toplevel,'UserData',DATA);
SaveCellList(DATA);
if length(F) == 1
    SetTrialList(DATA,C,DATA.currenttrial);
else
    SpoolSpikes(DATA.spoolfig,'excludelist', trials);
end

function [e,p] = cPoint(DATA)
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    
function SetTrialList(DATA, C, strial)
    F = findobj('type','figure','Tag',DATA.tag.spikes);
    if nargin < 2 || isempty(C)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    C = Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)};
    end
    tn = [DATA.Expt.Trials.Trial]';
    xcl = FindExcludedTrials(DATA, DATA.currentpoint(1),DATA.currentpoint(2), 1, C);
    tn(xcl) = tn(xcl).*-1;
    strial(strial > length(tn)) = 1;
    
    if length(F) == 1
        it = findobj(F,'Tag','ChooseTrial');
        set(it,'string',num2str(tn),'value',strial);
    end
        
function PlotAllCellMean(DATA, type)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);
    [eid, id] = sort(eid);
    cid = cid(id);
    clid = clid(id);
    [nr,nc] = Nsubplots(length(eid));
    SetFigure(DATA,DATA.tag.allmean);
    subplot(1,1,1);
    for j = 1:length(eid)
        mysubplot(nr,nc,j);
        hold off; 
        PlotMeanSpike(Clusters{eid(j)}{cid(j)},cid(j),clid(j),type,DATA);
        set(gca,'Xtick',[],'Ytick',[]);
        axdata.toplevel = DATA.toplevel;
        axdata.eid = eid(j);
        axdata.probe = cid(j);
        set(gca,'UserData',axdata,'uicontextmenu',AddContextMenu(gca,'subplot'));
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        c = get(gca,'Children');
        for k = 1:length(c)
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        end
        d = DistanceMeasure(Clusters{eid(j)}{cid(j)},1,DATA.mahaltype);
        h = text(xl(2),yl(1),sprintf('E%dP%d/%d:%.2f %s',eid(j),cid(j),clid(j),d,DATA.expnames{eid(j)}),'horizontalalignment','right','verticalalignment','bottom');
        if d < DATA.crit.mahal
            set(h,'color','r');
        end
        if j == 1
            title(sprintf('Cell%d',cellid));
        end
    end

function DATA = SetAllXYClustering(DATA, onoff) 
    F = SetFigure(DATA,DATA.tag.allxy);
    if onoff
    set(F, 'WindowButtonDownFcn',@XYButtonPressed);
    set(F, 'WindowButtonMotionFcn',@XYButtonDragged);
    set(F, 'WindowButtonUpFcn',@XYButtonReleased);
    set(F, 'WindowScrollWheelFcn',@ScrollWheel);
    set(F,'UserData',DATA.toplevel);
    set(F,'keypressfcn',@XYKeyPressed);
    else
        c = get(F,'children');
    end
    DATA.allclustering = onoff;
    
 function MarkTrialStarts(Expt, ticks, xcl)
     h = ishold;
     hold on;
     yl = get(gca,'ylim');
     yl(2) = yl(1) + diff(yl)/8;
     if isfield(Expt.Header,'timeoffset')
         toff = Expt.Header.timeoffset;
     else
         toff = 0;
     end
     for j = 1:length(Expt.Trials)
         if ticks
             t = Expt.Trials(j).Start(1);
         else
             t = Expt.Trials(j).Start(1)./10000;
         end
         t = t+toff;
         if ismember(j,xcl)
             plot([t t],yl,'m-');
         else
             plot([t t],yl,'k-');
         end
     end
    if h == 0
        hold off;
    end
    
function sdx = PlotXYSequence(DATA, pos, varargin)
    plottype = 'time';
    T = DATA.Expt.Trials;
    Expt = DATA.Expt;
    setfig = 1;

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'density',6)
            plottype = varargin{j};
        elseif strncmpi(varargin{j},'expt',4)
            j = j+1;
            Expt = varargin{j};
            T = Expt.Trials;
        elseif strncmpi(varargin{j},'noplot',6)
            plottype = 'noplot';
        end
        j = j+1;
    end
    
    if isstruct(pos) %passed  Cluster
        C = pos;
        e = C.exptid;
        p = C.probe;
        setfig = 0;
    else
        Clusters = getappdata(DATA.toplevel,'Clusters');
        e = pos(1);
        p = pos(2);
        C = Clusters{e}{p};
    end
    nc = unique(C.clst(:))';
    nt = length(T);
    if ~strcmp(plottype,'noplot')  && setfig
    SetFigure(DATA,DATA.tag.xyseq);
    hold off;
    end
%    need to use time or trial bins for this test, as rate changes 
    if isfield(Expt.Header,'timeoffset')
        toff = Expt.Header.timeoffset;
    else
        toff = 0;
    end
    allx = [];
    ally = [];
    for j = nc
        xcl = FindExcludedTrials(DATA,e,p,j-1, C);
        id = find(C.clst == j);
        if C.shape == 0
            r = CalcRadius(C, C.xy);
        else
            xy = xyrotate(C.xy(:,1),C.xy(:,2),C.angle);
            r = xy(:,1);
        end
        if strcmp(plottype,'time')
            plot(C.times(id)+toff,r(id),'.','color',DATA.colors{j},'buttondownfcn',@HitXYseq);
            hold on;
        end
        smw = round(nt/10);
        if smw > 10
            smw = 10;
        end
        for k = 1:nt-smw
            ts(k) = T(k).Start(1)./10000;
            if sum(ismember([k:k+smw],xcl) == 0)
            if diff(size(C.clst)) > 0
                id  = find(C.clst == j & C.times.*10000 > T(k).Start(1) & C.times.*10000 < T(k+smw).End(end));
            else
                id  = find(C.clst' == j & C.times.*10000 > T(k).Start(1) & C.times.*10000 < T(k+smw).End(end));
            end
            sds(k) = std(r(id));
            sms(k) = mean(r(id));
            else
                sds(k) = NaN;
                sms(k) = NaN;
            end
        end
        sid = find(~isnan(sds));
        sds = sds(sid);
        sms = sms(sid);
        ssds(j) = std(sds);
        msds(j) = mean(sds);
        sdx(j) = std(sds)./mean(sds);
        mdx(j) = std(sms)./mean(sms);
    end
if strcmp(plottype,'time')
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    if yl(2) > 0
        plot(ts(sid)+toff,sds .*yl(2)./max(sds),'ms');
    else
        plot(ts(sid)+toff,sds .*yl(1)./max(sds),'ms');
    end
    plot(ts(sid)+toff,sms,'r-');
    if C.shape ~= 0
        plot(minmax(ts)+toff,[C.crit C.crit],'r-');
    end
    MarkTrialStarts(Expt,0,xcl);
    title(sprintf('SDindex %.2f (%.2f/%.2f) CV %.2f',sdx(end),ssds(end),msds(end),mdx(end)));
    set(gca,'buttondownfcn',@HitXYseq);
elseif strcmp(plottype,'density')
    edges(1) = T(1).Start(1);
    for j = 2:length(T);
        edges(j) = (T(j).Start(1)+T(j-1).End(end))/2
    end
    edges(j+1) = T(j).End(end);
    [a,b] = histc(C.times,edges./10000);
    DensityPlot(b,r,'ynormal');
end

function HitXYseq(a,b)
    DATA = GetDataFromFig(a);
    SetTrialList(DATA,{},DATA.currenttrial);
    xy = get(gca,'currentpoint');
    if isfield(DATA.Expt.Header,'timeoffset')
        toff = DATA.Expt.Header.timeoffset;
    else
        toff = 0;
    end
    id = find([DATA.Expt.Trials.TrialStart] < (xy(1,1)-toff) .* 10000);
    if isempty(id)
        id = 1;
    end
    T = DATA.Expt.Trials(id(end));
    fprintf('Trial %d:  %d, id%d\n',id(end),T.Trial,T.id);
    SpoolSpikes(DATA,DATA.currentpoint,'ids',id(end))


function [eid, cid, clid] = FindCell(DATA, cellid, expts)
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);
    [eid, id] = sort(eid);
    cid = cid(id);
    clid = clid(id);
    if nargin > 2 %find cell for particualr expts
        id = find(ismember(eid, expts));
        eid = eid(id);
        cid = cid(id);
        clid = clid(id);
    end

    
function PlotExptCells(DATA, type)

function PlotAllCell(DATA, type)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');
    if isempty(cellid)
        cellid = 1;
    end
    colors = mycolors('spkcolors');
    
    ts = now;
    if DATA.profiling
        profile on;
    end
    
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);

    
    [eid, id] = sort(eid);
    cid = cid(id);
    AllSpikes = CheckAllSpikes(DATA, eid, cid);
    clid = clid(id);
    cid = 1+mod(cid-1,DATA.nprobes);
    [nr,nc] = Nsubplots(length(eid));
    if strncmp(type,'allspks',7)
        SetFigure(DATA,DATA.tag.allspikes);
    else
        SetFigure(DATA,DATA.tag.allxy);
    end
    subplot(1,1,1);
    for j = 1:length(eid)
        mysubplot(nr,nc,j);
        hold off; 
        if strcmp(type,'trighist')
            hist(Clusters{eid(j)}{cid(j)}.trighist);
        elseif strncmp(type,'allspks',7)
            C = Clusters{eid(j)}{cid(j)};
            if strcmp(type,'allspkswithmean')
                h = QuickSpikes(DATA, AllSpikes{eid(j),cid(j)},C,'showmean');
            else
                h = QuickSpikes(DATA, AllSpikes{eid(j),cid(j)},C);
            end
            if clid(j) ~= 1
                
                np = clid(j)+1;
                if length(h) >= np
                    set(h(np),'color','r');
                    set(h(2),'color',colors{np});
                elseif length(h) > 1
                    set(h(2),'color','k');
                else
                    set(h(1),'color','k');
                end
            end
            if DATA.show.linecontextmenus
                AddLineContextMenu(DATA, h, eid(j), cid(j));
            end
            dstr = datestr(AllSpikes{eid(j),cid(j)}.Header.ctime);
            dstr = dstr(1:6);
            h = AddTitle(DATA, C, 2);
           %h = title(sprintf('E%.1fP%d 1D%.1f 2D%.1f %s',C.exptno,cid(j),C.mahal(4),C.mahal(1),dstr));
        end
        set(gca,'Xtick',[],'Ytick',[]);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        c = get(gca,'Children');
        for k = 1:length(c)
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        end
    end
    if DATA.profiling
        fprintf('took %.2f\n',mytoc(ts));
        profile viewer;
    end
    
function CallAllVPcs(DATA, eid, pid)
    
    if DATA.usealltrials
        args = {'usealltrials'};
    else
        args = {};
    end
    Expts = getappdata(DATA.toplevel,'Expts');
    if ~isempty(Expts)
        args = {args{:}, 'Expts', Expts};
    end
    ts = now;
    AllFullV = getappdata(DATA.toplevel,'AllFullV');
    if strcmp(DATA.allvpcsmode,'fromspikes')
        name = sprintf('%s/Expt%dSpikes',DATA.name,DATA.exptid(eid));
        args = {args{:} 'tchan' pid};
        AllVPcs(name,args{:},'reapply','nocheck');
    elseif ~isempty(AllFullV)
        exs = CellToMat(AllFullV,'exptno');
        id = find(exs == DATA.exptid(eid));
        args = {args{:} 'tchan' pid};
        AllVPcs(AllFullV{id},args{:},DATA.allvpcsmode,'nocheck');
    else
        
        if strncmp(DATA.DataType,'Grid',4)
            name = sprintf('%s/Expt%d.p%dFullV.mat',DATA.name,DATA.exptid(eid),pid);
        else
            name = sprintf('%s/Expt%dFullV.mat',DATA.name,DATA.exptid(eid));
            args = {args{:} 'tchan' pid};
        end
        if exist(name)
            AllVPcs(name,args{:},DATA.allvpcsmode);
        end
    end
    fprintf('Calling AllVPcs took %.1f sec\n',mytoc(ts));

function SetCellFromLine(a,b, cluster, cell)
    axdata = get(gca,'UserData');
    if cluster == 0 && axdata.clnum > 0
        cluster = axdata.clnum;
    end
    if isfield(axdata,'toplevel;')
    DATA = get(axdata.toplevel,'UserData');
    else
    DATA = GetDataFromFig(a);
    end
    if strcmp(cell,'optimize')
        Clusters = getappdata(DATA.toplevel,'Clusters');
        DATA.cluster = Clusters{axdata.eid}{axdata.probe};
        C = OptimizeCluster(DATA);
        Clusters{axdata.eid}{axdata.probe} = C;
        SetFigure(DATA,DATA.tag.xyplot); hold off;
        PlotClusterXY(DATA,C,'shorttitle','tight');
        return;
    elseif strncmp(cell,'allxy',5) || strncmp(cell,'allspk',6)
        DATA.currentpoint(1) = axdata.eid;
        DATA.currentpoint(2) = axdata.probe;
        if strcmp(cell,'allxyprobe')
            DATA = PlotExptsProbe(DATA, 'AllprobeXY');
        elseif strcmp(cell,'allxyexpt')
            DATA = PlotAllProbeXY(DATA, 'AllexptXY');
        elseif strcmp(cell,'allspksexpt')
            DATA = PlotAllProbe(DATA, 'allspks');
        elseif strcmp(cell,'allspksprobe')
            DATA = PlotAllProbe(DATA, 'allprobespks');
        end
        set(DATA.toplevel,'UserData',DATA);
        return;
    elseif strcmp(cell,'duplicate')
        xc = getappdata(DATA.toplevel,'xcCellList');
        a = cluster;
        cluster = xc(a).cl;
        axdata.probe = xc(a).p;
        cell = -axdata.cell;
    elseif strcmp(cell,'spool')
        SpoolSpikes(DATA, [axdata.eid axdata.probe]);
        return;
    elseif strcmp(cell,'findduplicate')
        FindDuplicates(DATA, [axdata.eid axdata.probe axdata.clnum]);
        return;
    elseif strcmp(cell,'getfullv')
        CallAllVPcs(DATA, axdata.eid, axdata.probe)
        return;
    end
    if sum(DATA.selectprobe(:) >0) > 1
        Clusters = getappdata(DATA.toplevel,'Clusters');
        [a,b] = find(DATA.selectprobe > 0);
        for j = 1:length(a)
            DATA = SetCellEntry(DATA, Clusters{a(j)}{b(j)}, a(j),b(j),cluster, cell,'nosave');
        end
        SaveCellList(DATA);
    else
        DATA = SetCellEntry(DATA, [], axdata.eid,axdata.probe,cluster, cell);
    end
    DATA.currentpoint(2) = axdata.probe;
    if isfield(axdata,'celllabels') && ishandle(axdata.celllabels)
        delete(axdata.celllabels);
        axdata.celllabels = AddCellLabels(DATA, axdata.eid, axdata.probe);
        set(gca,'UserData',axdata);
    end
    set(DATA.toplevel,'UserData',DATA);
    PlotCellList(DATA,'showfig');

function DeleteCellFromLine(a,b, cluster, cell)
    axdata = get(gca,'UserData');
    DATA = get(axdata.toplevel,'UserData');
    
    DATA = DeleteCell(DATA, DATA.currentpoint(1),axdata.probe,cluster);
    set(DATA.toplevel,'UserData',DATA);
    PlotCellList(DATA,'showfig');
    
function SetCellFromSubplot(a,b, cell)
    axdata = get(gca,'UserData');
    DATA = get(axdata.toplevel,'UserData');
   
    DATA.CellList(DATA.currentpoint(1),axdata.probe,DATA.cellcluster) = cell;
    DATA.CellDetails.Quality(DATA.currentpoint(1),axdata.probe,DATA.cellcluster) = 4;
    PlotCellList(DATA,'showfig');

function DATA = PlotExptsProbe(DATA, type);
%plots all Expts for one probe
    oldf = gcf;
    DATA = ClearSelections(DATA,0,0);
    Clusters = getappdata(DATA.toplevel,'Clusters');
    expts = 1:length(DATA.exptid);
    nex = length(expts);
    [nr, nc] = Nsubplots(nex);
    p = DATA.currentpoint(2);
    allplot.type = type;
    if strmatch(type,'allprobespks')
        AllSpikes = CheckAllSpikes(DATA, expts, p);
        SetFigure(DATA,DATA.tag.allspikes);
        set(gcf, 'KeyPressFcn',{@KeyPressed,3},'Keyreleasefcn',{@KeyReleased, 3});
    elseif strmatch(type,'alltemplatespks')
        for j = 1:length(DATA.usepeaks)
            AllSpikes = CheckAllSpikes(DATA, expts, abs(DATA.usepeaks(expts(j))));
        end
        SetFigure(DATA,DATA.tag.allspikes);
        set(gcf, 'KeyPressFcn',{@KeyPressed,3},'Keyreleasefcn',{@KeyReleased, 3});
    elseif  strmatch(type,{'AllprobeXY' 'AllTemplateXY' 'AllExptXY'})
        SetFigure(DATA,DATA.tag.allxy);
        set(gcf, 'KeyPressFcn',{@XYKeyPressed,3},'Keyreleasefcn',{@KeyReleased, 3});
    end
%AllXY this probe comes here
    newf = gcf;
    setappdata(gcf,'allplot',allplot);
    fname  = get(gcf,'Name');
    set(gcf,'name',sprintf('Building Plot E%dP%d...',expts(1),p(1)));
    setappdata(0,'alt_is_down',0);
    setappdata(0,'control_is_down',0);
    drawnow;
    figure(newf);
    for j = 1:length(expts)
        mysubplot(nr,nc,j);
        hold off;
        eid = expts(j);
        axdata.type = 'probeexpt';
        axdata.toplevel = DATA.toplevel;
        axdata.eid = eid;
        axdata.probe = p;
        [Clusters, DATA] = CheckClusterLoaded(DATA, eid);

        if strmatch(type,{'allprobespks' 'alltemplatespks'})
            if strmatch(type,'alltemplatespks')
                p = abs(DATA.usepeaks(expts(j)));
            end
            C = Clusters{expts(j)}{p};
            adjid = [eid-1 eid+1];
            adjid = adjid(adjid > 0 & adjid <= size(AllSpikes,1));
            [a,b] = isacell(DATA,eid,p);
            [c,d] = isacell(DATA,adjid,p);
            id = find(c >0);
            if length(id)
                c = c(id(1));
                [e,d] = isacell(DATA,adjid(c),p);
            else
                c = 0;
            end
            if a > 0 || sum(strcmp(type,{'allprobespks' 'alltemplatespks'}))
               [h, x] = QuickSpikes(DATA, AllSpikes{eid,p},C,'showmean');
                axdata.lines = h;
                AddLineContextMenu(DATA, h, eid, p);
                ylim = get(gca,'ylim');
                xlim = get(gca,'xlim');
                if a > 0 || c > 0
                    if (a > 0)
                        th = title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,p,C.mahal(4),C.mahal(1)));
                        set(th,'color','b','fontweight','bold');
                    else
                        th=title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,p,C.mahal(4),C.mahal(1)));
                    end
                    axdata.celllabels = AddCellLabels(DATA, eid, p,'adjacent',adjid);
                else
                    th=title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,p,C.mahal(4),C.mahal(1)));
                end
                set(gca,'ButtonDownFcn',{@HitXYPlot, eid,p });
                TightPlot(gca);
            end
        elseif strmatch(type,{'AllprobeXY' 'AllTemplateXY'})
            if strmatch(type,'AllTemplateXY')
                p = DATA.usepeaks(eid);
                if p == 0
                    p = DATA.currentpoint(2);
                end
                p = abs(p);
            end
            xyplots = PlotClusterXY(DATA,Clusters{eid}{p},'shorttitle','tight','cellid',DATA.currentcluster);
            title(sprintf('E%dP%d:%.1f %.1f',DATA.exptid(eid),p,DistanceMeasure(Clusters{eid}{p},1,DATA.mahaltype),Clusters{eid}{p}.dropi(3)));
            missed = MissedCell(DATA,[eid p]);
            set(xyplots(1),'ButtonDownFcn',{@HitXYPlot, eid, p});
            axdata.lines = xyplots;
            axdata.allxy = 1;
            AddLineContextMenu(DATA, xyplots, eid, p);
            FastAxes(gca);
            if DATA.plottrighist
                AddTrigHist(DATA,Clusters{eid}{p},DATA.currentcluster);
            end

            if missed
                h = get(gca,'title');
                set(h,'color','g');
            end
            c = get(gca,'Children');
            c = setdiff(c,xyplots);
            for k = 1:length(c)
                set(c(k),'ButtonDownFcn',{@HitXYPlot, eid, p});
            end
            AddCellLabels(DATA, eid, p);
            set(gca,'ButtonDownFcn',{@HitXYPlot, eid, p},'UserData',axdata);
        end
        set(gca,'UserData',axdata);
    end
         set(gcf,'name',fname);
         %set(gcf,'renderer',DATA.renderer);
         drawnow;
%set(0,'CurrentFigure',oldf);
figure(oldf);
%gcf


function [S, new] = CheckAllSpikes(DATA, e,p)
    S = getappdata(DATA.toplevel, 'AllSpikes');
    sz = size(S);
    new = 0;
    CheckClusterLoaded(DATA,e);
    for j = 1:length(e)
        if sum(size(e) == size(p)) == 2
            if sz(1) < e(j) || sz(2) < p(j) || ~isfield(S{e(j),p(j)},'values')
                S =  AddSpikes(DATA, S, e(j),p(j));
                new = new+1;
            end
        else
            for k = 1:length(p)
                if sz(1) < e(j) || sz(2) < p(k) || ~isfield(S{e(j),p(k)},'values')
                    S =  AddSpikes(DATA, S, e(j),p(k));
                    new = new+1;
                end
            end
        end
    end
    if new
        setappdata(DATA.toplevel,'AllSpikes',S);
    end

    
function need = ClusterNeedsRefresh(C, I)
    
    need = 0;
    d = dir(I.loadname);
    if strfind(I.loadname,'AutoCluster')
        name = strrep(I.loadname,'AutoCluster','Cluster');
        if exist(name)
            need = 3;
        end
    end
    if d.datenum > I.loadtime
        need = 4;
    end
    if need > 0
        d = dir(strrep(I.loadname,'ClusterTimes','ClusterTimesDetails'));
        if isempty(d) || d.datenum < I.loadtime
            need = 2; %Details not modified = quick save.  
        end
    end

function [Clusters, DATA] = CheckClusterLoaded(DATA,e)
        
    Clusters = getappdata(DATA.toplevel,'Clusters');
    ClusterInfo = getappdata(DATA.toplevel,'ClusterInfo');
    S = getappdata(DATA.toplevel, 'AllSpikes');
    if isempty(Clusters)
        Clusters = {};
    end
    E = getappdata(DATA.toplevel,'Expts');
    if isempty(E)
        E = {};
    end
%Only get GMdata if there are new Clusters. Getting/Freeing GM data is slow    
    GM = {};

    newc = 0;
    newe = 0;
    newsp = 0;

for j = 1:length(e)    
    eid = DATA.exptid(e(j));
    if e(j) > length(Clusters) || isempty(Clusters{e(j)})
        need = 1;
    else
        need = ClusterNeedsRefresh(Clusters{e(j)},ClusterInfo{e(j)});
    end
    
    quickreload = 1;
    if need == 2 && quickreload
        need = 0;
    end
    if need > 0
        if need == 2 %quick Cluster, so rebuild Cluster.clst
%if only AutoClusterTimes existed, then quick save still creates ClusterTimes. So change name now    
            ClusterInfo{e(j)}.loadname = strrep(ClusterInfo{e(j)}.loadname,'AutoClusterTimes','ClusterTimes');
            fprintf('Loading modified %s\n',ClusterInfo{e(j)}.loadname);
            a = load(ClusterInfo{e(j)}.loadname);
            details = ClusterInfo{e(j)};
            details.loadtime = now;
            
            C = Clusters{e(j)};
            D = a.Clusters;
            clfields = {'space' 'xyr' 'sign' 'crit' 'bestspace' 'gmdprime'};
            for p = 1:length(D)
                C{p} = CopyFields(C{p},D{p},clfields);
                clst = ones(size(C{p}.clst));
                [t,tid] = intersect(C{p}.times,D{p}.times);
                clst(tid) = 2;
                for c = 1:length(D{p}.next)
                    if isfield(D{p}.next{c},'times')
                        if length(C{p}.next) < c
                            C{p}.next{c} = D{p}.next{c};
                        else
                            C{p}.next{c} = CopyFields(C{p}.next{c},D{p}.next{c},clfields);
                        end
                        [t,tid] = intersect(C{p}.times,D{p}.next{c}.times);
                        clst(tid) = c+2;
                    end
                end
                C{p}.clst = clst;
            end
            newc = 1;
        else
            if need > 1  %should have loaded something before
                fprintf('Loading Cluster %s\n',ClusterInfo{e(j)}.loadname);
            else
                fprintf('Loading Cluster for Expt %d\n',e(j));
            end
            [C, FullV, details] = LoadCluster(DATA.name,eid,'rawxy','alltimes');
%Check Spikes. If they are loaded for this E,P, then they need reloading
%since the cluster did.
            sz = size(S);
            if sz(1) >= e(j)
                for p  = 1:length(C)
                    if sz(2) >= p && isfield(S{e(j),p},'values')
                        S =  AddSpikes(DATA, S, e(j),p);
                        newsp = newsp+1;
                    end
                end
            end
        end
        if isempty(GM)
            GM = getappdata(DATA.toplevel,'GMfits');
            if isempty(GM)
                GM = {};
            end
        end
        
        for k = 1:length(C)
            C{k}.exptid = e(j);
            if ~isfield(C{k},'bestspace')
                C{k}.bestspace(1) = 0;
                C{k}.bestspace(2) = NaN;
            end
            if ~isfield(C{k},'sign')
                C{k}.sign = 0;
            end
            if ~isfield(C{k},'crit')
                C{k}.crit = 0;
            end
            for cl = 1:length(C{k}.next)
                if isfield(C{k}.next{cl},'space')
                    C{k}.next{cl}.exptid = e(j);
                end
            end
        end
        ClusterInfo{e(j)} = details;
        DATA = ReadClustersIntoData(DATA,C,e(j));
        [Clusters{e(j)},GM{e(j)}] = CondenseCluster(C);
        newc = 1;
    end
    if e(j) > length(E) || isempty(E{e(j),1})
        if exist(DATA.exptname)  %one expt file
            [DATA, E] = LoadExpts(DATA);           
        else
            [E{e(j),1}, DATA] = LoadExpt(DATA,eid);
        end
        newe = 1;
    end
end
    
if newe
    E = SetExptTimeOffset(E);
    setappdata(DATA.toplevel,'Expts',E)
end
if newsp
    setappdata(DATA.toplevel,'AllSpikes',S);
end
if newc
    setappdata(DATA.toplevel,'Clusters',Clusters)
    setappdata(DATA.toplevel,'ClusterInfo',ClusterInfo)
    setappdata(DATA.toplevel,'GMfits',GM)
end
if newc || newe
    set(DATA.toplevel','UserData',DATA);
    DATA.modified = 1;
end

function S = AddSpikes(DATA, S, e,p, varargin)
    verbose = 1;
    setdata = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'setdata',3)
            setdata = 1;
        end
        j = j+1;
    end
   
    mnk = GetMonkeyName(DATA.name);
    if ~isfield(DATA,'exptid')
        DATA.exptid = DATA.CellDetails.exptids;
    end
    [a,b] = fileparts(DATA.name);
    if isempty(b) %name is a directory
        [a,b] = fileparts(a);
    end
    b = regexprep(b,'[MG]([0-9]*)(.*)','$0');
    [c,d] = fileparts(a);
    xs = '';
    if isempty(S) && isappdata(DATA.toplevel,'AllSpikes')
        S = getappdata(DATA.toplevel,'AllSpikes');
    end
    if rem(DATA.exptid(e),1) > 0.001
        xs = 'a';
    end
    name = [DATA.name '/Spikes/' mnk b '.p' num2str(p)  't' num2str(floor(DATA.exptid(e))) xs '.mat'];
    if verbose > 0
        fprintf('Reading %s at %s\n',name,datestr(now));
    end
    S{e,p} = ReadSpikeFile(name);
    S{e,p}.probe = p;
    if setdata
        setappdata(DATA.toplevel,'AllSpikes',S);
    end
            
     
    
function PlotAllExpt(DATA, type)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    [nr,nc] = Nsubplots(length(Clusters));
    SetFigure(DATA,DATA.tag.allexpt);
    cr = [];
    if strcmp(type,'exptdpim')
        for e = 1:length(Clusters)
            dprange = [];
            for k = length(Clusters{e}):-1:1
                dprange(k,:) = minmax(Clusters{e}{k}.MeanSpike.dp(:));
            end
            crange(e,:) = minmax(dprange(:));
        end
        cr = [min(crange(:,1)) max(crange(:,2))];
    end
    for e = 1:length(Clusters)
        mysubplot(nr,nc,e);
        hold off; 
        if strmatch(type,{'exptim' 'exptdpim'}) 
            h = PlotProbeMeans(Clusters{e},type,'caxis',cr);
            crange(e,:) = caxis;
            set(h,'buttondownfcn',{@HitExptPlot,'expt', e})
        end
        set(gca,'buttondownfcn',{@HitExptPlot,'expt', e},'xtick',[],'ytick',[]);
    end
    if isempty(cr)
    cr = [min(crange(:,1)) max(crange(:,2))];
    for e = 1:length(Clusters)
        mysubplot(nr,nc,e);
        caxis(cr);
    end
    end

function PlotAllProbeMean(DATA, type, varargin)
    if strcmp(type,'AllMeanIm')
        type = 'imageonly';
    elseif strcmp(type,'AllMean')
        type = 'lineonly';
    elseif strcmp(type,'AllExptIm')
        type = 'exptim';
    end
        
    Clusters = getappdata(DATA.toplevel,'Clusters');
    eid = DATA.currentpoint(1);
    [nr,nc] = Nsubplots(length(Clusters{eid}));
    SetFigure(DATA,DATA.tag.allxy);
    cmenu = AddCellContextMenu(DATA, 'subplot');
    axdata.toplevel = DATA.toplevel;
    for j = 1:length(Clusters{eid})
        mysubplot(nr,nc,j);
        hold off; 
        PlotMeanSpike(Clusters{eid}{j},0,1,type);
        set(gca,'Xtick',[],'Ytick',[]);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        axdata.probe = j;
        axdata.eid = eid;
        set(gca,'UIContextMenu',cmenu, 'UserData', axdata);
    end

function PlotAllExptProbeMean(DATA, type, varargin)
    if strcmp(type,'AllMeanIm')
        type = 'imageonly';
    elseif strcmp(type,'AllExptMeanIm')
        type = 'smallimage';
    elseif strcmp(type,'AllMean')
        type = 'lineonly';
    elseif strcmp(type,'AllExptIm')
        type = 'exptim';
    end
        
    Clusters = getappdata(DATA.toplevel,'Clusters');
    eid = DATA.currentpoint(1);
    SetFigure(DATA,DATA.tag.allxy);
    cmenu = AddCellContextMenu(DATA, 'subplot');
    axdata.toplevel = DATA.toplevel;
    
    if length(DATA.selectexpts)
        expts = DATA.selectexpts;
    else
        expts = 1:length(Clusters);
    end
    if length(DATA.proberange)
        probelist = DATA.proberange;
    else
       probelist = 1:DATA.nprobes;
    end
    nex = length(expts);
    for eid = 1:length(expts)
    for j = 1:length(Clusters{eid})
        mins(eid,j) = min(Clusters{expts(eid)}{j}.MeanSpike.ms(:));
        maxs(eid,j) = max(Clusters{expts(eid)}{j}.MeanSpike.ms(:));
    end
    end
        
    clim = [min(mins(:)) max(maxs(:))];
    ClearPlot;
    for e = 1:length(expts)
    for j = 1:length(probelist)
        eid = expts(e);
        p = probelist(j);
        mysubplot(nex,length(probelist),j + (e-1)*length(probelist));
        hold off; 
        PlotMeanSpike(Clusters{eid}{p},0,1,type);
        set(gca,'Xtick',[],'Ytick',[],'clim',clim);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        delete(h);
        axdata.probe = j;
        axdata.eid = eid;
        set(gca,'UIContextMenu',cmenu, 'UserData', axdata);
    end
    drawnow;
    end

    
function cmenu = AddCellContextMenu(DATA, type)
    cmenu = uicontextmenu;
    
    if strcmp(type,'subplot')
    for j = 1:50
        uimenu(cmenu,'label',sprintf('Cell %d',j),'Callback',{@SetCellFromSubplot, j});
    end
    else
    
    uimenu(cmenu,'label',sprintf('Spike%d',j-1),'foregroundcolor',DATA.colors{j});
    for k = 1:50
        c = uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, j-1,  k});
        if k == 1
            set(c,'separator','on');
        end
    end
    end
    
function cmenu = AddContextMenu(DATA, type)
    cmenu = uicontextmenu;
    
    if strcmp(type,'subplot')
        uimenu(cmenu,'label','spool','Callback',{@SetCellFromLine, j-1,  'spool'});
        uimenu(cmenu,'label','->FullV','Callback',{@SetCellFromLine, j-1,  'getfullv'});
        uimenu(cmenu,'label','cut cluster1 (&Line)','Callback',{@XYCluster, 'lines', 'useaxdata'});
        uimenu(cmenu,'label','cut ellipse&1','Callback',{@XYCluster, 'ellipse1', 'useaxdata'});
        uimenu(cmenu,'label','cut ellipse&2','Callback',{@XYCluster, 'ellipse2', 'useaxdata'});
        uimenu(cmenu,'label','cut ellipse&3','Callback',{@XYCluster, 'ellipse3', 'useaxdata'});
        uimenu(cmenu,'label','cut ellipse&4','Callback',{@XYCluster, 'ellipse4', 'useaxdata'});
        uimenu(cmenu,'label','No cutting','Callback',{@XYCluster, 'nocut'});
        uimenu(cmenu,'label','save (&Quick)','Callback',{@XYCluster, 'quicksave'});
        uimenu(cmenu,'label','Optimize','Callback',{@XYCluster, 'optimize'});
        uimenu(cmenu,'label','test line','Callback',{@XYCluster, 'testline'});
    elseif strcmp(type,'cellplot')
        tm = uimenu(cmenu,'label','Ex0 P0:','Tag','Title');
        uimenu(cmenu,'label','spool','Callback',{@SetCellFromLine, 1,  'spool'});
        uimenu(cmenu,'label','->AllVPcs','Callback',{@SetCellFromLine, 1,  'getfullv'});
        uimenu(cmenu,'label','Find Duplicates','Callback',{@SetCellFromLine, 1,  'findduplicate'});
        pm = uimenu(cmenu,'label','Plots');
        uimenu(pm,'label','All XY This Probe','Callback',{@SetCellFromLine, 1,  'allxyprobe'});
        uimenu(pm,'label','All XY This Expt','Callback',{@SetCellFromLine, 1,  'allxyexpt'});
        uimenu(pm,'label','All Spks This Probe','Callback',{@SetCellFromLine, 1,  'allspksprobe'});
        uimenu(pm,'label','All Spks This Expt','Callback',{@SetCellFromLine, 1,  'allspksexpt'});
        pm = uimenu(cmenu,'label','Other');
        uimenu(pm,'label','Comment','Callback',{@TagMenu, 'comment'});
        uimenu(pm,'label','Flip Line Crit','Callback',{@TagMenu, 'flipline'});
        cells = unique(DATA.CellList(:));
        cells = cells(cells > 0);
        if max(cells) > 20
            nc = max(cells)+2;
        else
            nc = 20;
        end
        for k = 1:nc
            c = uimenu(tm,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, 0,  k});
            if ismember(k,cells)
                if isfield(DATA,'cellcolors') && k <= length(DATA.cellcolors)
                set(c,'foregroundcolor',DATA.cellcolors{k});
                else
                set(c,'foregroundcolor','r');
                end
            end
            if k == 1
                set(c,'separator','on');
            end
        end
    else
    
    uimenu(cmenu,'label',sprintf('Spike%d',j-1),'foregroundcolor',DATA.colors{j});
    for k = 1:50
        c = uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, j-1,  k});
        if k == 1
            set(c,'separator','on');
        end
    end
    end
    
function cmenu = AddLineContextMenu(DATA, h, e, p, varargin)

    setcl = 0;
    setdup = 0;
    pid = 0;
    starth = 2; %default is to ignore first line
    j = 1;
    
    while j <= length(varargin)
        if strncmpi(varargin{j},'cellnum',5)
            j = j+1;
            xc = getappdata(DATA.toplevel,'xcCellList');
            pid = varargin{j};
            setcl = xc(pid).cl;
            a = xc(pid).p;
            starth = 1;
        elseif strncmpi(varargin{j},'duplicate',5)
            j = j+1;
            setdup = varargin{j};
        end
        j = j+1;
    end
    cells = unique(DATA.CellList(:));
    cells = cells(cells > 0);
    for j = starth:length(h)
        if h(j) > 0 && ishandle(h(j))
            cmenu = uicontextmenu;
            if setcl
                k = DATA.CellList(e,p,setcl);
            elseif j > size(DATA.CellList,3)+1 || e > size(DATA.CellList,1)
                k = 0;
            else
                k = DATA.CellList(e,p,j-1);
            end
            if k > 0
                uimenu(cmenu,'label',sprintf('Spike%d Cell%d',j-1,k),'foregroundcolor',get(h(j),'color'));
                c = uimenu(cmenu,'label','clear','Callback',{@DeleteCellFromLine, j-1,  k});
            else
                uimenu(cmenu,'label',sprintf('Spike%d',j-1),'foregroundcolor',get(h(j),'color'));
            end
            if setdup
                uimenu(cmenu,'label',sprintf('Duplicate of %d',setdup),'Callback',{@SetCellFromLine, pid,  'duplicate'});
            end
            uimenu(cmenu,'label','&spool','Callback',{@SetCellFromLine, j-1,  'spool'});
            uimenu(cmenu,'label','&optimize','Callback',{@SetCellFromLine, j-1,  'optimize'});
            uimenu(cmenu,'label','->FullV','Callback',{@SetCellFromLine, j-1,  'getfullv'});
        for k = 1:50
            c = uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, j-1,  k});
            if ismember(k,cells)
                set(c,'Foregroundcolor','g');
            end
            if k == 1
                set(c,'separator','on');
            end
        end
        set(h(j),'uicontextmenu',cmenu);
        end
    end

function h = PlotProbeMeans(C,type, varargin)
    
    mscale = [0 0.5];
    maxmahal = 5;
    cr = [];
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'caxis',3)
            j = j+1;
            cr = varargin{j};
            mscale(1) = cr(1);
            mscale(2) = maxmahal/(cr(2)-cr(1));
            mscale(2) = (cr(2)-cr(1))/maxmahal;
        end
        j = j+1;
    end
    for j = 1:length(C)
        if strcmp(type,'exptim')
            nv = size(C{j}.MeanSpike.ms,2);
            im(j,1:nv) = C{j}.MeanSpike.ms(j,:);
        elseif strcmp(type,'exptdpim')
            row = ceil(size(C{j}.MeanSpike.dp,1)/2);
            nv = size(C{j}.MeanSpike.dp,2);
            im(j,1:nv) = C{j}.MeanSpike.dp(j,:);
            dp = min([maxmahal max(C{j}.mahal([1 4]))]);
            im(j,end) = (dp+mscale(1)) .*mscale(2);
            im(j,end) = (dp.*mscale(2))+mscale(1);
        end
    end
    h = imagesc(im);

    if length(cr)
        caxis(cr);
    end
    
function DATA = PlotAllProbe(DATA, type)
%Plots All SPks, for probes or expts
    Clusters = getappdata(DATA.toplevel,'Clusters');
    AllSpikes = {};
    axdatatype = [];
    allprobes = 0;
    if ~strcmp(type,'selectspks')
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    end

    if DATA.profiling
        profile on;
    end
    eid = DATA.currentpoint(1);
    %'exptspks', 'allspks' plot all probes for one expt
     if strmatch(type,{'allspks' 'allexptspks' 'allcellspks' 'spoolall' 'exptspks'})
         AllSpikes = CheckAllSpikes(DATA, eid, 1:DATA.nprobes);
         np =length(Clusters{eid});
         pid = 1:np;
         expts = eid;
         allprobes = 1;
     elseif strmatch(type,{ 'allprobespks' 'AllProbeMean'})
         np =length(DATA.exptid);
         expts = 1:np;
         if strmatch(type,{ 'allprobespks'})
         AllSpikes = CheckAllSpikes(DATA, expts, DATA.currentpoint(2));
         else
         end
         pid = ones(size(expts)) .* DATA.currentpoint(2);
         axdatatype = 'probeexpt';
     elseif strmatch(type,{'selectspks'})
         if isfield(DATA,'usepeaks') && sum(DATA.usepeaks > 0)
            expts = find(DATA.usepeaks > 0);
            pid = DATA.usepeaks(expts);
         else
             [expts, pid] = find(DATA.selectprobe > 0);
         end
         np = length(expts);
         for j = 1:np
             AllSpikes = CheckAllSpikes(DATA, expts(j), pid(j));
         end
     end
     %get Clusters again in case changed
    Clusters = getappdata(DATA.toplevel,'Clusters');
    Expts = getappdata(DATA.toplevel,'Expts');
    if DATA.profiling
        profile viewer;
    end
    DATA.voffset = SetVOffset(DATA, AllSpikes,expts(1));

    
    DATA.Expt = Expts{eid,1};
    if np == 24 && strcmp(DATA.ArrayConfig.type,'12x2')
        nr = 4;
        nc = 6;
    else
        [nr,nc] = Nsubplots(np);
    end
     if strcmp(type, 'spoolall')
         DATA.plotspk.allprobes =1;
            SpoolAllProbes(DATA, eid, AllSpikes(eid,:),Clusters{eid});
            return;
     elseif strcmp(type, 'spooleverything')
         DATA.plotspk.allprobes =1;
         for j = 1:length(DATA.exptid)
             AllSpikes = CheckAllSpikes(DATA, j, 1:DATA.nprobes);
            stopped = SpoolAllProbes(DATA, j, AllSpikes(j,:),Clusters{j});
            if stopped
                return;
            end
         end
         return;
     elseif strcmp(type, 'spooleverycell')
         for j = 1:length(DATA.exptid)
             cid = isacell(DATA, j, 1:DATA.nprobes);
             AllSpikes = CheckAllSpikes(DATA, j, cid);
            stopped = SpoolAllProbes(DATA, j, AllSpikes(j,cid),Clusters{j});
            if stopped
                return;
            end
         end
         return;
     elseif strcmp(type, 'spoolcells')
         cid = isacell(DATA, eid, 1:DATA.nprobes);
             AllSpikes = CheckAllSpikes(DATA, eid, find(cid));
            SpoolAllProbes(DATA, eid, AllSpikes(eid,cid),Clusters{eid});
            return;
     end
     if strmatch(type,{'AllProbeMean'})
         mtype = DATA.plotmeantype;
     else
         mtype = type;
     end
     if strmatch(type,{'allspks' 'allexptspks' ...
              'allcellspks' 'allprobespks' 'selectspks' 'AllProbeMean'})
         F = SetFigure(DATA,DATA.tag.allspikes);
         subplot(1,1,1); %not clf - that wipes menus etc
     else
         SetFigure(DATA,DATA.tag.allxy);
     end

     cmenu = uicontextmenu;
    for j = 1:50
        uimenu(cmenu,'label',sprintf('Cell %d',j),'Callback',{@SetCellFromSubplot, j});
    end
    axdata.toplevel = DATA.toplevel;
    if allprobes && isfield(DATA.ArrayConfig,'X')
        plotpos = SetPlotPos(DATA, np, nr, nc);
        scalebar = 0.005;
    else
        plotpos = 1:np;
        scalebar = 1;
    end

     for j = 1:np
         mysubplot(nr,nc,plotpos(j),'tight','width',0.95);
        hold off; 
        p = pid(j);
        if strmatch(type,{'allspks' 'allexptspks' 'allcellspks' 'selectspks' 'exptspks' 'allprobespks'})
            if sum(strcmp(type,{'selectspks' 'allprobespks'}))
                C = Clusters{expts(j)}{pid(j)};
                eid = expts(j);
            else
                C = Clusters{eid}{j};
            end
            adjid = [eid-1 eid+1];
            adjid = adjid(adjid > 0 & adjid <= size(AllSpikes,1));
            [a,b] = isacell(DATA,eid,pid(j));
            [c,d] = isacell(DATA,adjid,pid(j));
            id = find(c >0);
            if length(id) && length(adjid) > max(c)
                c = c(id(1));
                [e,d] = isacell(DATA,adjid(c),p);
            else
                c = 0;
            end
            if a > 0 || sum(strcmp(type,{'allspks' 'selectspks' 'allprobespks'}))
            [h, x] = QuickSpikes(DATA, AllSpikes{eid,pid(j)},C,'showmean','scalebar',scalebar);
            axdata.lines = h;
            AddLineContextMenu(DATA, h, eid, pid(j));
            if isfield(x,'meanh')
                AddLineContextMenu(DATA, x.meanh, eid, pid(j));
            end
            if isfield(AllSpikes{eid,pid(j)},'Header')
            dstr = datestr(AllSpikes{eid,pid(j)}.Header.ctime);
            dstr = dstr(1:6);
            else
                dstr = '';
            end
            
            if a > 0 || c > 0
                th=title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,p,C.mahal(4),C.mahal(1)));
                axdata.celllabels = AddCellLabels(DATA, eid, pid(j),'adjacent',adjid);
            else
                if c > 0
                    PC = Clusters{adjid(c)}{j};
                    axdata.celllabels = AddCellLabels(DATA, adjid(c), j,'parentheses');
                    hold on;
                    for k = 1:min([length(d) length(h)-1])
                        if h(k+1) > 0 & ishandle(h(k+1))
                            color = get(h(k+1),'color');
                            color = 1-color;
                            if k > 1
                                offset = std(PC.next{k-1}.MeanSpike.ms(j,:))./4;
                                plot(PC.next{k-1}.MeanSpike.ms(j,:)-offset,'--','color',color,'linewidth',2);
                            else
                                offset = std(PC.MeanSpike.ms(j,:))./4;
                                plot(PC.MeanSpike.ms(j,:)-offset,'--','color',color,'linewidth',2);
                            end
                        end
                    end
                end
                th=title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,p,C.mahal(4),C.mahal(1)));
            end
            nid = find(~ismember(d,b) & d > 0);
            if diff(size(nid)) < 0
                nid = nid';
            end
                if c > 0
            for k = nid;
               hold on;
               PC = Clusters{adjid(c)}{p};
               mspid = GetMeanSpikeProbe(PC, p);
               if k < length(h) & h(k+1) > 0 & ishandle(h(k+1))
                            color = get(h(k+1),'color');
                            color = 1-color;
                            if k > 1
                                offset = std(PC.next{k-1}.MeanSpike.ms(mspid,:))./4;
                                plot(PC.next{k-1}.MeanSpike.ms(mspid,:)-offset,'--','color',color,'linewidth',2);
                            else
                                offset = std(PC.MeanSpike.ms(mspid,:))./4;
                                plot(PC.MeanSpike.ms(mspid,:)-offset,'--','color',color,'linewidth',2);
                            end
                        end
            end
                end
            end
        else
            eid = expts(j);
            h = PlotMeanSpike(Clusters{expts(j)}{p},0,1,mtype);
            axdata.lines = h;
        end
        axdata.probe = p;
        axdata.eid = eid;
        axdata.type = axdatatype;
        set(gca,'Xtick',[],'Ytick',[]);
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid, p});
        set(gca,'UserData',axdata,'UIContextMenu',AddContextMenu(DATA, 'subplot'));
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl(j,:) = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(j,2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
     end
    
    
     if DATA.plotspk.onescale
         c = get(gcf,'Children');
         ylim = minmax(yl(:));
         for j = 1:length(c)
             if strcmp(get(c(j),'type'),'axes')
                 set(c(j),'ylim',ylim);
             end
         end
     end
    if DATA.plotspk.showcellmeans
            for j = 1:length(Clusters{eid})
                [a,b] = isacell(DATA,eid,j);
                b = b(b>0);
                for k = 1:length(b)
                    ms = Clusters{eid}{j}.MeanSpike.ms;
                    amp = std(ms');
                    id = find(amp > amp(j)./5);
                    for p = id
                        mysubplot(nr,nc,p);
                        hold on;
                        plot(ms(p,:),'linewidth',2,'color',DATA.colors{b(k)+1});
                    end
                end
            end
    end

function mspid = GetMeanSpikeProbe(C, p)
    if size(C.MeanSpike.ms,1) == 1
        mspid = 1;
    else
        mspid = p;
    end
               
function plotpos = SetPlotPos(DATA,np, nr, nc)
    if isfield(DATA.ArrayConfig,'X')
        nx = length(unique(DATA.ArrayConfig.X));
        if ismember(nx,[np 1]) %simple linear array
            plotpos = 1:np;
        else
            if length(unique(DATA.ArrayConfig.X)) == 2 && np == 24
                plotpos = [1 7 2 8 3 9 4 10 5 11 6 12 13 19 14 20 15 21 16 22 17 23 18 24];
            else
                for j = 1:np
                    plotpos(j) = DATA.ArrayConfig.Y(j) + (DATA.ArrayConfig.X(j)-1).*nr;
                end
            end
        end
    else
        plotplos = 1:np;
    end
    
function ReplotXY(a,b,eid, probe,cid)
    DATA = GetDataFromFig(a);  
    
    if ishandle(a)
    cmenu = get(a,'UIContextMenu');
    else
        cmenu = [];
    end
    callback = get(gca,'ButtonDownFcn');
    axdata = get(gca,'UserData');
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    if bt == 2  && ~isempty(cmenu)
        return;
    end
    if DATA.elmousept.down > 0
        return;
    end
    if cid ~= DATA.currentcluster
        newcluster = 1;
    else
        newcluster = 0;
    end
    DATA.currentcluster = cid;
    Clusters = getappdata(DATA.toplevel,'Clusters');
    fprintf('Replotting E%dP%dC%d\n',eid,probe,cid);
    PlotClusterXY(DATA,Clusters{eid}{probe},'shorttitle','cellid',cid,'tight');
    if DATA.plottrighist
        AddTrigHist(DATA,Clusters{eid}{probe},DATA.currentcluster);
    end
    set(gca,'Xtick',[],'Ytick',[],'ButtonDownFcn',callback,'UserData',axdata);
    AddCellLabels(DATA, eid, probe, 'NW');
    if newcluster
        set(DATA.toplevel,'UserData',DATA);
        F = gcf;
        if DATA.plothist
            SetFigure(DATA,DATA.tag.hist,'front');
            PlotClusterHistogram(DATA, Clusters{eid}{probe}, DATA.refitgm,'cluster', DATA.currentcluster);
            figure(F);
        end
    end
    
function DATA = PlotAllProbeXY(DATA,varargin)
    %plots all probes for one expt
    oneprobe = 0;
    selectprobes = [];
    axdatatype = [];
    plotexpt = [];
    j =1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'oneprobe',7)
            oneprobe = 1;
        elseif strncmpi(varargin{j},'expt',4)
            j = j+1;
            plotexpt = varargin{j};
        elseif strncmpi(varargin{j},'select',6)
            oneprobe = 2;
            if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            selectprobes = varargin{j};
            end
        end
        j = j+1;
    end
    
    
    Clusters = getappdata(DATA.toplevel,'Clusters');
    eid = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if length(plotexpt)
        pid = 1:DATA.nprobes;
        expts = ones(size(pid)).*plotexpt;
        np = DATA.nprobes;
    elseif oneprobe == 1
        np = length(Clusters);
        expts = 1:np;
        pid = ones(size(expts)) .* DATA.currentpoint(2);
    elseif oneprobe > 1
        if isempty(selectprobes) 
            if isfield(DATA,'usepeaks') && sum(DATA.usepeaks>0) > 1
                expts = find(DATA.usepeaks > 0);
                pid = DATA.usepeaks(expts);
            else
                [expts, pid] = find(DATA.selectprobe > 0);
            end
        else
        [expts, pid] = find(selectprobes > 0);
        end
        np = length(expts);
    else %All XY one expt
        [Clusters, DATA] = CheckClusterLoaded(DATA, eid);
        if isfield(DATA.ArrayConfig,'X')
            np = length(DATA.ArrayConfig.X);
        else
            np = length(Clusters{eid});
        end
        pid = 1:np;
        expts = ones(size(pid)) .* eid;
    end
    if strcmp(DATA.ArrayConfig.type,'12x2')
        nr = 4;
        nc = 6;
    else
        [nr,nc] = Nsubplots(np);
    end
    SetFigure(DATA,DATA.tag.allxy);
    set(gcf,'Name','Building all XY...');
    drawnow;
    plotpos = SetPlotPos(DATA, np, nr, nc);
    fprintf('Building All XY\n');
    ClearPlot();
    for j = 1:np
        if isfield(DATA.ArrayConfig,'X')
            mysubplot(nr,nc,plotpos(j),'tight');            
        else
            mysubplot(nr,nc,j);
        end
        hold off; 
        if oneprobe == 1
            xyplots = PlotClusterXY(DATA,Clusters{j}{p},'vshorttitle','tight');
            missed = MissedCell(DATA,[j p]);
        elseif oneprobe > 1
            xyplots = PlotClusterXY(DATA,Clusters{expts(j)}{pid(j)},'vshorttitle','tight');
            missed = MissedCell(DATA,[expts(j) pid(j)]);
        else
            xyplots = PlotClusterXY(DATA,Clusters{eid}{j},'vshorttitle','tight');
            missed = MissedCell(DATA,[eid j]);
        end
        if ~isempty(xyplots)
            set(xyplots(1),'ButtonDownFcn',{@HitXYPlot, eid, j});
        end
        set(gca,'Xtick',[],'Ytick',[]);
        if DATA.plottrighist
            AddTrigHist(DATA,Clusters{expts(j)}{pid(j)},DATA.currentcluster);
        end

        if missed
            h = get(gca,'title');
            set(h,'color','g');
        end
        c = get(gca,'Children');
        c = setdiff(c,xyplots);
        for k = 1:length(c)
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid, pid(j)});
        end
        AddCellLabels(DATA, expts(j), pid(j));
                axdata.type = axdatatype;
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid, pid(j)});
            
    end
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    set(gcf,'Name',DATA.tag.allxy);

function handles = AddCellLabels(DATA, eid, probe, varargin)

    handles = [];
    adjid = [];
    yl = get(gca,'Ylim');
    xl = get(gca,'Xlim');
    braces = 0;
    position = [1 0];
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'adjacent',7)
            j = j+1;
            adjid = varargin{j};
        elseif strncmpi(varargin{j},'parenth',7)
            braces = 1;
        elseif strncmp(varargin{j},'NW',2)
            position = [0 1];
        end
        j = j+1;
    end

    if position(2) > 0.5
        vstep = -diff(yl)/12;
    else
        vstep = diff(yl)/12;
    end
        
    position(1) = xl(1)+ diff(xl) * position(1);
    position(2) = yl(1)+ diff(yl) * position(2);

    if size(DATA.CellList,1) >= eid
        cid = find(DATA.CellList(eid,probe,:) > 0);
        nh = 0;
        for c = 1:length(cid)
            str = sprintf('Cell %d',DATA.CellList(eid,probe,cid(c)));
            if braces
                str = ['(' str ')'];
            end
            t = text(position(1),position(2)+vstep*(c-1),str,...
                'fontweight','bold','VerticalAlignment','bottom','HorizontalAlignment','right');
            if DATA.plot.density
                set(t,'color','w');
            else
                set(t,'color',DATA.colors{cid(c)+1},'ButtonDownFcn',{@ReplotXY, eid, probe, cid(c)});
            end
            handles(c) = t;
            nh = c;
        end
        aid = find(DATA.CellList(eid,probe,:) == 0 & sum(DATA.CellList(adjid,probe,:),1) > 0);
        for c = 1:length(aid)
            cellid = DATA.CellList(adjid,probe,aid(c));
            cellid = cellid(cellid > 0);
            str = sprintf('(E%dCell %d)',adjid(1),cellid(1));
            t = text(xl(2),yl(1)+vstep*(nh+c-1),str,...
                'fontweight','bold','VerticalAlignment','bottom','HorizontalAlignment','right');
            if DATA.plot.density
                set(t,'color','w');
            else
                set(t,'color',DATA.colors{aid(c)+1},'ButtonDownFcn',{@ReplotXY, eid, probe, aid(c)});
            end
            handles(c+nh) = t;
        end

    end
        
function missed = MissedCell(DATA, pos)    
    e = pos(1);
    p = pos(2);
    fitd = 0;
    if isfield(DATA,'GaussFitdp') && size(DATA.GaussFitdp,1) >= e
        fitd = DATA.GaussFitdp(e,p,2);
        if sum(DATA.gmfitpos(e,p,:)) < 2
            fitd = 0;
        end
    end
    d = max(DATA.mahal(e,p,[1 3]));
    if isacell(DATA,e,p) == 0 && ((d > 2 && abs(fitd) > 1) || abs(fitd) > 3) && DATA.peakdiff(e,p) == 0 && DATA.dropi(e,p,1) > 1
        missed = 1;
    else
        missed = 0;
    end

    
function HitExptPlot(src, b, type, e)
    DATA = GetDataFromFig(src);    

    DATA.currentpoint(1) = e;
    DATA = ClearSelections(DATA,0,0);
    if DATA.plotspks
        PlotAllProbe(DATA,'allspks');
    end
    set(DATA.toplevel,'UserData',DATA);

 function HitXYPlot(src, b, e,p)
    DATA = GetDataFromFig(src);    
    cmenu = get(src,'UIContextMenu');
    axdata = get(gca,'UserData');
    if isfield(axdata,'allxy') && DATA.allclustering == 1
        return;
    end
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    if bt == 2  && ~isempty(cmenu)
        return;
    end
    if ~isfield(DATA,'selectprobe')
        DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    end
    if e > size(DATA.selectprobe,1) ||  p > size(DATA.selectprobe,2)
        DATA.selectprobe(e,p) = 1;
    end
    DATA.selectprobe(e,p) = ~DATA.selectprobe(e,p);
    if DATA.selectprobe(e,p) && bt == 3 %shift select
        if isfield(axdata,'type')
            fg = get(gca,'Parent');
            c = get(fg,'Children');
            types = get(c,'type');
            axid = strmatch('axes',types);
            if strcmp(axdata.type,'probeexpt')
                for j = 1:length(axid)
                    a = get(c(axid(j)),'UserData');
                    exlist(j) = a.eid;
                end
                [b, a] = sort(exlist);
                axid = axid(a); %sorte by exlist
                id = find(DATA.selectprobe(:,p) > 0);
                a = max(id(id < e));
                if isempty(a)
                    a = min(id(id > e));
                else
                    DATA.selectprobe(a:e,p) = 1;
                    ax = gca;
                    for j = a:e
                        set(c(axid(j)),'color',[0 0 0]);
                    end
                    axes(ax);
                end
            end
        end
    end
    if DATA.selectprobe(e,p)
        set(gca,'color',[0 0 0]);
        if bt ~= 3
        DATA = ShowData(DATA,e,p,'oneprobe');
        end
    else
        set(gca,'color',[1 1 1]);
    end
    DATA.currentpoint = [e p];
    SetFigure(DATA,DATA.tag.all);
%    DATA.markh = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),'color','w');
    set(DATA.toplevel,'UserData',DATA);
    
function TightPlot(ax)
    t = get(ax,'title');
    ylim = get(ax,'ylim');
    xlim = get(ax,'xlim');
    set(t,'position',[mean(xlim) ylim(2)],'HorizontalAlignment','center','VerticalAlignment','top');
    set(ax,'Xtick',[],'Ytick',[]);
    
function PlotAllCellXY(DATA)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');
    
    ts = now;
    if DATA.profiling
        profile on;
    end
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);

    
    [eid, id] = sort(eid);
    cid = cid(id);
    clid = clid(id);
    cid = 1+mod(cid-1,DATA.nprobes);
    [nr,nc] = Nsubplots(length(eid));
    SetFigure(DATA,DATA.tag.allxy);
    subplot(1,1,1);
    colors = mycolors('spkcolors');
    for j = 1:length(eid)
        mysubplot(nr,nc,j);
        hold off; 
        plots = PlotClusterXY(DATA,Clusters{eid(j)}{cid(j)},'shorttitle','cellid',clid(j),'allxy');
        if clid(j) ~= 1 && length(plots) > 1
            np = clid(j)+1;
            if length(plots) >= np
                set(plots(np),'color','r');
                set(plots(2),'color',colors{np});
            else
                set(plots(2),'color','k');
            end
        end
        if DATA.plottrighist
            AddTrigHist(DATA,Clusters{eid(j)}{cid(j)}, DATA.currentcluster);
        end
        set(gca,'Xtick',[],'Ytick',[]);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        color = 'b';
        if max(Clusters{eid(j)}{cid(j)}.mahal([1 3])) < 3
            color = 'r';
        elseif max(Clusters{eid(j)}{cid(j)}.mahal([1 3])) < 2
            color = 'k';
        end
           
        if length(plots) > 1
            ph = plots(2:end);
        else
            ph = [];
        end
        set(h,'fontweight','normal','fontsize',12,'color',color);
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        c = get(gca,'Children');
        for k = 1:length(c)
            if ~ismember(c(k),ph);
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
            end
        end
        AddCellLabels(DATA, eid(j), cid(j));
    end
    if DATA.profiling
        fprintf('Took %.2f\n',mytoc(ts));
        profile viewer;
    end

    
    function h= DrawLine(E,varargin)

        if ~isfield(E,'pos') || isfield(E,'crit')  %pos can get out of data when save from plotclusters
            xy= xyrotate([E.crit E.crit],get(gca,'ylim'),-E.angle);
            E.pos = [xy(1,1), xy(1,2) xy(2,1) xy(2,2);];
        end
x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
if length(E.pos) > 7
    xp = E.pos([5 7]);
    yp = E.pos([6 8]);
else
xp = [mean(x) + diff(y)/4 mean(x)-diff(y)/4];
yp = [mean(y) - diff(x)/4 mean(y)+diff(x)/4];
end
if isfield(E,'h') & ishandle(E.h(1)) 
    set(E.h(1),'Xdata',x,'Ydata',y);
    h = E.h;
    if length(E.h) > 1 && ishandle(E.h(2))
        set(E.h(2),'Xdata',xp,'Ydata',yp);
    end
else
    hold on;
    h(1) = plot(real(x),real(y),varargin{:});
    h(2) = plot(xp,yp,'r:');
    hold off;
end


function [nx, ny] = Data2Norm(x,y)
    
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    
    nx = (x-xl(1))./diff(xl);
    ny = (y-yl(1))./diff(xl);
    
function h= testDrawLine(E,varargin)

        if ~isfield(E,'pos') || isfield(E,'crit')  %pos can get out of data when save from plotclusters
            xy= xyrotate([E.crit E.crit],get(gca,'ylim'),-E.angle);
            E.pos = [xy(1,1), xy(1,2) xy(2,1) xy(2,2);];
        end
x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
[x,y] = Data2Norm(x,y);
if length(E.pos) > 7
    xp = E.pos([5 7]);
    yp = E.pos([6 8]);
else
xp = [mean(x) + diff(y)/4 mean(x)-diff(y)/4];
yp = [mean(y) - diff(x)/4 mean(y)+diff(x)/4];
end
if isfield(E,'lh') & ishandle(E.lh(1)) 
    set(E.lh(1),'Position',[x(1) y(1) x(2) y(2)]);
    h = E.lh;
    if length(E.h) > 1 && ishandle(E.h(2))
        set(E.h(2),'Xdata',xp,'Ydata',yp);
    end
else
    hold on;
    h(1) = annotation('line',real(x),real(y),varargin{:});
    %h(2) = plot(xp,yp,'r:');
    hold off;
end

function h= testDrawEllipse(E,varargin)

if E.shape(1) == 1 || E.shape(1) == 2
    h = testDrawLine(E,varargin{:});
    return;
end


function h= DrawEllipse(E,varargin)

if E.shape(1) == 1 || E.shape(1) == 2
    h = DrawLine(E,varargin{:});
    return;
end
if ~isfield(E,'xyr')  || E.down 
E.xyr(1) = mean(E.pos([1 3]));
E.xyr(2) = mean(E.pos([2 4]));
E.xyr(3) = abs(diff(E.pos([1 3])))/2;
E.xyr(4) = abs(diff(E.pos([2 4])))/2;
end
yl = get(gca,'ylim');
xl = get(gca,'xlim');
x = get(gca,'position');
fp = get(gcf,'position');
ar = abs(diff(x([1 3])).* diff(fp([1 3])) ./   (diff(x([2 4])) .*diff(fp([2 4]))));

a = E.xyr(3); %x radius in nomalized units
b = E.xyr(4);
sn = 0;
cn = 1;
x = linspace(0,a);
if isfield(E,'aspectratio') && E.aspectratio > 0
    b  = b ./E. aspectratio;
    y =  sqrt(b.^2 - (x.*b/a).^2);
else
    y =  sqrt(b.^2 - (x.*b/a).^2);
end

sn = sin(E.angle);
cn = cos(E.angle);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+mean(E.xyr(1));
if isfield(E,'aspectratio') && E.aspectratio > 0
    y = yr.*E.aspectratio+mean(E.xyr(2));
else
    y = yr+mean(E.xyr(2));
end    

if isfield(E,'h') & ishandle(E.h)
    set(E.h,'Xdata',real(x),'Ydata',real(y));
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end

function PlotAllCellSpikes(DATA)
    square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];
    markh = NaN;;
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');
    [eid, cid, clid] = FindCell(DATA, cellid);

    SetFigure(DATA,DATA.tag.celllist);
    spkcolors = mycolors('spkcolors');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    AllSpikes = CheckAllSpikes(DATA,eid,cid);
    for j = 1:length(eid)
        SetFigure(DATA,DATA.tag.celllist);
        if ishandle(markh)
            delete(markh);
        end
        markh = plot(square(1,:)+cid(j),square(2,:)+eid(j),'g-','linewidth',2);
        if isempty(AllSpikes)
            SpoolSpikeFile(DATA, eid(j), cid(j));
        else
            colors = spkcolors;
            SetFigure(DATA,DATA.tag.xyplot,'front');
            hold off;
            if clid(j) ~= 1
                np = clid(j)+1;
                colors{np} = spkcolors{2};
                colors{2} = spkcolors{np};
            end
            PlotClusterXY(DATA,Clusters{eid(j)}{cid(j)}, 'colors', colors);
            drawnow;
            stopped = SpoolSpikes(DATA, [eid(j) cid(j)], 'colors', colors);
            if stopped
                return;
            end
        end
    end
    
function DATA = SetCellEntry(DATA, C,  e, p, c, cellid, varargin)    

    savelist = 1;
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'nosave',6)
            savelist = 0;
        end
        j = j+1;
    end
    if isempty(C)
        Clusters = getappdata(DATA.toplevel,'Clusters');
        C = Clusters{e}{p};
    end
    quality = 4;
    if c < 0
        id = find(DATA.CellList(e,p,:) == cellid);
        if length(id) == 1
            quality = 4 + c;
            c = id;
        else
            return;
        end
    end
    if c <= size(DATA.CellList,3)
        oldcell = DATA.CellList(e,p,c);
    else
        oldcell = 0;
    end
    if cellid > 0
        [a,b] = find(squeeze(DATA.CellList(e,:,:)) == cellid);
        DATA.CellList(e,a,b) = 0; %remove this cell from other clusters this expt
    end
    DATA.CellList(e,p,c) = cellid;
    DATA.CellChanges = cat(1,DATA.CellChanges,[e p cellid c now oldcell]);
    if c > 1
        if c > length(C.next)+1
            return;
        end
        C = C.next{c-1};
    end
    if isfield(C,'excludetrialids') && length(C.excludetrialids)
        % find trial #s that match the ids excluded
        DATA.CellDetails.excludetrials{e,p,c} = union(DATA.CellDetails.excludetrials{e,p,c}, C.excludetrialids);
    end
    DATA.CellDetails.Quality(e,p,c) = quality; %default
    if savelist
        SaveCellList(DATA);
    end

function SetCellNumber(a,b, fcn)
DATA = GetDataFromFig(a);
[Clusters, DATA] = CheckClusterLoaded(DATA, 1:length(DATA.exptid));
modifier = 0;
cellid = 0;

    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);

if sum(strcmp(fcn,{'set' 'duplicates'}))
    it = findobj(get(a,'parent'),'Tag','ClusterModifier');
    if ~isempty(it)
        modifier = get(it,'value');
        strs = get(it,'string');
        modstr = deblank(strs(modifier,:));
        if modifier == 4
            DATA.CellDetails.note(DATA.currentpoint(1),DATA.currentpoint(2)) = get(it,'value');
            DATA.CellDetails.notestr = get(it,'string');
        end
        if modifier >  4
        set(it,'value',1);
        end
    end
    it = findobj(get(a,'parent'),'Tag','CellNumberId');
        if isempty(it)
            return;
        end
        cellid = get(it,'value');
        if strcmp(fcn,'duplicates')
            cellid = -cellid;
        end
    if sum(DATA.selectprobe(:)) > 0 &&  ismember(modifier, [1 2 3 4])
        for j = 1:size(DATA.selectprobe,1)
            p = find(DATA.selectprobe(j,:) > 0);
            if length(p) > 1
                errordlg(sprintf('%d probes selected in expt %d',length(p),j));
            elseif length(p)
                DATA = SetCellEntry(DATA, Clusters{j}{p}, j,p, modifier, cellid,'nosave');
            end
        end
        PlotCellList(DATA,'showfig');
        set(DATA.toplevel,'UserData',DATA);
        SaveCellList(DATA);
        return;
    end

    
    if ~isempty(it)
        cellid = get(it,'value');
        DATA.currentcell = cellid;
        if strcmp(modstr,'Quality-poor');
            DATA = SetCellEntry(DATA, Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)}, DATA.currentpoint(1),DATA.currentpoint(2), -2, cellid);
        elseif strcmp(modstr,'Quality-suboptimal');
            DATA = SetCellEntry(DATA, Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)}, DATA.currentpoint(1),DATA.currentpoint(2), -1, cellid);
        elseif modifier == 5
            id = find(DATA.CellList == cellid);
            DATA.CellList(id) = 0;
            DATA.CellDetails.Templates(cellid,:) = DATA.templateid;
            id = find(DATA.usepeaks > 0);
            for j = id;
                oldcell = DATA.CellList(j,DATA.usepeaks(j),DATA.cellcluster);
                DATA.CellList(j,DATA.usepeaks(j),DATA.cellcluster) = cellid;
                DATA.CellDetails.Quality(j,DATA.usepeaks(j),DATA.cellcluster) = 4;
                DATA.CellChanges = cat(1,DATA.CellChanges,[j DATA.usepeaks(j) cellid DATA.cellcluster now oldcell]);
            end
            %        id = find(DATA.usepeaks == 0);
            %       deletes = sum(DATA.CellList(id,:) == cellid);
            SaveCellList(DATA);
        elseif ismember(modifier, [1 2 3 4])
            DATA = SetCellEntry(DATA, Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)}, DATA.currentpoint(1),DATA.currentpoint(2), modifier, cellid);
        end
        set(DATA.toplevel,'UserData',DATA);
    end
    PlotCellList(DATA,'showfig');
elseif strcmp(fcn,'delete')
    it = findobj(get(a,'parent'),'tag','CellCluster');
    cl = get(it,'value');
    it = findobj(get(a,'parent'),'Tag','CellNumberId');
    if length(it) == 1
        cellid = get(it,'value');
    end
    DATA = DeleteCell(DATA, DATA.currentpoint(1),DATA.currentpoint(2),cl);
    [a,b] = find(DATA.selectprobe);
    for j  = 1:length(a)
        id = find(ismember(DATA.CellList(a(j),b(j),:),cellid));
        if length(id) == 1
            DATA.cellcluster = id;
            DATA = DeleteCell(DATA, a(j), b(j), DATA.cellcluster);
        end
    end
    SaveCellList(DATA);
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    PlotCellList(DATA, 'showfig');
    set(DATA.toplevel,'UserData',DATA);
elseif strcmp(fcn,'change')
    DATA.currentcell = get(a,'value');
    CellChanged(DATA);
    cl = find(DATA.CellList(e,p,:) == DATA.currentcell);
    if cl > 0
        it = findobj('Tag','CellCluster');
        set(it,'value',cl);
    end
    set(DATA.toplevel,'UserData',DATA);
end

function SaveCellList(DATA)
    cellfile = [DATA.name '/CellList.mat'];
    CellList = DATA.CellList;
    CellDetails = DATA.CellDetails;
    CellChanges = DATA.CellChanges;
    if isfield(DATA,'fitjumps')
        CellDetails.fitjumps = DATA.fitjumps;
    end
    save(cellfile, 'CellList','CellDetails','CellChanges');


function DATA = PlotCellList(DATA, varargin)
plotmahal = 0;
force = 0;
reload = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'showfig',5) %Force figure creation
        force = 1;
    elseif strncmpi(varargin{j},'reload',5) %reload from disk
        reload = 1;
    end
    j = j+1;
end


if force
    f = SetFigure(DATA, DATA.tag.celllist);
else
    f = findobj('Tag',DATA.tag.celllist,'type','figure');
    if isempty(f)
        return;
    else
        figure(f);
    end
end
if reload
        DATA = LoadCellFile(DATA);
end

if isempty(DATA.CellList)
    return;
end
DATA.fig.celllist = f;
colors = mycolors;
subplot(1,1,1);
set(f,'UserData',DATA.toplevel);
hold off;
set(gca,'position',[0.03  0.1 0.94 0.85]);
if plotmahal
im = zeros(size(DATA.CellList));
id = find(sum(DATA.CellList,3) > 0);
X = squeeze(DATA.mahal(:,:,1));
im(id) = X(id);
end

if isfield(DATA,'nclusters')
    ncls = DATA.nclusters;
    if size(ncls,1) < size(DATA.CellList,1)
        ncls(size(DATA.CellList,1),:) = 1;
    end
else
    ncls = ones(size(DATA.CellList));
end

for j = 1:size(DATA.CellList,1)
    for k = 1:size(DATA.CellList,2)
        id = find(DATA.CellList(j,k,:) > 0);
        if max(id) > ncls(j,k)
            ncls(j,k) = max(id);
        end
    end
end

nc = max(ncls); %max # clusters in any one expt, for each probe
%so nc has length = to number of probes
offset = DATA.clusteroffset;
if plotmahal
    imh = imagesc(im,'buttondownfcn',{@HitImage, 1});
    caxis([0 5]);
else
    for j = 1:size(DATA.CellList,2)*2
        if size(DATA.CellList,3) > offset
        CellIm(:,j) = DATA.CellList(:,ceil(j/2),1+offset);
        else
            CellIm(:,j) = 0;
        end
    end
    id = find(nc >1);
    if size(DATA.CellList,3) > offset+1
    for j = 1:length(id)
        tid = find(ncls(:,id(j)) > 1);
        CellIm(tid,id(j)*2) = DATA.CellList(tid,id(j),offset+2);
    end
    end
        
CellIm(CellIm<0) = -1;
imh = imagesc([0.75 size(DATA.CellList,2)+0.25],[1 size(DATA.CellList,1)],CellIm,'buttondownfcn',{@HitImage, 3});
end
hold on;
showlines = 2;  %may reactivate one day....

%Get cell colors before making context menu
if showlines == 2
    cells = unique(DATA.CellList);
    cells = cells(cells > 0);
    cm = colormap(gca);
    yl = get(gca,'ylim');
    for j = 1:length(cells)
        id = find(DATA.CellList == cells(j));
        [a,y,b] = ind2sub(size(DATA.CellList),id);
        h = text(mean(y)-1+mean(b)/3,yl(1),sprintf('%d',cells(j)));
        set(h,'verticalalignment','bottom');
        cid = round(j .*length(cm)./length(cells));
        set(h,'color',cm(cid,:));
        DATA.cellcolors{cells(j)} = cm(cid,:);
    end
end

cmenu = AddContextMenu(DATA,'cellplot');
set(imh,'uicontextmenu',cmenu);


if DATA.plot.cellgrid > 0 && DATA.markcell.grid
for p = 1:DATA.plot.cellgrid:size(DATA.CellList,2)
    plot([p-0.5 p-0.5],get(gca,'ylim'),'k:','buttondownfcn',{@HitImage, 3});
end
for p = 1:DATA.plot.cellgrid:size(DATA.CellList,1)
    plot(get(gca,'xlim'),[p-0.5 p-0.5],'k:','buttondownfcn',{@HitImage, 3});
end
end
[a, eid] = intersect(DATA.CellDetails.exptids,DATA.exptid);
cells = unique(DATA.CellList(eid,:,1+offset));
cells = cells(cells > 0);


if showlines == 1
    for j = 1:length(cells)
        [x,y] = find(DATA.CellList(eid,:,1+offset) == cells(j));
        for k = 1:length(y)
            if DATA.nclusters(x(k),y(k)) > 1
                plot([y(k)-0.25 y(k)-0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
            else
                plot([y(k) y(k)], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
            end
        end
        [a,b] = min(x);
        h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
        set(h,'color',colors{j},'buttondownfcn',{@HitImage, 3});
    end
    
    if size(DATA.CellList,3) > offset+1
        cells = unique(DATA.CellList(:,:,2+offset));
        cells = cells(cells > 0);
        for j = 1:length(cells)
            [x,y] = find(DATA.CellList(:,:,2+offset) == cells(j));
            for k = 1:length(y)
                plot([y(k)+0.25 y(k)+0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
            end
            [a,b] = min(x);
            h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
            set(h,'color',colors{j},'buttondownfcn',{@HitImage, 3});
        end
    end
else
end




iscellim = sum(DATA.CellList,3) > 0;
if size(DATA.CellList,3) >= 2+offset
    iscellimb = squeeze(DATA.CellList(:,:,2+offset)) > 0;
else
    iscellimb = zeros(size(iscellim));
end
%reduce size in case all clsuters not loaded yet
if isfield(DATA,'dropi')
    iscellim = iscellim(1:size(DATA.dropi,1),:); 
    iscellimb = iscellimb(1:size(DATA.dropi,1),:); 
end
if DATA.markcell.candidates
    [x,y] = find(isnan(DATA.CellList(:,:,1)));
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','g');
    end
    if isfield(DATA,'missing')
    [x,y] = find(DATA.missing);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','y');
    end
    end
end
if DATA.markcell.goodmu
    [x,y] = find(DATA.marked ==4);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','g');
    end
end
if DATA.markcell.goodmu
    [x,y] = find(DATA.marked ==2);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','g');
    end
end
if DATA.markcell.dropi > 0 && isfield(DATA,'dropi')
    DATA.dropi(DATA.dropi == 0) = NaN;
    [x,y] = find(squeeze(DATA.dropi(:,:,1+offset)) < DATA.markcell.dropi & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','m');
    end
    if size(DATA.dropi,3) > 1+offset
    [x,y] = find(squeeze(DATA.dropi(:,:,2+offset)) < DATA.markcell.dropi & iscellimb);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'box2','color','m');
    end
    end
end
if DATA.markcell.mahal > 0
    X = max(DATA.mahal(:,:,[1 3]),[],3);
    [x,y] = find(X < DATA.markcell.mahal & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','r');
    end
end
if DATA.markcell.ellipses > 0
    [x,y] = find(DATA.cutshape ==0 & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','m');
    end
end
if DATA.markcell.tagged > 0 && isfield(DATA,'tagged');
    [x,y] = find(DATA.tagged);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','y');
    end
    for j = 1:length(DATA.Comments)
        DrawBox(DATA.Comments(j).ex,DATA.Comments(j).p,3,'color','g');
    end
end
if DATA.markcell.readmethod > 0
    Clusters = getappdata(DATA.toplevel,'Clusters');
    DATA.readmethod = CellToMat(Clusters,'exptreadmethod');
    [x,y] = find(DATA.readmethod ==0);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','y');
    end
end
if DATA.markcell.quick > 0
    Clusters = getappdata(DATA.toplevel,'Clusters');
    mstatus = CellToMat(Clusters,'manual');
    quickstatus = CellToMat(Clusters,'quick');
    [x,y] = find(quickstatus ==1);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','y');
    end
    [x,y] = find(mstatus == 2);
    for j = 1:length(x)
        DrawBox(x(j),y(j),3,'color','m');
    end
end

if  isfield(DATA,'fitjumps')
    DATA.fitjumps(isnan(DATA.fitjumps)) = 0;
    drift = [0 cumsum(DATA.fitjumps)];
    p = DATA.currentpoint(2)+drift(DATA.currentpoint(1));
    plot(p-drift,1:length(drift),'w-','buttondownfcn',{@HitImage, 3});
end

if DATA.show.ed
    probesep = median(DATA.probesep(DATA.probesep > 0));
    epos = (DATA.electrodedepth - DATA.electrodedepth(DATA.currentpoint(1))) * 1000./probesep;
    p = DATA.currentpoint(2);
    plot(p-epos,1:length(epos),'w-','buttondownfcn',{@HitImage, 3});
end
if DATA.markcell.expnames
    for j = 1:length(DATA.exptid)
        h = text(0.5,j,DATA.expnames{j},'horizontalalignment','left','buttondownfcn',{@HitImage, 3},'color','w');
    end
end
if DATA.plotcells.showmahal > 0
    m = (max(DATA.mahal(:,:,[1 3]),[],3) > DATA.plotcells.showmahal);
    dp = (squeeze(DATA.GaussFitdp(:,:,2)) > DATA.plotcells.showfitdp & sum(DATA.gmfitpos,3) == 2);
    peakok = (abs(DATA.peakdiff) < 1);
    dropok = (squeeze(DATA.dropi(:,:,1)) > 1);
    [aid, bid] = find(m & peakok & dropok);
    for j = 1:length(aid)
        DrawBox(aid(j),bid(j),3,'color','r');
    end
    [cid, did] = find(dp & peakok & dropok);
    for j = 1:length(cid)
        DrawBox(cid(j),did(j),3,'color','g');
    end
    [aid, bid] = find(dp & m & peakok & dropok);
    for j = 1:length(aid)
        DrawBox(aid(j),bid(j),3,'color','w');
    end
end

if sum(DATA.selectprobe(:)) > 0
    [a,b] = find(DATA.selectprobe > 0);
    for j = 1:length(a)
        h = DrawBox(a(j),b(j),3, 'color','w','uicontextmenu',cmenu);
    end
elseif isfield(DATA,'currentpoint')
    h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),3,'color','w','uicontextmenu',cmenu);
end

it = findobj(gcf,'Tag','CellNumberId');
if length(it) == 1
    set(it,'value',DATA.currentcell);
end
MarkExpts(DATA, DATA.markexpts);
set(gca,'uicontextmenu', cmenu);
set(gcf, 'KeyPressFcn',{@KeyPressed,3},'Keyreleasefcn',{@KeyReleased, 3});




function DATA = InitInterface(DATA)

    scrsz = get(0,'Screensize');
    hspace = 0.005;
    cntrl_box = figure('Position', [10 scrsz(4)-310 360 200],...
        'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.name,'menubar','none');
    DATA.toplevel = cntrl_box;
    listh = 0.3;
    lst = uicontrol(gcf, 'Style','listbox','String', DATA.strings,...
            'Callback', {@PlotClusters, 'setentry'},...
            'units','norm', 'Position',[hspace 0 0.7 listh]);
    DATA.lstui = lst;
    nr = 12;
    nc = 4;
    bp = [0.8+hspace listh+0.1./nr 1./(nc)-hspace 1./nr];
    uicontrol(gcf,'style','pop','string','plain|probe|cutspace', ...
        'Callback', {@Update}, 'Tag','ColorScheme',...
        'units', 'norm', 'position',bp,'value',1);
    
    bp = [hspace listh+0.1./nr 1./(nc)-hspace 1./nr];
    bp(1) = hspace;
    bp(3) = 0.1;
    bp(2) = listh+ 0.1./nr;
    bp(3) = 0.05;
    uicontrol(gcf,'style','text','units','norm','string','E','position',bp);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.12;
    uicontrol(gcf,'style','pop','string',num2str([1:50]'), ...
        'Callback', {@SetExpt, 'set'}, 'Tag','ExptList',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.045;
    uicontrol(gcf,'style','pushbutton','string','\/', ...
        'Callback', {@NextButton, 'd'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','pushbutton','string','/\', ...
        'Callback', {@NextButton, 'u'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);

    bp(1) = hspace;
    bp(2) = bp(2)+1./nr;
    bp(3) = 0.05;
    uicontrol(gcf,'style','text','units','norm','string','P','position',bp);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.12;
    uicontrol(gcf,'style','pop','string',num2str([1:24]'), ...
        'Callback', {@SetProbe, 'set'}, 'Tag','ProbeList',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.045;
    uicontrol(gcf,'style','pushbutton','string','<', ...
        'Callback', {@NextButton, 'l'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','pushbutton','string','>', ...
        'Callback', {@NextButton, 'r'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = bp(2)+1./nr;
    bp(1) = hspace;
    bp(3) = 1./nc;
    uicontrol(gcf,'style','pushbutton','string','spool', ...
        'Callback', {@SpoolCurrentSpikes}, 'Tag','SpoolButton',...
        'units', 'norm', 'position',bp,'value',1);


    bp(2) = listh+ 0.1./nr;
    bp(1) = 0.29;
    bp(3)=1./nc;
    uicontrol(gcf,'style','pop','string','quality|times|xcorr|ptscatter|mahaldip|bestspace|bestdip|bmcmahal|cov(shape)|xc(shape)', ...
        'Callback', {@PlotClusters, 'setplot'}, 'Tag','plottype',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','pop','string','mahal|fitdp|dips|mahaln|mahal1-n|mahal+var|man-auto|drop-mahal|Fit-mahal|muamp-spkvar|xcorr|probexcorr|probeco|spkshapecorr|spkshape|spkshapeim||spksize|Evec|BuildTimes|BuildDates|PcGms|Exclusions|Tagged|TriggerFind|CandidateCells|SpkRate|EventRate|shapexc', ...
        'Callback', {@PlotAllClusters}, 'Tag','plotalltype',...
        'units', 'norm', 'position',bp,'value',1);
    hspace = 0.005;
    bp(1) = hspace;
    bp(1) = 0.29;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','spool', ...
        'Callback', {@Update}, 'Tag','SpoolSpikes',...
        'units', 'norm', 'position',bp,'value',DATA.spoolspikes);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','spk xy', ...
        'Callback', {@Update}, 'Tag','SpikeXY',...
        'units', 'norm', 'position',bp,'value',DATA.showspkxy);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','spkmean', ...
        'Callback', {@Update}, 'Tag','SpikeMean',...
        'units', 'norm', 'position',bp,'value',DATA.showspkmean);
    bp(1) = hspace;
    bp(1) = 0.29;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','Histogram', ...
        'Callback', {@Update}, 'Tag','PlotHistogram',...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Show GM', ...
        'Callback', {@Update}, 'Tag','ShowGM',...
        'units', 'norm', 'position',bp,'value',DATA.plot.showgm);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','refit GM', ...
        'Callback', {@Update}, 'Tag','RefitGM',...
        'units', 'norm', 'position',bp,'value',DATA.refitgm);
    bp(1) = bp(1)+bp(3)+hspace;
    if strcmp(DATA.plotexpttype,'none')
        plotexpt = 0;
    else
        plotexpt = 1;
    end

    bp(1) = hspace;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','AllXY', ...
        'Callback', {@Update}, 'Tag','PlotAllXY',...
        'units', 'norm', 'position',bp,'value',DATA.plotallxy);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','xyseq', ...
        'Callback', {@Update}, 'Tag','PlotXYseq',...
        'units', 'norm', 'position',bp,'value',DATA.plotxyseq);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','TrigHist', ...
        'Callback', {@Update}, 'Tag','PlotTrigHist',...
        'units', 'norm', 'position',bp,'value',DATA.plottrighist);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Quickspks', ...
        'Callback', {@Update}, 'Tag','PlotSpks',...
        'units', 'norm', 'position',bp,'value',DATA.plotspks);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','AllVPcs', ...
        'Callback', {@Update}, 'Tag','AllVPcs',...
        'units', 'norm', 'position',bp,'value',DATA.show.allvpcs);

    
    
    bp(1) = hspace;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','text','string','Exclude:', ...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Peak Mismatch', ...
        'Callback', {@Update}, 'Tag','PeakMismatch',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.offpeak);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Peak Match', ...
        'Callback', {@Update}, 'Tag','PeakMatch',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.onpeak);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','NonCell', ...
        'Callback', {@Update}, 'Tag','NonCell',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.noncell);
    
    bp(1) = hspace;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','text','string','Label', ...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Ex#', ...
        'Callback', {@Update}, 'Tag','LabelExptno',...
        'units', 'norm', 'position',bp,'value',DATA.show.exptno);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','Name', ...
        'Callback', {@Update}, 'Tag','LabelExptName',...
        'units', 'norm', 'position',bp,'value',DATA.show.exptname);
    bp(1) = bp(1)+bp(3)+hspace;
    uicontrol(gcf,'style','checkbox','string','elec depth', ...
        'Callback', {@Update}, 'Tag','LabelEd',...
        'units', 'norm', 'position',bp,'value',DATA.show.ed);
    bp(1) = hspace;
    bp(3) = 0.4;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','text','string','Distance Measure/Crit',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.15;
        uicontrol(gcf,'style','pop','string','1D|2D|ND|2Gauss|Dip', ...
        'Callback', {@GuiMenu, 'setmahaltype'}, 'Tag','setmahaltype',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+hspace;
    bp(3) = 0.15;
    uicontrol(gcf,'style','edit','string',num2str(DATA.crit.mahal), ...
        'Callback', {@GuiMenu, 'setmahalcrit'}, 'Tag','setmahalcrit',...
        'units', 'norm', 'position',bp,'value',1);

    hm = uimenu(cntrl_box,'Label','&File','Tag','XClusterMenu');
    uimenu(hm,'Label','&Load All','Callback',{@LoadAll, 'force'});
    uimenu(hm,'Label','&ReLoad New','Callback',{@LoadAll, 'ifnew'});
    uimenu(hm,'Label','Use Current as &Template','Callback',{@PlotClusters, 'followcorr'});
    uimenu(hm,'Label','Check Clusters','Callback',{@GuiMenu, 'checkclusters'});
    uimenu(hm,'Label','Reexamine raw &V','Callback',{@GuiMenu, 'CallAllVPcs'},'accelerator','a');
    uimenu(hm,'Label','Load All FullV','Callback',{@LoadFullVs, 'fullv'});
    uimenu(hm,'Label','&Close','Callback',{@PlotClusters, 'close'});
    hm = uimenu(cntrl_box,'Label','Options','Tag','OptionMenu');
    uimenu(hm,'Label','&Density','Callback',{@OptionMenu, 1});
    uimenu(hm,'Label','2D GM Clustering','Callback',{@OptionMenu, 8});
    uimenu(hm,'Label','Plot Cell','Callback',{@OptionMenu, 2});
    uimenu(hm,'Label','Plot AllXY This Expt','Callback',{@OptionMenu, 4});
    uimenu(hm,'Label','Plot AllXY This Probe','Callback',{@OptionMenu, 11});
    uimenu(hm,'Label','Plot MeanSpike (line) This Expt','Callback',{@OptionMenu, 5});
    uimenu(hm,'Label','Plot MeanSpike(IM) This Expt','Callback',{@OptionMenu, 6});
    uimenu(hm,'Label','Plot Expt All Probes','Callback',{@OptionMenu, 3});
    uimenu(hm,'Label','Exclude Selected Trials','Callback',{@OptionMenu, 7});
    uimenu(hm,'Label','Exclude Trials Cl2','Callback',{@OptionMenu, 'excludecl2'});
    uimenu(hm,'Label','Exclude Trials Cl3','Callback',{@OptionMenu, 'excludecl3'});
    uimenu(hm,'Label','Exclude Trials Cl4','Callback',{@OptionMenu, 'excludecl4'});
    uimenu(hm,'Label','clear exclusions Cl1','Callback',{@OptionMenu, 'usealltrialscl1'});
    uimenu(hm,'Label','clear exclusions Cl2','Callback',{@OptionMenu, 'usealltrialscl2'});
    uimenu(hm,'Label','Use Auto Clusters','Callback',{@OptionMenu, 'useautoclusters'});
    uimenu(hm,'Label','xcorr selected','Callback',{@OptionMenu, 9});
    uimenu(hm,'Label','xcorr all','Callback',{@OptionMenu, 10});
    uimenu(hm,'Label','xcorr cells','Callback',{@OptionMenu, 'xcorrcells'});
    uimenu(hm,'Label','xcorr cell all Expts','Callback',{@OptionMenu, 'xcorrallcells'});
    uimenu(hm,'Label','xcorr all Expts/Probes','Callback',{@OptionMenu, 'xcorrallprobes'});
    uimenu(hm,'Label','Find Candidate Connections','Callback',{@OptionMenu, 'findconnect'});
    uimenu(hm,'Label','Smoothed Image','Callback',{@OptionMenu, 12});
    uimenu(hm,'Label','CellFinder','Callback',{@OptionMenu, 13});
    sm = uimenu(hm,'Label','PlotAll');
    uimenu(sm,'Label','Expts for 1 Probe','Callback',{@OptionMenu, 'setallprobeplot'});
    uimenu(sm,'Label','Expts for 1 Cell','Callback',{@OptionMenu, 'setallcellplot'});
    sm = uimenu(hm,'Label','Layout/Setting');
    uimenu(sm,'Label','Save Layout','Callback',{@OptionMenu, 'savelayout'});
    uimenu(sm,'Label','Save default layout','Callback',{@OptionMenu, 'savedefaultlayout'});
    uimenu(sm,'Label','Load layout','Callback',{@OptionMenu, 'loadlayout'});
    uimenu(sm,'Label','Windows to Front','Callback',{@OptionMenu, 'tofront'});
    uimenu(sm,'Label','Load Settings','Callback',{@OptionMenu, 'loadconfig'});
    uimenu(sm,'Label','Save Settings','Callback',{@OptionMenu, 'saveconfig'});
    uimenu(sm,'Label','Save Default Settings','Callback',{@OptionMenu, 'savedefaultconfig'});
    sm = uimenu(hm,'Label','Settings');
%    uimenu(sm,'Label','Save','Callback',{@OptionMenu, 'saveconfig'});
%    sm = uimenu(hm,'Label','Spooling','Tag','SpoolOptions');
    uimenu(sm,'Label','Used Saved Cluster Codes','Callback',{@OptionMenu, 'usesavedcodes'});
    uimenu(sm,'Label','Show Excluded Trials','Callback',{@OptionMenu, 'showexcluded'});
    uimenu(sm,'Label','Superimpose Cell Meanss','Callback',{@OptionMenu, 'showcellmeans'});
    sm = uimenu(hm,'Label','AllVPcs mode','tag','allvpcsmode');
    uimenu(sm,'Label','from spk file','Callback',{@OptionMenu, 'allvpcsmode'},'tag','fromspikes');
    uimenu(sm,'Label','usecluster','Callback',{@OptionMenu, 'allvpcsmode'},'tag','usecluster');
    uimenu(sm,'Label','reapply','Callback',{@OptionMenu, 'allvpcsmode'},'tag','reapply');

    
    hm = uimenu(cntrl_box,'Label','&Plots','Tag','PlotMenu','accelerator','p');
    sm = uimenu(hm,'Label','Expts','Tag','ExptPlotMenu');
    uimenu(sm,'Label','Mean Tuning','Callback',{@PlotMenu, 'expt', 'means'});
    uimenu(sm,'Label','Mean (RC too)','Callback',{@PlotMenu, 'expt', 'rcmeans'});
    uimenu(sm,'Label','Trial Counts','Callback',{@PlotMenu, 'expt', 'trialcounts'});
    uimenu(sm,'Label','Combine','Callback',{@PlotMenu, 'expt', 'combine'});
    uimenu(sm,'Label','None','Callback',{@PlotMenu, 'expt', 'none'});
    sm = uimenu(hm,'Label','&Cells','Tag','CellPlotMenu');
    AddPlotMenu(sm, 'cells');
    sm = uimenu(hm,'Label','&Probes','Tag','ProbePlotMenu');
    AddPlotMenu(sm, 'probes');

    sm = uimenu(hm,'Label','MeanSpike','Tag','ProbePlotMenu');
    uimenu(sm,'Label','Image','Callback',{@PlotMenu, 'mean', 'imageonly'});
    uimenu(sm,'Label','Image+lines','Callback',{@PlotMenu, 'mean', 'DpLines'});
    uimenu(sm,'Label','Mean Lines','Callback',{@PlotMenu, 'mean', 'meanlines'});
    uimenu(sm,'Label','Lines (with dprime)','Callback',{@PlotMenu, 'mean', 'lineonly'});
    uimenu(sm,'Label','Image','Callback',{@PlotMenu, 'mean', 'MeanIm'});
    uimenu(sm,'Label','MU Image','Callback',{@PlotMenu, 'mean', 'TwoIm'});
    sm = uimenu(hm,'Label','CrossCorr (&x)','Tag','xCorrPlotMenu');
    uimenu(sm,'Label','Cells all &Expts','Callback',{@OptionMenu, 'xcorrallcells'});
    uimenu(sm,'Label','Cells in &current Expt','Callback',{@OptionMenu, 'xcorrcells'},'accelerator','C');
    uimenu(sm,'Label','&Selected Probes','Callback',{@OptionMenu, 9},'accelerator','c');
    uimenu(sm,'Label','All c&lusters this expt','Callback',{@OptionMenu, 'allclusterxcorr'},'accelerator','e');
    uimenu(sm,'Label','All Expts/&Probes','Callback',{@OptionMenu, 'xcorrallprobes'});
    uimenu(sm,'Label','xcorr &all','Callback',{@OptionMenu, 10});
    uimenu(sm,'Label','Find Candidate Connections','Callback',{@OptionMenu, 'findconnect'});
    uimenu(sm,'Label','&Recalculate All Cells','Callback',{@OptionMenu, 'recalcxcorrcell'});
    sm = uimenu(hm,'Label','&Other','Tag','OtherPlotMenu');
    uimenu(sm,'Label','Shape/Size','Callback',{@PlotMenu, 'cells' 'shapesize'});

    uimenu(hm,'Label','Current Pop &Summary','Callback',{@OptionMenu, 'popplotl'});

    sm = uimenu(cntrl_box,'Label','&Trials','Tag','TrialMenu');
    uimenu(sm,'Label','&Exclude Selected Trials','Callback',{@OptionMenu, 7});
    uimenu(sm,'Label','Exclude Trials Cl&2','Callback',{@OptionMenu, 'excludecl2'});
    uimenu(sm,'Label','Exclude Trials Cl&3','Callback',{@OptionMenu, 'excludecl3'});
    uimenu(sm,'Label','Exclude Trials Cl&4','Callback',{@OptionMenu, 'excludecl4'});
    uimenu(sm,'Label','&clear exclusions Cl1','Callback',{@OptionMenu, 'usealltrialscl1'});
    uimenu(sm,'Label','clear exclusions Cl2','Callback',{@OptionMenu, 'usealltrialscl2'});
    uimenu(sm,'Label','&Show Excluded Trials','Callback',{@OptionMenu, 'showexcluded'});
    sm = uimenu(cntrl_box,'Label','Tag','Tag','TagMenu');
    uimenu(sm,'Label','Possible Cell','Callback',{@TagMenu, '?cell'});
    uimenu(sm,'Label','AllVPcs - Threshold','Callback',{@TagMenu, 'threshold'});
    uimenu(sm,'Label','AllVPcs  > 1 cell','Callback',{@TagMenu, 'morecells'});
    uimenu(sm,'Label','AllVPCs - improve','Callback',{@TagMenu, 'improve'});
    uimenu(sm,'Label','AllVPCs - error','Callback',{@TagMenu, 'error'});
    %    hm = uimenu(sm,'Label','Comment','Tag','TagMenu');
    uimenu(sm,'Label','Add Comment Manually','Callback',{@TagMenu, 'comment'},'accelerator','M');
%    uimenu(sm,'Label','Add Comment (test)','Callback',{@TagMenu, 'comment'},'accelerator','M');
    uimenu(sm,'Label','Comment: Poor Stability','Callback',{@TagMenu, 'poor stability'});
    uimenu(sm,'Label','Comment: Poor Isolation','Callback',{@TagMenu, 'poor isolation'});
    uimenu(sm,'Label','Comment: Dropping Spikes','Callback',{@TagMenu, 'dropping spikes'});
%    uimenu(sm,'Label','Comment: Same As another probe','Callback',{@TagMenu, 'repeat probe'});
    uimenu(sm,'Label','Clear Selected','Callback',{@TagMenu, 'clear'});
    uimenu(sm,'Label','Print Tags','Callback',{@TagMenu, 'print'});
    DATA.tagstrings = {'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear'};

    DATA.fig.top = DATA.toplevel;
    setappdata(0,'control_is_down',0);
    setappdata(0,'alt_is_down',0);
    [DATA.fig.clusters, isnew] = SetFigure(DATA,DATA.tag.clusters);
    tmp.parentfig = DATA.toplevel;
    set(DATA.fig.clusters,'Name','Clusters','UserData',tmp);
    set(cntrl_box,'UserData',DATA);
    
function AddPlotMenu(sm, type)
        
    if strcmp(type,'probes')
    uimenu(sm,'Label','&Probe All Spks +xy','Callback',{@PlotMenu, 'probes', 'probeall'});
    uimenu(sm,'Label','&Expt All Spks +xy','Callback',{@PlotMenu, 'probes', 'exptall'},'accelerator','e');
    uimenu(sm,'Label','Spool All Probes','Callback',{@PlotMenu, 'probes', 'spoolall'});
    uimenu(sm,'Label','Spool All Probes All Expts','Callback',{@PlotMenu, 'probes', 'spooleverything'});
    uimenu(sm,'Label','Spool Probes With Cells','Callback',{@PlotMenu, 'probes', 'spoolcells'});
    uimenu(sm,'Label','A&ll Spks This Expt','Callback',{@PlotMenu, 'probes', 'allspks'});
    uimenu(sm,'Label','&All Spks This Probe','Callback',{@PlotMenu, 'probes', 'allprobespks'});
    uimenu(sm,'Label','All Spks Template Peaks','Callback',{@PlotMenu, 'probes', 'alltemplatespks'});
    uimenu(sm,'Label','&Selected Spks','Callback',{@PlotMenu, 'probes', 'selectspks'});
    uimenu(sm,'Label','All Cell Spks This Expt','Callback',{@PlotMenu, 'probes', 'allcellspks'});
    uimenu(sm,'Label','Selected XY','Callback',{@PlotMenu, 'probes', 'SelectXY'});
    uimenu(sm,'Label','All &XY This Expt','Callback',{@PlotMenu, 'probes', 'AllXY'},'accelerator','x');
    uimenu(sm,'Label','All X&Y This Probe','Callback',{@PlotMenu, 'probes', 'AllprobeXY'},'accelerator','y');
    uimenu(sm,'Label','All XY Template Pekas','Callback',{@PlotMenu, 'probes', 'AllTemplateXY'});
    uimenu(sm,'Label','Mean Spikes','Callback',{@PlotMenu, 'probes', 'AllMean'});
    uimenu(sm,'Label','All Probes Mean Spike','Callback',{@PlotMenu, 'probes', 'AllProbeMean'});
    uimenu(sm,'Label','Mean Spikes Image','Callback',{@PlotMenu, 'probes', 'AllMeanIm'});
    uimenu(sm,'Label','Mean Spike All Probes and Expts','Callback',{@PlotMenu, 'probes', 'AllExptMeanIm'});
    uimenu(sm,'Label','Spike dprime All Probes and Expts','Callback',{@PlotMenu, 'probes', 'AllExptdpIm'});
    elseif strcmp(type,'cells')
        uimenu(sm,'Label','All&XY','Callback',{@PlotMenu, 'cells', 'AllXY'});
        uimenu(sm,'Label','All&Mean','Callback',{@PlotMenu, 'cells', 'AllMean'});
        uimenu(sm,'Label','AllMean&Im','Callback',{@PlotMenu, 'cells', 'AllMeanIm'});
        uimenu(sm,'Label','All&TrigHist','Callback',{@PlotMenu, 'cells', 'TrigHist'});
        uimenu(sm,'Label','Spool All - This expt','Callback',{@PlotMenu, 'cells', 'spoolall'});
        uimenu(sm,'Label','Spool All - All expts','Callback',{@PlotMenu, 'cells', 'spooleverycell'});
        uimenu(sm,'Label','Spool selected, all &expts','Callback',{@PlotMenu, 'cells', 'spoolone'});
        uimenu(sm,'Label','AllSpks Cell','Callback',{@PlotMenu, 'cells', 'allspks'});
        uimenu(sm,'Label','All Cells Spks Expt','Callback',{@PlotMenu, 'probes', 'cellspks'});
        uimenu(sm,'Label','Spool All Cells','Callback',{@PlotMenu, 'cells', 'spoolall'});
        uimenu(sm,'Label','Rate &Sequence','Callback',{@PlotMenu, 'cells', 'rateseqone'});
        uimenu(sm,'Label','&Rate Sequence All Cells','Callback',{@PlotMenu, 'cells', 'rateseqall'});
    end
    
function AddExptList(hm, callbacklabel, DATA)
    
    KnownExpts = { 'image.&OT' 'image.or'; ...
        '&grating.OXM' 'grating.orXme'; ...
        '&rls.OTRC' 'rls.orRC'; ...
        '&image.ORBW' 'image.orXobP'; ...
        '&Plaid psych' 'checker.pRP'; ...
        'checker.OT' 'checker.or'; ...
        'rds.&AC' 'rds.dxXce'; ...
        'rds.&FACRC' 'rds.dxXceXFrRC'; ...
        'image.SFOB' 'image.sfXob'; ...
        'image.SZOB' 'image.szXob'};
    h = uimenu(hm,'Label','None','Callback',{@PlotMenu, callbacklabel 'none'});
    if strncmp(callbacklabel,'select',5)
        set(h,'Label','All');
    end
    if isfield(DATA,'expnames')
        names = unique(DATA.expnames);
        for j = 1:length(names)
            nexpts = sum(strcmp(names{j},DATA.expnames));
            id = find(strcmp(names{j},KnownExpts(:,2)));
            if length(id) == 1
                uimenu(hm,'Label',sprintf('%s(%d)',KnownExpts{id,1},nexpts),'Callback',{@PlotMenu, callbacklabel, names{j}});
            else
                uimenu(hm,'Label',sprintf('%s(%d)',names{j},nexpts),'Callback',{@PlotMenu, callbacklabel, names{j}});
            end
        end
    end
    uimenu(hm,'Label','relist','Callback',{@PlotMenu, callbacklabel 'relist'});

    
function DATA = LoadXcorrFiles(DATA)    
    for j = 1:length(DATA.strings)
    end
    
function DATA = LoadFullVs(a,b, type, varargin)
    
    DATA = GetDataFromFig(a);
    AllFullV = LoadAllFullV(DATA.name,'verbose');
    setappdata(DATA.toplevel,'AllFullV',AllFullV);
    DATA.allvpcsmode = 'reapply';
    set(DATA.toplevel,'UserData',DATA);
    
function DATA = LoadAll(a,b, type, varargin)
    useman = 1;
    loadexpts = 1;
    AllExpts = {};
    DATA = GetDataFromFig(a);
    msgmode = 0;
    expts = [];
    exargs = {};
    ClusterInfo = {};
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'expts',5)
            j = j+1;
            expts = varargin{j};
        elseif strncmpi(varargin{j},'nowarn',5)
            msgmode = 2;
        elseif strncmpi(varargin{j},'usealltrials',7)
            exargs = {exargs{:} varargin{j}};
        end
        j = j+1;
    end
    
    
    
    if isfield(DATA,'toplevel')
    Clusters = getappdata(DATA.toplevel,'Clusters');
    if isempty(Clusters)
        Clusters = {};
    end
    end
    if strcmp(b,'auto')
        useman = 0;
    end
    np = 0;
    ManExpt = [];
    NewExpt = [];
    AutoExpt = [];
    FullVData = {};
    bid = strmatch('Copy of',DATA.strings);
    id = setdiff(1:length(DATA.strings),bid);
    DATA.strings = DATA.strings(id);
    
    for j = 1:length(DATA.strings)
        if regexp(DATA.strings{j},'Expt[0-9]*a')
            offset = 0.1;
        else
            offset = 0;
        end
        if strfind(DATA.strings{j},'AutoClusterTimes.mat')
            k = sscanf(DATA.strings{j},'Expt%d');
            AutoExpt(j) = k+offset;
            ManExpt(j) = 0;
        elseif strfind(DATA.strings{j},'NewClusterTimes.mat')
            k = sscanf(DATA.strings{j},'Expt%d');
            NewExpt(j) = k+offset;
        elseif useman && strfind(DATA.strings{j},'ClusterTimes.mat')
            k = sscanf(DATA.strings{j},'Expt%d');
            ManExpt(j) = k+offset;
            AutoExpt(j) = 0;
        end
    end
    
    outname = [DATA.name '/ExptList.mat'];
    if exist(outname,'file')
        Ex = load(outname);
        if isfield(Ex.ExptList,'usealltrials') && Ex.ExptList.usealltrials > 0
            fprintf('Setting usealltrials\n');
            DATA.usealltrials = 1;
        end
    end

    %need to think about cases where there was just the autocut and a 
    %manual cut is created
    if strcmp(type,'ifnew') %Already set
        AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
        GMfits = getappdata(DATA.toplevel,'GMfits');
        ClusterInfo = getappdata(DATA.toplevel,'ClusterInfo');

        if isempty(GMfits)
            GMfits = {};
        end
            
        
        d = dir(DATA.name);
        mid = [];
        for j  = 1:length(d)
            if strfind(d(j).name,'AutoClusterTimes.mat')
            elseif strfind(d(j).name,'ClusterTimes.mat') & isempty(strfind(d(j).name,'OldClusterTimes.mat'))
                id = strmatch(d(j).name,DATA.strings);
                if isempty(id) || d(j).datenum > DATA.lastcheck%new file
                    if length(id) == 1
                        DATA.strings{id} = d(j).name;
                        AutoClusters{id} = Clusters{id};
                    else
                        id = length(DATA.strings)+1;
                    end
                    exptno = sscanf(d(j).name,'Expt%d');
                    eid = find(DATA.exptid == exptno);
                    
                    fprintf('Loading %s (expt %d)\n',d(j).name, eid);
                    [C, FullV, details] = LoadCluster(DATA.name, exptno,'rawxy','alltimes');
                    ClusterInfo{eid}.loadname = details.loadname;
                    ClusterInfo{eid}.loadtime = details.loadtime;
                    for c = 1:length(C)
                        if size(AllSpikes,1) >= eid && size(AllSpikes,2) >= c && ~isempty(AllSpikes{eid,c})
                            if C{c}.savetime(1) > AllSpikes{eid,c}.Header.ctime %spikes are out of date
                                AllSpikes{eid,c} = [];
                            end
                        end
                    end
                    [Clusters{eid}, GMfits{eid}] = CondenseCluster(C);
                end
            end
        end
        cid = 1:length(DATA.strings);
        setappdata(DATA.toplevel,'AllSpikes',AllSpikes);
        setappdata(DATA.toplevel,'Clusters',Clusters);
        setappdata(DATA.toplevel,'ClusterInfo',ClusterInfo);
        setappdata(DATA.toplevel,'GMfits',GMfits);
        loadexpts = 0;
    elseif useman == 0
        DATA.strings = DATA.strings(find(AutoExpt));
        cid = 1:length(DATA.strings);
        mid = [];
    elseif length(ManExpt)
        %cid find expts where AutoExpt is 0 (= ManExpt set), or AutoExpt is
        %not in ManExpt;
        cid = find(~ismember(AutoExpt,unique(ManExpt(ManExpt > 0))));
        cid = setdiff(cid, find(NewExpt > 0));
        if ~isempty(expts)
            [gex, cid] = intersect(ManExpt(cid),expts);
        end
        
        mid = find(ManExpt > 0);
    else
        cid = 1:length(AutoExpt);
        mid = [];
    end
    
    DATA.lastcheck = now;
    reloaded = [];
    ts = now;
    AllExpts = {};
    DATA.DataType = 'Default';
    DATA.idsorted = 0;
    if loadexpts > 0 && exist(DATA.exptname,'file')
        [DATA, AllExpts] = LoadExpts(DATA,exargs{:});
    end
    if strncmp(DATA.DataType,'Grid', 4)
        DATA.nprobeplot = 0; %Don't plot adjoining probes
        DATA.probesperfile = 1;
    end
    if isempty(cid) %no clustertimes data yet
        DATA = AddError(DATA,'No ClusterTimes Files in %s (Use ''loadauto'' to get AutoClusterTimes)',DATA.datadir);
        return;
    end
    monk = GetMonkeyName(DATA.name);
    for j = 1:length(cid)
        k = cid(j);
        if strfind(DATA.strings{k},'Times.mat') && isempty(strfind(DATA.strings{k},'NewClusterTimes'))
            d = dir([DATA.name '/' DATA.strings{k}]);
            if d.datenum > DATA.dates(k) || strcmp(type,'force')
                DATA.dates(j) = d.datenum;
                if ismember(k,mid) && DATA.loadautoclusters
                    autocl = 0;
                    [AutoClusters{j}, autodetails] = LoadClusters([DATA.name '/' strrep(DATA.strings{k},'Cluster','AutoCluster')]);
                else %if no manual cluster, use the AutoCluster as main....
                    [AutoClusters{j}, autodetails] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                    autocl = 1;
                end
                ClusterInfo{j} = autodetails;
 
                if strcmp(type,'ifnew') %Only Change modified Clusters, so hat don't lose bits from AutoClusters
                    details.exptno = Clusters{k}{1}.exptno;
                    if 0 %%dont need this any more
                        fprintf('Loading %s\n',DATA.strings{k});
                        [C, details] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                        for k = 1:length(C)
                            %if only cluster 2 has changed, ctime for cluster 1 does not. But savetime does.
                            if isfield(C{k},'xy') & (C{k}.ctime > Clusters{j}{k}.ctime || C{k}.savetime(1) > Clusters{j}{k}.savetime(1))
                                Clusters{j}{k} = C{k};
                                reloaded(j,k)=1;
                            end
                        end
                    end
                else
                    fprintf('Loading %s\n',DATA.strings{k});
                    [Clusters{j}, details] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                    if isfield(details,'FullVData')
                        FullVData{j} = details.FullVData;
                    else
                        fprintf('No FullV Data for %s\n',DATA.strings{j});
                    end
                end
                if isfield(details,'loadtime') %Not using AutoClusters
                    ClusterInfo{j} = details;
                end
                np = max([np length(Clusters{j})]);
                exptno(j) = details.exptno;
                for k = 1:length(Clusters{j})
                    if ~isfield(Clusters{j}{k},'mahal')
                        s = sprintf('Cluster %d Ex %d missing',k,j);
                        if msgmode == 0
                            s = questdlg(sprintf('Cluster %d Ex %d missing',k,j),'ClusterError','OK','Ignore Others','Ignore Others');
                            if strcmp(s,'Ignore Others')
                                msgmode = 1;
                            end
                        else
                            if msgmode == 1
                            errordlg(sprintf('Cluster %d Ex %d missing',k,j),'ClusterError','modal');
                            end
                            mycprintf('errors','Cluster %d Ex %d missing\n',k,exptno(j));
                        end
                        if size(AutoClusters,1) >= j && size(AutoClusters,2) >= k && ~isempty(AutoClusters{j}{k})
                            Clusters{j}{k} = AutoClusters{j}{k};
                            Clusters{j}{k}.auto = 1;
                        else
                            Clusters{j}{k}.auto = NaN;
                        end
                    elseif length(AutoClusters{j}) >= k && AutoClusters{j}{k}.savetime(1) > Clusters{j}{k}.savetime(1) && Clusters{j}{k}.auto == 1
                        Clusters{j}{k} = AutoClusters{j}{k};
                    end
                    if autocl
                        Clusters{j}{k}.auto = 1;
                    end
                    if length(AutoClusters{j}) < k
                         A = [];
                    else
                        A = AutoClusters{j}{k};
                    end
                        
                    if isfield(Clusters{j}{k},'nspks')  %absent if bad probe/no cut
                    if k >1 && length(Clusters{j}{k}.probe) ==1 && Clusters{j}{k}.probe < k 
                        Clusters{j}{k}.probe = k;
                    end
                    if ~isfield(Clusters{j}{k},'clst')
                        DATA = AddError(DATA,'Missing clst for %d,%d',exptno(j),k);
                    end
                    if ~isfield(Clusters{j}{k},'xy') && isfield(A,'xy')
                        Clusters{j}{k}.xy = A.xy;
                        if ~isfield(Clusters{j}{k},'clst') 
                            Clusters{j}{k}.clst = A.clst;
                            Clusters{j}{k}.times = A.times;
                        end
                        %next test should not be true, but is sometimes. Eg
                        %211 Ex13 P24.  clst should not be in Clusters.
                        if length(Clusters{j}{k}.times) < length(Clusters{j}{k}.clst) && Clusters{j}{k}.auto == 1
                            Clusters{j}{k}.clst = A.clst;
                            Clusters{j}{k}.times = A.times;
                        end
                    end
                    end
                    if ~isfield(Clusters{j}{k},'exptno') || isempty(Clusters{j}{k}.exptno)
                        Clusters{j}{k}.exptno = details.exptno;
                    end
                    if ~isfield(Clusters{j}{k},'auto')
                        Clusters{j}{k}.auto = 1;
                    end
                    if ~isfield(Clusters{j}{k},'sign')
                        Clusters{j}{k}.sign= 0;
                    end
                    if isfield(Clusters{j}{k},'t')
                        fprintf('E%.1fP%d has t\n',Clusters{j}{k}.exptno,k);
                    end
                    if ~isfield(Clusters{j}{k},'exptno')
                        Clusters{j}{k}.exptno = exptno(j);
                    end
                end
                smrname = regexprep(DATA.name,['/M([0-9]*)(.*)'],['$0/' monk 'M$1']);
                exfile = [smrname '.' sprintf('%d',floor(exptno(j))) 'idx.mat'];
                if exist(exfile,'file') && loadexpts > 0
                    precount = 0;
                    fprintf('Loading %s\n',exfile);
                    [Trials, Expts] = APlaySpkFile(exfile,'noerrs',exargs{:});
                    if isempty(Expts)
                        Expt.Header.exptno = exptno(j);
                    else
                        Expt = FillTrials(Expts{1},'ed');
                        Expt.Header.errs = Trials.errs;
                        Expt.Header.exptno = exptno(j);
                        DATA.electrodedepth(j) = mean([Expt.Trials.ed]);
                    end
                    if isfield(Trials.Comments,'Peninfo') & isfield(Trials.Comments.Peninfo,'probesep')
                        DATA.probesep(j) = Trials.Comments.Peninfo.probesep;
                    end
                    DATA.trialids{j} = [Expt.Trials.id];
                    Expt.Header.expname = Expt2Name(Expt);
                    AllExpts{j,1} = Expt;
                end
            else
                DATA.electrodedepth(j) = NaN;
                exptno(j) = Clusters{j}{1}.exptno;
            end
            end
    end

    fprintf('Load Time %.1f\n',mytoc(ts));
       
    for j = 1:length(Clusters)
        Clusters{j} = FixClusters(Clusters{j});
    end
    DATA = LoadComments(DATA, DATA.name);
   if strcmp(type,'ifnew')     && length(reloaded) && DATA.spikesloaded
      [a,b] = find(reloaded);
      LoadSelectedSpikes(DATA,a,b);
   end
    [exptno, id] = sort(exptno);
    DATA.strings = DATA.strings(cid(id));
    DATA.dates = DATA.dates(id);
    DATA.exptid = exptno;
    if ~isfield(DATA,'trialids') && isdir(DATA.name) %No expts loaded. ? Dir name but only one expt
        [a,b] = fileparts(DATA.name);
        monkey = GetMonkeyName(DATA.name);
        matfile = [DATA.name '/' monkey b '.mat'];
        if exist(matfile,'file')
            DATA.exptname = matfile;
            [DATA, AllExpts] = LoadExpts(DATA);
        end
    end
    
    if DATA.idsorted == 0
        DATA.trialids = DATA.trialids(id);
    end
    if length(FullVData) == length(id)
        FullVData = FullVData(id);
    end
    Clusters = Clusters(id);
    ClusterInfo = ClusterInfo(id);
    if length(AllExpts)
        if isempty(DATA.exptname) %%expts in separate files. Named like Clusterfiles
            Expts = AllExpts(id,:);
        else  %only use Expts that have Clusters loaded
            Expts = AllExpts(exptno,:);
        end
        
        for j = 1:size(Expts,1)
            DATA.expnames{j} = Expt2Name(Expts{j,1},'addsuff');
            DATA.electrodedepth(j) = mean(GetEval(Expts{j,1},'ed'));
        end
        setappdata(DATA.toplevel,'Expts',Expts);
        ExptList.expnames = DATA.expnames;
    else
        ExptList.expnames = [];
    end
    ExptList.exptid = DATA.exptid;
    ExptList.usealltrials = DATA.usealltrials;
    if ~isempty(ExptList.expnames) %can happen with relist
        outname = [DATA.name '/ExptList.mat'];
        save(outname,'ExptList');
    end
   if ~strcmp(type,'ifnew')   
    DATA.nprobes = np;
    DATA.voffset = [1:DATA.nprobes] .*2;

   end

    set(DATA.toplevel,'UserData',DATA);
    if DATA.useautoclusters
        setappdata(DATA.toplevel,'AutoClusters',AutoClusters);
    end
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    it = findobj(DATA.toplevel,'Tag','ExptList');
    if length(it) == 1
        set(it, 'string',num2str(DATA.exptid'));
    end
    set(DATA.toplevel,'UserData',DATA);
    [Clusters, GMfits] = CondenseClusters(Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
    setappdata(DATA.toplevel,'ClusterInfo',ClusterInfo);
    if isappdata(DATA.toplevel,'GMfits')
        X = getappdata(DATA.toplevel,'GMfits');
        for j = 1:length(GMfits)
            if ~isempty(GMfits{j})
                X{j} = GMfits{j};
            end
        end
        setappdata(DATA.toplevel,'GMfits',X);
    else
        setappdata(DATA.toplevel,'GMfits',GMfits);
    end
    setappdata(DATA.toplevel,'FullVData',FullVData);
    DATA = LoadCellFile(DATA);
    DATA = LoadExtra(DATA,0);
    DATA = CalcDistances(DATA);
    
    it = findobj(DATA.toplevel,'Tag','ProbeList');
    set(it,'string',num2str([1:DATA.nprobes]'))
    fprintf('Done\n');

    DATA = CheckExpts(DATA,'read');
    DATA = CheckExpts(DATA,'errs');
    if strcmp(type,'loadspikes')
        fprintf('Loading Spikes...');
        LoadAllSpikes(DATA,'quiet');
        DATA.spikesloaded = 1;
        set(DATA.toplevel,'UserData',DATA);
        fprintf('\n');
    end
    CheckClusters(DATA, 'fittimes');
    
    function [Expt, DATA]  = LoadExpt(DATA, e)    
        mnk = GetMonkeyName(DATA.name);
        smrname = regexprep(DATA.name,[mnk '/([MG][0-9]*)(.*)'],['$0/' mnk '$1']);
%        smrname = regexprep(smrname,'jbe/G([0-9]*)(.*)','$0/jbeG$1');
        exfile = [smrname '.' sprintf('%d',floor(e)) 'idx.mat'];
        if exist(exfile,'file')
            precount = 0;
            fprintf('Loading %s\n',exfile);
            [Trials, Expts] = APlaySpkFile(exfile,'noerrs');
            if isempty(Expts)
                Expt.Header.exptno = DATA.exptid(e);

            else
                if length(Expts) > 1
                    for j = 1:length(Expts)
                        nt(j) = length(Expts{j}.Trials);
                    end
                    fprintf([sprintf(' %d',nt) '\n']);
                    [a,b] = max(nt);
                    Expt = FillTrials(Expts{b},'ed');
                else
                    Expt = FillTrials(Expts{1},'ed');
                end
                Expt.Header.errs = Trials.errs;
                Expt.Header.exptno = e;
                DATA.electrodedepth(e) = mean([Expt.Trials.ed]);
            end
            if isfield(Trials.Comments,'Peninfo') & isfield(Trials.Comments.Peninfo,'probesep')
                DATA.probesep(e) = Trials.Comments.Peninfo.probesep;
            end
            DATA.expnames{e} = Expt2Name(Expt,'addsuff');
            DATA.exped(e) = mean(GetEval(Expt,'ed'));
            DATA.trialids{e} = [Expt.Trials.id];
        end

function [DATA, Expts]  = LoadExpts(DATA, varargin)
    
    
    if DATA.usealltrials
        varargin = {varargin{:} 'usealltrials'};
    end
    if exist(DATA.exptname)
        fprintf('Loading %s\n',DATA.exptname);
        [Trials, Expts] = APlaySpkFile(DATA.exptname,'noerrs','nospikes',varargin{:});
        for j = 1:length(Expts)
            Expts{j} = FillTrials(Expts{j},'ed');
            Expts{j}.Header.exptno = j;
            DATA.electrodedepth(j) = mean([Expts{j}.Trials.ed]);
            DATA.trialids{j} = [Expts{j}.Trials.id];
            if isfield(Expts{j}.Header,'DataType')
                DATA.DataType = Expts{j}.Header.DataType;
            end
        end
        if isfield(Trials.Comments,'Peninfo') & isfield(Trials.Comments.Peninfo,'probesep')
            DATA.probesep(j) = Trials.Comments.Peninfo.probesep;
        end
        Expts = Expts';
        DATA.idsorted = 1;
    else
        exptno = DATA.exptid;
        for j = 1:length(exptno)
            e = floor(exptno(j));
            [Expts{e,1} , DATA] = LoadExpt(DATA, e);
        end
        Expts = SetExptTimeOffset(Expts);
    end
    setappdata(DATA.toplevel,'Expts',Expts);
        

function Clusters = FixClusters(Clusters)
    for k = 1:length(Clusters)
        if ~isfield(Clusters{k},'next')
            Clusters{k}.next = {};
        elseif ~iscell(Clusters{k}.next)
            last = rmfields(Clusters{k}.next,'next');
            Clusters{k} = rmfield(Clusters{k},'next');
            Clusters{k}.next{1} = last;
        end
    end
    
function DATA = CheckExpts(DATA, type)
    Expts = getappdata(DATA.toplevel,'Expts');
 for j = 1:length(Expts)
     if strcmp(type,'errs')  &&  isfield(Expts{j}.Header,'errs')
         for k = 1:length(Expts{j}.Header.errs)
             DATA = AddError(DATA,'%d: %s\n',j,Expts{j}.Header.errs{k});
         end
     elseif strcmp(type,'method')
         DATA = AddError(DATA,'Ex %d Method %d\n',j,Expts{j}.Header.ReadMethod);
     end
 end
 
function res = CheckSpikeFiles(DATA, type)
    
Clusters = getappdata(DATA.toplevel,'Clusters');
AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
res.prefix = DATA.name;
nact = 1;
if strcmp(type,'spkid')
Im = zeros(size(AllSpikes));

for j = 1:size(AllSpikes,1)
    for k = 1:size(AllSpikes,2)
        S = AllSpikes{j,k};
        C = Clusters{j}{k};
        a = length(C.clst);
        b = length(S.times);
        if b ~= a
            fprintf('Ex %d(%.0f), Probe %d %d Spikes but Cluster has %d\n',j,C.exptno,k,b,a);
            if C.auto == 1
                res.acts(nact).args = {'reclassify' 'autocut' 'savespikes'};
            else
                res.acts(nact).args = {'reclassify' 'savespikes'};
            end
            res.acts(nact).exptid = j;
            res.acts(nact).exptno = floor(C.exptno);
            res.acts(nact).probes = k;
            if rem(C.exptno,1) > 0
            res.acts(nact).name = sprintf('Expt%.0faFullV.mat',C.exptno);
            else
            res.acts(nact).name = sprintf('Expt%.0fFullV.mat',C.exptno);
            end
            
            nact = nact+1;
            Im(j,k) = 1;
        end
            
    end
end
SetFigure(DATA,DATA.tag.clusters);
imagesc(Im);
elseif strcmp(type,'build')
     for j = 1:length(Clusters)
         for k = 1:length(Clusters{j})
             C = Clusters{j}{k};
             Clusters{j}{k}.probe = k;
             need(j,k) = CheckSpikeFile(DATA,Clusters{j}{k});
             if need(j,k)
                 if C.auto == 1
                     res.acts(nact).args = {'reclassify' 'autocut' 'savespikes'};
                 else
                     res.acts(nact).args = {'reclassify' 'savespikes'};
                 end
                 res.acts(nact).exptid = j;
                 res.acts(nact).exptno = floor(C.exptno);
                 res.acts(nact).probes = k;
                 if rem(C.exptno,1) > 0
                     res.acts(nact).name = sprintf('Expt%.0faFullV.mat',C.exptno);
                 else
                     res.acts(nact).name = sprintf('Expt%.0fFullV.mat',C.exptno);
                 end
                 nact = nact+1;
             end
         end
     end
     fprintf('%d Files need building\n',sum(need(:)));
end

function DATA = CheckClusterLineSign(DATA)
    
Clusters = getappdata(DATA.toplevel,'Clusters');
Expts = getappdata(DATA.toplevel,'Expts');
sgn = CellToMat(Clusters,'sign');
shp = CellToMat(Clusters,'shape');
angles = CellToMat(Clusters,'angle');

if isfield(DATA,'tagged')
    oldsum = sum(DATA.tagged(:) ==2);
    DATA.tagged(DATA.tagged ==2) = 0;
else
    oldsum = 0;
end
p = find(sum(sgn < 0) & sum(sgn > 0));
for j =1:length(p)
    [a,b] = Counts(sgn(:,p(j)));
    if a(b == -1) > a(b == 1)
        id = find(sgn(:,p(j)) == 1);
        gid = find(sgn(:,p(j)) == -1);
    else
        id = find(sgn(:,p(j)) == -1);
        gid = find(sgn(:,p(j)) == 1);
    end
    shape = prctile(shp(gid,p(j)),50);
    angle = prctile(angles(gid,p(j)),50);
    for k = 1:length(id)
        if shp(id(k),p(j)) == shape && cos(angle - angles(id(k),p(j))) > 0
            mycprintf('blue','E%d P%d Cluster Line Inverted\n',id(k),p(j));
            DATA.tagged(id(k),p(j)) = 2;
        end
    end     
end

if isfield(DATA,'tagged') && sum(DATA.tagged(:) == 2) || oldsum > 0
    DATA.markcell.tagged = 1;
    PlotCellList(DATA);
    set(DATA.toplevel,'UserData',DATA);
end

 function CheckClusters(DATA, type)
    
Clusters = getappdata(DATA.toplevel,'Clusters');
Expts = getappdata(DATA.toplevel,'Expts');
if iscell(type)
    for j = 1:length(type)
        CheckClusters(DATA,type{j});
    end
    return;
end

 for j = 1:length(Clusters)
     nt = length(Expts{j,1}.Trials);
     for k = 1:length(Clusters{j})
         C = Clusters{j}{k};
         if strcmp(type,'exclusions') && isfield(C,'excludetrialids')
             if isfield(C,'restricttimerange')
                 fprintf('Row%dP%d %d/%d excluded Trials (%.1f-%.1f)\n',j,k,length(C.excludetrialids),nt,C.restricttimerange(1),C.restricttimerange(2));
             else
                 fprintf('Row%dP%d %d/%d excluded Trials\n',j,k,length(C.excludetrialids),nt);
             end
         elseif strcmp(type,'nclusters')
             if isfield(C,'next') & length(C.next) > 1
                 fprintf('Row%dP%d %d clusters\n',j,k,length(C.next)+1);
             end
         elseif strcmp(type,'fitspace')
             CheckFitDim(C);
         elseif strcmp(type,'fittimes')
             if isfield(C,'savetime') && C.savetime(end)-C.savetime(1) < -0.1
                 cprintf('blue','Expt%d(Row%d)P%d Fit is old\n',DATA.exptid(j),j,k);
             end
         end
     end
 end

function DATA = ConvertExclusion(DATA)
    
    
    Expts = getappdata(DATA.toplevel,'Expts');
    if isempty(Expts) && size(DATA.CellDetails.excludetrials,1) < size(DATA.CellList,1)
        DATA.CellDetails.excludetrials{size(DATA.CellList,1),DATA.nprobes,6} = [];
        return;
    end
    if size(DATA.CellDetails.excludetrials,1) < length(Expts)
        DATA.CellDetails.excludetrials{length(Expts),DATA.nprobes,6} = [];
    end
    if isfield(DATA.CellDetails,'version') & DATA.CellDetails.version > 1.1
        CheckExclusion(DATA);
        return;
    end
    XC = DATA.CellDetails.excludetrials;
    for j = 1:size(DATA.CellDetails.excludetrials,1)
        for k = 1:size(DATA.CellDetails.excludetrials,2)
            for c = 1:size(DATA.CellDetails.excludetrials,3)
                X = DATA.CellDetails.excludetrials{j,k,c};
                if length(X)
                    if max(X) > length(Expts{k,j}.Trials)
                        fprintf('Too Many Excluded E%dP%dC%d: %d vs %d\n',j,k,c,length(X),length(Expts{k,j}.Trials));
                        X = X(X < length(Expts{k,j}.Trials));
                        ids = [Expts{k,j}.Trials(X).id];
                        XC{k,j,c} = ids;
                    else
                        ids = [Expts{k,j}.Trials(X).id];
                        XC{k,j,c} = ids;
                    end
                end
            end
        end
    end
    n =1;
    DATA.CellDetails.version = DATA.version;
    if size(DATA.CellDetails.excludetrials,1)
    DATA.CellDetails.excludetrials = XC;
    end
%    set(DATA.toplevel,'UserData',DATA);
    

function DATA = CheckExclusion(DATA)
    
    nerr = 0;
    Expts = getappdata(DATA.toplevel,'Expts');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    if isempty(Expts)
        return;
    end
    for j = 1:size(DATA.CellDetails.excludetrials,1)
        for k = 1:size(DATA.CellDetails.excludetrials,2)
            for c = 1:size(DATA.CellDetails.excludetrials,3)
%                X = DATA.CellDetails.excludetrials{j,k,c};
                X =[]; % don't have AllExpts loaded any more
                if length(X)
                    id = find(ismember([Expts{j,k}.Trials.id],X));
                    if length(id) < length(X)
                        fprintf('Excluded E%dP%dC%d: %d/%d\n',j,k,c,length(id),length(X));
                        nerr = nerr+1;
                    end
                end
            end
        end
    end
    n =1;
    fprintf('Checking Exclusion lists: %d errors\n',nerr);
    cells = unique(DATA.CellList);
    cells = cells(cells > 0);
    
    for j = 1:length(cells)
        [eid, pid, cid] = FindCell(DATA,cells(j));
        [eid, id] = intersect(eid, DATA.exptid);
        pid = pid(id);
        cid = cid(id);
        for k = 1:length(eid)
            C = Clusters{eid(k)}{pid(k)};
            if sum(C.clst == cid(k)+1) == 0
                fprintf('No Spikes in Cluster for Cell %d in E%dP%dC%d\n',cells(j),eid(k),pid(k),cid(k));
            end
            nc = sum(DATA.CellList(eid(k),pid(k),:) == cells(j));
            if nc > 1
                fprintf('cell %d defined twice in E%dP%d\n',cells(j),eid(k),pid(k));
            end
        end
    end
%    set(DATA.toplevel,'UserData',DATA);

function name = FileName(DATA, ex, probe, type)
    
   smrname = regexprep(DATA.name,'lem/M([0-9]*)','$0/lemM$1');
   if strmatch(type,'Expt')
       name = [smrname '.' num2str(ex) 'idx.mat'];
   elseif strmatch(type,'FullV')
       if DATA.probesperfile == 1
           name = [DATA.name '/Expt' num2str(ex) '.p' num2str(probe) 'FullV.mat'];
       else
           name = [DATA.name '/Expt' num2str(ex) 'FullV.mat'];
       end
   end

    
function DATA = CompareShapes(DATA, type)    
    spk = 1;
    ids = [];
    if ~isfield(DATA,'nprobes')
        DATA.nprobes = 24;
    end
    SetFigure(DATA,DATA.tag.clusters);
    clf;
    if strmatch(type,{'exptlines' 'exptimage'})
        DATA = PlotShapes(DATA, type);
    else
        Clusters = getappdata(DATA.toplevel,'Clusters');
        for j = 1:DATA.nprobes
            subplot(4,6,j);
            if type == 1
                CompareShape(DATA, Clusters, j);
            elseif type == 2
                PlotShape(DATA,j);
            end
            drawnow;
        end
    end


function SpikeDistances(DATA, type)
spk = 1;
ids = [];
SetFigure(DATA, DATA.tag.clusters,'front');
clf;
for j = 1:length(DATA.strings)
    [sizes(j,:,:), ds(j,:,:), quality(j,:)] = SpikeDistance(DATA, j);
end
dval = unique(ds(:));
for j = 1:length(dval)
    id = find(ds == dval(j));
    ms(j) = mean(sizes(id));
end
hold off;
plot(dval, ms,'o-');

for q = [-2 2 4 6]
[rows, cols]= find(quality.* sign(q) > q);
for j = 1:length(rows)
    d(j,:) = ds(rows(j),cols(j),:);
    s(j,:) = sizes(rows(j),cols(j),:);
end
for j = 1:length(dval)
    id = find(d == dval(j));
    ms(j) = mean(s(id));
end
hold on;
plot(dval, ms,'ro-');
end


function [sizes, d, quality] = SpikeDistance(DATA, eid)
spk = 1;
ids = [];
SetFigure(DATA, DATA.tag.clusters,'front');
clf;
C = Clusters{eid};
for j = 1:length(C)
    sizes(j,:) = std(C{j}.MeanSpike.ms')./std(C{j}.MeanSpike.ms(j,:));
    d(j,:) = abs([1:length(C)]-j);
    quality(j) = C{j}.mahal(1);
end

function CompareProbes(DATA, type)
spk = 1;
ids = [];
DATA.nprobes = 24;
SetFigure(DATA, DATA.tag.clusters,'front');
clf;
for j = 1:length(DATA.strings)
    x = CompareProbesShape(DATA, j, type);
    if length(x.shapexc) == DATA.nprobes
    C(j,:,:) = x.shapexc;
    xca(j,:) = x.xca;
    xcb(j,:) = x.xcb;
    end
    drawnow;
end
SetFigure(DATA, DATA.tag.all,'front');
subplot(1,2,1);
Cmean = squeeze(mean(C));
imagesc(Cmean);
subplot(1,2,2);
hold off;
plot(mean(xca));
hold on;
plot(mean(xcb),'r');
set(gca,'xlim',[1 DATA.nprobes]);
[E,V] = eig(Cmean);
SetFigure(DATA,DATA.tag.eigenvectors,'front')
hold off; 
colors = mycolors;
for j = 0:5
    plot(E(:,end-j),'color',colors{j+1});
    hold on;

end

function xcs = ShiftXcorr(allshape, a, b, npts)
    nprobes = round(size(allshape,2)./npts);
    for k = -7:7
        if k >= 0
            apts = 1:((nprobes-k) .* npts);
            bpts = 1+(k .* npts):size(allshape,2);
        else
            bpts = 1:((nprobes+k) .* npts);
            apts = 1-(k .* npts):size(allshape,2);
        end
        xc = corrcoef(allshape(a,apts),allshape(b,bpts));
        xcs(k+8) = xc(1,2);
    end

function DATA = PlotShapes(DATA, type)
    
   [Clusters, DATA] = CheckClusterLoaded(DATA, 1:length(DATA.exptid));
   if strmatch(type,{'exptlines' 'exptimage'})
       for j = 1:length(DATA.exptid)   nx = length(DATA.exptid);

           for k = 1:DATA.nprobes;
               spt(j,k)  =  Clusters{j}{k}.spts(1);
               espt(j,k)  =  Clusters{j}{k}.spts(end);
           end
       end
       minspt = min(spt(:));
       fnpts = 1+max(espt(:))-minspt;  %total # pts needed to cover all spike ranges
       for j = 1:length(DATA.exptid)
           nv = 1;
           fnv = 1;
           for k = 1:DATA.nprobes;
               V = Clusters{j}{k}.MeanSpike.ms;
               allshape(j,nv:nv+size(V,2)-1) = V(k,:);
               df = spt(j,k)-minspt;
               fixshape(j,[fnv:fnv+size(V,2)-1]+df) = V(k,:);
               spts = Clusters{j}{k}.spts -  minspt+1;
%               fixshape(j,k,spts) = V(k,:);
               nv = nv+size(V,2);
               fnv = fnv+fnpts;
               npts = size(V,2);
               nvpts(j,k) = npts;
           end
       end
       aspt = prctile(spt(:),50);
       SetFigure(DATA,DATA.tag.clusters);
       fixed(1) = 0;
       if strcmp(type,'exptimage')
%build an image showing spike shape for each probe/expt
%each row is an expt. Spikes from each probe are concatenated
%calculate correleation coefficients wtih shifts between adjacent rows
%(xcs)
           for j = 2:length(DATA.exptid)
               xcs(j-1,:) = ShiftXcorr(fixshape, j, j-1,fnpts);
               if max(xcs(j-1,:)) < 0.
                   fixshape(j,:) = fixshape(j-1,:);
                   fixed(j) = 1;
               end
           end
           imagesc(xcs);           
%calculate correleation coefficients wtih shifts between all row pairs
%and find teh shfit that maximizes this (mxcs). 
           for j = 2:length(DATA.exptid)
           for k = 1:j
               x = ShiftXcorr(fixshape, j, k,fnpts);
               [a, b] =max(x);
                mxcs(j,k) = b-8;
                mxval(j,k) = a;
           end
           end
           tvals = [];
           jump = 1
           for j = 2:nx
               newvals = mxcs(j,jump:j-1);
               tvals = [tvals newvals];
               tmean(j) = sum(newvals) /sqrt(length(newvals));
               if sum(newvals) > sqrt(length(newvals));
                   jump = j;
                   jumps(j) = mean(newvals);
                   tvals = [];
               end
           end
           jumps(nx) = 0;
           id = find(jumps > 0);
           id = [1 id nx];
           imagesc(mxcs);
           last = 1;
           for j = 2:length(id)-1
               js = mxcs(id(j):id(j+1),id(j-1):id(j)-1);
               jump(j-1) = mean(js(:));
               for k = 3:j
                   js = mxcs(id(j):id(j+1),id(k-2):id(k-1)-1);
                   xc = mxval(id(j):id(j+1),id(k-2):id(k-1)-1);
                   alljump(j-1,k-1) = mean(js(xc > 0.7));
               end
               last = id(j);
           end
           for j = find(fixed)
               xcs(j,:) = ShiftXcorr(fixshape,j+1, j, npts);
           end
           [a,b] = max(xcs');
           c = cumsum(b-8);
           c = cumsum(jumps);
           DATA.driftguess = cumsum(jumps);
           DATA.shiftmatrix(1,:,:) = mxcs;
           DATA.shiftmatrix(2,:,:) = mxval;
           set(DATA.toplevel,'UserData',DATA);
           M = mxcs;
           if sum(mxval(:)>0.85) > prod(size(M))./4 %half are zero
               mcrit = 0.85;
           else
               fprintf('!!!Warning. Spike shape cross correlations low. Drift Estimate may be poor\n');
               mcrit = prctile(mxval(:),75);  %half are zeros anyway, so this is top 50%
           end
           M(mxval(:) < mcrit) = NaN;
% M is a matrix of displacement estimates between all expt pairs
% Now find the sigle drift estimate that best fits this matrix
           [P, details] = FitDriftMatrix(M,'maxiter', 10000);
           DATA.fitjumps = details.jumps;
           c = cumsum(details.jumps);
           SaveExtras(DATA);
           imagesc(fixshape);
           hold on;
           xoff = 910;
           plot(xoff -fnpts.*c,1:length(c));
           b = diff(b);
       elseif strcmp(type,'exptimage')
           imagesc(allshape);
       else
           clf;
           colors = mycolors;
           for j = 1:length(DATA.exptid)
               plot(allshape(j,:)+j*0.5,'color',colors{1+mod(j-1,5)});
               hold on;
           end
       end
   else
       spk = 1;
       ids = [];
       SetFigure(DATA,DATA.tag.clusters);
       clf;
       for j = 1:DATA.nprobes
           subplot(4,6,j);
           Shape(DATA,j);
           drawnow;
       end
   end

    function PlotShape(DATA, spk)
   [Clusters, DATA] = CheckClusterLoaded(DATA, 1:length(DATA.exptid));
  for j = 1:length(DATA.strings)
      if length(Clusters{j}) >= spk
          nv = size(Clusters{j}{spk}.MeanSpike.ms,2);
          shape(1:nv,j) = Clusters{j}{spk}.MeanSpike.ms(spk,:);
          ci(j) = Clusters{j}{spk}.mahal(1);
      end
  end
  a = minmax(shape(:));
  b = minmax(ci);
  shape(nv+1,1:length(ci)) = (ci .* diff(a)/diff(b))+a(1);
  imagesc(shape);

function xc = ShapeCorr(P,Q)
    xc = corrcoef(P.MeanSpike.ms(:),Q.MeanSpike.ms(:));
    xc = xc(1,2);
            
function result = CompareProbesShape(DATA, ex, type)  

C = Clusters{ex};
for j = 1:length(C)
    for k = 1:j
        if type == 2
            xc = corrcoef(C{j}.MeanSpike.ms(:),C{k}.MeanSpike.ms(:));
            shapes(j,k) = xc(1,2);
        else
            shapes(j,k) = sum(C{j}.MeanSpike.ms(:).*C{k}.MeanSpike.ms(:));
        end
    shapes(k,j) = shapes(j,k);
    end
end
subplot(1,2,1);
imagesc(shapes);
for j = 1:size(shapes,1)-1
    a(j)= shapes(j,j+1);
end
for j = 2:size(shapes,1)-1
    b(j)= shapes(j-1,j+1);
end
subplot(1,2,2);
hold off; 
plot(a);
hold on;
plot(b,'r');

result.shapexc = shapes;
result.xca = a;
result.xcb = b;

function CellFinderTest(DATA, type)

    
X = sum(DATA.CellList,3);
X = isnan(X) + X > 0;
[a,b] = find(X > 0);
for j  = 1:length(a)
    DATA = TrackCluster(DATA,b(j),a(j),'UsePts',X);
    cellps(j,:) = DATA.usepeaks;
end
plot(cellps);
MakeCellId(cellps);

function DATA = FindMissing(DATA)

    for e = 1:length(DATA.expnames)
        for p = 1:DATA.nprobes;
            DATA.missing(e,p) = MissedCell(DATA,[e p]);
        end
    end
    
    
function DATA =  FillCellList(DATA, mode)

    Clusters = getappdata(DATA.toplevel,'Clusters')
    if strcmp(mode,'frommark')
        for j = 1:length(Clusters)
            for k = 1:length(Clusters{j})
                C = Clusters{j}{k};
                if isfield(C,'marked')
                    marked(j,k) = C.marked;
                else
                    marked(j,k) = 1;
                end
                if isfield(C,'exptno')
                    exptids(j) = C.exptno;
                end
            end
        end
        if size(marked,1) > 10
            crit = 5;
        else 
            crit = size(marked,1)/2;
        end
        id = find(sum(marked==2) >= crit);
        for j = 1:length(id)
            gid = find(marked(:,id(j)) ==2);
            DATA.CellList(gid, id(j),1) = j;
        end
        if ~isfield(DATA.CellDetails,'exptids')
            DATA.CellDetails.exptids = exptids;
        end
        PlotCellList(DATA);
        s = questdlg('Use This CellList?','Cell List','OK','No','No');
        if strcmp(s,'OK')
             SaveCellList(DATA);
            set(DATA.toplevel,'UserData',DATA);
        end

    end
    
function CellFinder(DATA)
        
mid = DATA.mahaltype;
%First make a smoothed map of mahal distances, and find local maxima
if DATA.mahaltype == 5
    X = squeeze(DATA.mahal(:,:,1));
else
    X = squeeze(DATA.mahal(:,:,mid));
end

[a,b,G] = Gauss2D(0.8,[-3:1:3]);
Y = conv2(X,G,'same');
xm = diff(sign(diff(Y)));
ym = diff(sign(diff(Y,1,2)),1,2);
[a,b] = find(xm(:,2:end-1) < 0 & ym(2:end-1,:) < 0);
mxs = xm(:,2:end-1) + ym(2:end-1,:);

%add in any squares where mahal distance > 3
[c, d] = find(DATA.mahal(:,:,1) > 3 | DATA.mahal(:,:,2) > 3);
a = cat(1,a+1 ,c);
b = cat(1, b+1, d);

for j  = 1:length(a)
    DATA = TrackCluster(DATA,b(j),a(j));
    cellps(j,:) = DATA.usepeaks;
end
plot(cellps);
MakeCellId(cellps);


function [CellId, details] = MakeCellId(cellps)
    
 %cellps = Nstarts X N expts array, each cell being the best matching probe
 %for that expt and that starting point. 
for j = 1:size(cellps,2)
    cim(:,j) = hist(cellps(:,j),[0:24]);
end
 imagesc(cim(2:end,:));
%cmid is frequency with which probes are selected, for each expt

id = find(cim > 0);
[a,cid] = sort(cim(id),'descend');
cid = id(cid); %nonzero elements in descending order
[a,b] = ind2sub(size(cim),cid);

CellId = zeros(size(cim,2),size(cim,1)-1);
CellCounts = CellId;
%b is list of expts
%a is list of probes +1 - if a == 1 means probe was 0  
nc = 0;
for j = 1:length(a)
    if a(j) > 1
%find all runs that have matching probe at this expt
    id = find(cellps(:,b(j)) == a(j)-1);
%
c = hist(cellps(id,:),[0:24]);
[d,e] = max(c(2:end,:));
did = find(d>0);
    x = sub2ind(size(CellId),did,e(did));
    eid = find(d(did) > CellCounts(x));
    x = x(eid);
    if length(x) > 0
%check to see if this is already defined    
    if sum(CellId(x)) == 0
        oldcell = 0;
    else
        oldcell =  mode(CellId(x(CellId(x) > 0)));
        if sum(CellId(x) > 0 & CellId(x) ~= oldcell) > sum(CellId(x) == oldcell)/3
            oldcell = 0;
        end
    end
    if oldcell == 0 
        nc = nc+1;
        CellId(x) = nc;
    else
        CellId(x) = oldcell;
    end
    CellCounts(x) = d(did(eid));
    end
    end
end
for j = 1:size(CellId,1)
    [a,b] = Counts(CellId(j,:));
    id = find(a >1 & b > 0);
    for k = 1:length(id)
        aid = find(CellId(j,:) == b(id(k)));
        [a,c] = max(CellCounts(CellId(j,aid)));
        CellId(j,aid) = 0;
        CellId(j,aid(c)) = b(id(k));
    end
end
hold off;
imagesc(CellId);
fprintf('%d Cells\n',length(unique(CellId(:)))-1);
    
function DATA = TrackCluster(DATA, spk, ex, varargin)
    j = 1;
    usegrid = [];
    maxdistance = 4;
    minmahal = 2;
    while j <= length(varargin)
        if strncmpi(varargin{j},'usepts',6)
            j = j+1;
            usegrid = varargin{j};
            minmahal = 0;
        end
        j = j+1;
    end
Clusters = getappdata(DATA.toplevel,'Clusters');
src = Clusters{ex}{spk};
SetFigure(DATA,DATA.tag.templatesrc);
PlotMeanSpike(Clusters{ex}{spk},spk,1);
SetFigure(DATA,DATA.tag.clusters);

if length(usegrid)
    [a,b] = find(usegrid > 0);
    for j = 1:length(a)
          xc(a(j),b(j)) = CalcTemplateXcorr(src, Clusters{a(j)}{b(j)});
    end
else
    for j = 1:length(Clusters);
      for k = 1:length(Clusters{j});
          xc(j,k) = CalcTemplateXcorr(src, Clusters{j}{k});
      end
    end
end
hold off;
imagesc(xc,'buttondownfcn',{@HitImage, 1});
set(gcf,'KeyPressFcn',{@KeyPressed, 2}, 'Keyreleasefcn',{@KeyReleased, 3});
set(gca,'UserData',DATA.toplevel);
hold on;

if isfield(DATA,'probesep')
%calculate electode pos relative to refernence exp, in multiples of probe
%spacing.
    probesep = median(DATA.probesep);
    epos = (DATA.electrodedepth - DATA.electrodedepth(ex)) * 1000./probesep;
else
    probesep = 0;
    epos = zeros(size(DATA.electrodedepth));
end

xid = spk-2:spk+2;
xid = xid(xid > 0 & xid <= size(xc,2));
for j = ex:size(xc,1)
    if isempty(xid)
        peaks(j) = NaN;
    else
    [rmax(j), peaks(j)] = max(xc(j,xid));
    peaks(j) = peaks(j)+xid(1)-1;
    end
    if j < size(xc,1) && peaks(j) > minmahal
    xid = peaks(j)-2:peaks(j)+2;
    xid = xid(xid > 0 & xid <= size(xc,2));
    xid = xid(abs(xid - spk - epos(j+1)) < maxdistance);
    end
end
xid = spk-2:spk+2;
xid = xid(xid > 0 & xid <= size(xc,2));
for j = ex:-1:1
    [rmax(j), peaks(j)] = max(xc(j,xid));
    peaks(j) = peaks(j)+xid(1)-1;
    xid = peaks(j)-2:peaks(j)+2;
    if j > 1
        ediff = epos(j-1)-epos(j);
        xid = xid(abs(xid - spk - epos(j-1)) < maxdistance);
        if isempty(xid)
            xid = peaks(j)-2:peaks(j)+2;
            if abs(ediff) > 2
%                xid = xid - round(ediff);
            end
        end
    end
    xid = xid(xid > 0 & xid <= size(xc,2));    
end
%[rmax,peaks] = max(xc');
plot(spk-epos,1:size(xc,1),'w-')
plot(peaks,1:size(xc,1),'buttondownfcn',{@HitImage, 1});
square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];

%use 1D or 2D mahal, whichever is better
mahald = max(DATA.mahal(:,:,[1 3]),[],3);

for j = 1:length(peaks)
    if isnan(peaks(j))
        DATA.usepeaks(j) = 0;
    elseif mahald(j,peaks(j)) > DATA.crit.mahal && xc(j,peaks(j)) > DATA.xccrit
        plot(square(1,:)+peaks(j),square(2,:)+j,'w-','buttondownfcn',{@HitImage, 1},'linewidth',2);
        DATA.usepeaks(j) = peaks(j);
    else
        DATA.usepeaks(j) = -peaks(j); %best score, but not cell
    end
end

if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'Templates')
    for j = 1:size(DATA.CellDetails.Templates,1)
        if ex == DATA.CellDetails.Templates(j,2) && spk == DATA.CellDetails.Templates(j,1)
            cellid = j;
            for k = 1:size(DATA.CellList,1)
                p = find(DATA.CellList(k,:,1) == cellid);
                if length(p) && p ~= DATA.usepeaks(k)
                    plot(square(1,:)+p,square(2,:)+k,'w--','buttondownfcn',{@HitImage, 1},'linewidth',2);
                end
            end
        end
    end
end
DATA.steptype = 2;
DATA.markcc = NaN;
DATA.templatexc = xc;
title(sprintf('Tempate Ex %d(%d) P%d',ex,Clusters{ex}{spk}.exptno,spk));
set(DATA.toplevel,'UserData',DATA);


function h = DrawBoxes(DATA, imtype, varargin)

   [e,p] = find(DATA.selectprobe);
   for j = 1:length(e)
       h(j) = DrawBox(e(j),p(j),imtype, varargin{:});
   end
   
function h = DrawBox(ex, p, imtype, varargin)

    lx = -0.5;
    j = 1;
    passon = {};
    while j <= length(varargin)
        if strcmp(varargin{j},'box2')
            lx = 0;
        else
            passon = {passon{:} varargin{j}};
        end
        j = j+1;
    end
    
    if length(ex) > 1 || length(p) > 1
        square = [min(p)+lx max(p)+0.5 max(p)+0.5 min(p)+lx min(p)+lx; min(ex)-0.5 min(ex)-0.5 max(ex)+0.5 max(ex)+0.5 min(ex)-0.5];
    else
        square = [lx 0.5 0.5 lx lx; -0.5 -0.5 0.5 0.5 -0.5];
        square = [square(1,:)+p; square(2,:)+ex];
    end
    h = plot(square(1,:),square(2,:),'k-','buttondownfcn',{@HitImage, imtype},'linewidth',2);
    if length(passon)
        set(h,passon{:});
    end
        
        
        
function xc = CalcTemplateXcorr(tmpl, C)
src = tmpl.MeanSpike.ms;
b = C.MeanSpike.ms;
[x, ai, bi] = intersect(tmpl.spts,C.spts);
chspk = tmpl.probe-1:tmpl.probe+1;
bspk = chspk + C.probe -tmpl.probe;
probes = find(bspk > 0 & bspk <= 24 & chspk > 0 & chspk <= 24);
b = b(bspk(probes),bi);
a = src(chspk(probes),ai);
c = corrcoef(a(:),b(:));
xc = c(1,2);

 function EstimateDrift(DATA, where)

     Clusters = getappdata(DATA.toplevel,'Clusters');
     for e = 1:length(Clusters)-1
         for p = 3:length(Clusters{e})-2
             for j = -1:1;
                 a = Clusters{e}{p}.MeanSpike.ms(p-1:p+1,:);
                 b = Clusters{e+1}{p+j}.MeanSpike.ms(p+j-1:p+j+1,:);
                 l = min([size(a,2) size(b,2)]);
                 a = a(:,1:l);
                 b = b(:,1:l);
                 c = corrcoef(a(:),b(:));
                 xc(e,p,j+2) = c(1,2);
             end
         end
     end
     peaks = squeeze(sum(xc,2));
     [a,b] = max(peaks');
     dpos = cumsum(b-2);
     if strcmp(where,'cell')
         SetFigure(DATA,DATA.tag.celllist)
         plot(DATA.nprobes/2+dpos,[1:length(Clusters)-1]+0.5,'w:');
     end

function CompareShape(DATA, Clusters, spk)
ids = [];
  for j = 1:length(DATA.strings)
      if length(Clusters{j}) >= spk && ~isempty(Clusters{j}{spk}) && size(Clusters{j}{spk}.MeanSpike.ms,1) == 24 ...
              && GoodCluster(Clusters{j}{spk});
          ids = [ids j];
      end
  end
  
  for j = 1:length(ids)
      for k = 1:length(ids)
          a = Clusters{ids(j)}{spk}.MeanSpike.ms;
          b = Clusters{ids(k)}{spk}.MeanSpike.ms;
          l = min([size(a,2) size(b,2)]);
          a = a(:,1:l);
          b = b(:,1:l);
          c = corrcoef(a(:),b(:));
          xc(j,k) = c(1,2);
      end
  end
imagesc(xc);
[E,V] = eig(xc);
allid = [];
allxc = [];
for j = 1:size(E,1)
    p = E(:,j);
    [a,b] = max(abs(p));
    id = find(xc(b,:) > 0.8);
    allid = [allid id];
    allxc = [allxc xc(b,:)];
end


function DATA = MarkCurrentCluster(DATA)
if strmatch(DATA.plot.alltype,{'template' 'mahal+var'})
    SetFigure(DATA, DATA.tag.all);
    square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];
    hold on;
    if ishandle(DATA.markh)
        delete(DATA.markh);
    end
    cmenu = AddContextMenu(DATA,'cellplot');
    DATA.markh = plot(square(1,:)+DATA.currentpoint(2),square(2,:)+DATA.currentpoint(1),'k-','linewidth',2);
    set(DATA.markh,'buttondownfcn',{@HitImage,1},'uicontextmenu',cmenu);
else 
    DATA.markh = [];
end

function DATA = ClearSelections(DATA, force, setcurrent)
%wipe the list of selected probes
%Unless force is set, this may depend on GUI settings/state
    if sum(setcurrent) > 0
        DATA.selectprobe(DATA.currentpoint(1),DATA.currentpoint(2)) = 1;
        if length(setcurrent) == 2
            DATA.currentpoint = setcurrent;
        end
    end
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);

function DATA = SetProbe(mn,b, dir)
    DATA = GetDataFromFig(mn);
    nex = length(DATA.strings);
    
    p = get(mn,'value');
    DATA.currentpoint(2) = p;
    ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2),'oneprobe');
    if DATA.plotallxy
        PlotAllProbeXY(DATA,'oneprobe');
    end
    if DATA.plotspks
        DATA = ClearSelections(DATA,0,0);
        PlotAllProbe(DATA,'allprobespks');
    end
    DATA = MarkCurrentCluster(DATA);
set(DATA.toplevel,'UserData',DATA);

function DATA = SetExpt(mn,b,  fcn)
    DATA = GetDataFromFig(mn);
    nex = length(DATA.strings);
    e = get(mn,'value');
    DATA.currentpoint(1) = e;
    ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2),'oneprobe');
    PrintComments(DATA,e,1:DATA.nprobes);
    if DATA.plotallxy
        PlotAllProbeXY(DATA,'expt',DATA.currentpoint(1));
    end
    DATA = ClearSelections(DATA,0,0);
    if DATA.spoolspikes
        PlotAllProbe(DATA,'spoolall');
    elseif DATA.plotspks
        PlotAllProbe(DATA,'allspks');
    end
    str = sprintf('Ex %.1f: %s',DATA.exptid(DATA.currentpoint(1)),DATA.expnames{DATA.currentpoint(1)});
    fprintf('%s\n',str);
    SetFigure(DATA, DATA.tag.all);
    title(str);
    Expts = getappdata(DATA.toplevel,'Expts');
    DATA.Expt = Expts{e,1};
    SetTrialList(DATA,{}, DATA.currenttrial);
    stopped = 0;
     DATA = MarkCurrentCluster(DATA);

set(DATA.toplevel,'UserData',DATA);


function DATA = NextButton(mn,b, dir)
    DATA = GetDataFromFig(mn);
   nex = length(DATA.strings);
if dir == 'r' %step 1 rihgt
    if DATA.currentpoint(2) < DATA.nprobes
    DATA.currentpoint(2) = DATA.currentpoint(2)+1;
    else
        DATA.currentpoint(2) = 1;
        DATA.currentpoint(1) = DATA.currentpoint(1)+1;
    end
    DATA = ClearSelections(DATA,0,0);
elseif dir == 'l' %step 1 left
    if DATA.currentpoint(2) > 1
    DATA.currentpoint(2) = DATA.currentpoint(2)-1;
    else
        DATA.currentpoint(2) = DATA.nprobes;
        if DATA.currentpoint(1) > 1
        DATA.currentpoint(1) = DATA.currentpoint(1)-1;
        end
    end
    DATA = ClearSelections(DATA,0,0);
elseif dir == 'd' %step 1 down - if looking at Templatesscores, step trough these
    if isfield(DATA,'usepeaks') && DATA.steptype ==2
        if DATA.currentpoint(1) < nex
            j = DATA.currentpoint(1)+1;
            id = find(DATA.usepeaks(j:end) > 0);
            if ~isempty(id)
                j = DATA.currentpoint(1)+id(1);
                DATA.currentpoint(1) = j;
                DATA.currentpoint(2) = DATA.usepeaks(j);
            end
        else
            id = find(DATA.usepeaks> 0);
            if ~isempty(id)
                DATA.currentpoint(1) = id(1);
            end
        end
    else
        if DATA.currentpoint(1) < nex
            DATA.currentpoint(1) = DATA.currentpoint(1)+1;
        else
            DATA.currentpoint(1) = nex;
            DATA.currentpoint(2) = DATA.currentpoint(2)+1;
        end
    end
        
        
elseif dir == 'u' %step 1 up
    if isfield(DATA,'usepeaks') && DATA.steptype ==2
        if DATA.currentpoint(1) > 1
            j = DATA.currentpoint(1)-1;
            id = find(DATA.usepeaks(1:j) > 0);
            if ~isempty(id)
                DATA.currentpoint(1) = id(end);
                DATA.currentpoint(2) = DATA.usepeaks(id(end));
            end
        else
            id = find(DATA.usepeaks> 0);
            if ~isempty(id)
                DATA.currentpoint(1) = id(1);
            else
                DATA.currentpoint(1) = 1;
            end
        end
    else
        if DATA.currentpoint(1) > 1
    DATA.currentpoint(1) = DATA.currentpoint(1)-1;
    else
        DATA.currentpoint(1) = 1;
        if DATA.currentpoint(2) > 1
            DATA.currentpoint(2) = DATA.currentpoint(2)-1;
        end
        end
    end
    DATA = ClearSelections(DATA,0,0);
end
DATA = ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2),'oneprobe');
if DATA.plotallxy
    PlotAllProbeXY(DATA);
    drawnow;
end
if ismember(dir,'ud')  %moving expts up and down with arrow keys, spoole/show all probes
    e = DATA.currentpoint(1);
    PrintComments(DATA,e,1:DATA.nprobes);
    if DATA.spoolspikes
        PlotAllProbe(DATA,'spoolall');
    elseif DATA.show.watchexptallspks
        PlotAllProbe(DATA,'allspks');
    end
    DATA = ClearSelections(DATA, 0, 0);
    str = sprintf('Ex %.1f: %s',DATA.exptid(DATA.currentpoint(1)),DATA.expnames{DATA.currentpoint(1)});
    fprintf('%s\n',str);
    SetFigure(DATA, DATA.tag.all);
    title(str);
    Expts = getappdata(DATA.toplevel,'Expts');
    DATA.Expt = Expts{e,1};
    SetTrialList(DATA,{},DATA.currenttrial);
end
stopped = 0;
set(DATA.toplevel,'UserData',DATA);


function DATA = SpoolCurrentSpikes(mn,b)
    DATA = GetDataFromFig(mn);
    AllSpikes = CheckAllSpikes(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
    DATA = get(DATA.toplevel,'UserData');
    
    if isempty(AllSpikes)
        a = SpoolSpikeFile(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
        if isfield(a,'svfig')
            DATA.spoolfig = a.svfig;
            set(DATA.toplevel,'UserData',DATA);
        end
    else
        Clusters = getappdata(DATA.toplevel,'Clusters');
        SetFigure(DATA,DATA.tag.xyplot,'front');
        hold off;
        plots = PlotClusterXY(DATA,Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)});
        drawnow;
        SetTrialList(DATA,Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)},DATA.currenttrial);
        SpoolSpikes(DATA, DATA.currentpoint);
    end
    
    
    
function res = SpoolSpikeFile(DATA,e,p)
        
    args = {};
    Clusters = getappdata(DATA.toplevel,'Clusters');
    Expts = getappdata(DATA.toplevel,'Expts');
    if DATA.NewCut.exptid == e && isfield(DATA.NewCut,'clst')
        C = DATA.NewCut;
    else
    C = Clusters{e}{p};
    end
    xs = '';
    if rem(C.exptno,1) > 0.001
        xs = 'a';
    end
    [a,b] = fileparts(DATA.name);
    [c,d] = fileparts(a);
    if DATA.nprobespool > 0
        chspk = p-DATA.nprobespool:p+DATA.nprobespool;
        chspk = chspk(chspk > 0 & chspk <=24 & chspk ~= p);
%        chspk = [p+1];
        for j = 1:length(chspk)
            args = {args{:} 'spkfile' [DATA.name '/Spikes/' b '.p' num2str(chspk(j))  't' num2str(floor(C.exptno)) xs '.mat']};
        end
    end
    name = [DATA.name '/Spikes/' b '.p' num2str(p)  't' num2str(floor(C.exptno)) xs '.mat'];
    xcl = FindExcludedTrials(DATA,e,p,1, C);
    fprintf('Spooling Expt %d(row %d) Probe %d\n',C.exptno,e,p);
    if max(xcl) > length(Expts{e,p}.Trials)
        res = PlotSpikeFile(name,'Trials',Expts{e,p},args{:},C);
    else
        res = PlotSpikeFile(name,'Trials',Expts{e,p},'exclude',xcl,args{:},C);
    end

function [x, details] = FindExcludedTrials(DATA,e,p, cluster, C)
    details.xid = [];
    x = [];
    if ~isfield(DATA,'trialids') || cluster < 1
        return;
    end
    if isfield(C,'missingtrials')
        x = find(ismember(DATA.trialids{e},C.missingtrials));
    end
    if isfield(C,'excludetrialids')
        y = find(ismember(DATA.trialids{e},C.excludetrialids));
        x = [x y];
        details.xid = C.excludetrialids;
    end
    details.nt = length(DATA.trialids{e});
    details.fraction = 0;
    
    if ~isfield(DATA.CellDetails,'excludetrials')
        DATA.CellDetails.excludetrials=[];
    end
    XC = DATA.CellDetails.excludetrials;
  
        
    if size(XC,1) < e || size(XC,2) < p || size(XC,3) < cluster
        return;
    end
    y = find(ismember(DATA.trialids{e},XC{e,p, cluster}));
    x = [x y];
    details.xid = DATA.trialids{e}(x);
    details.fraction = length(x)./details.nt;
    return;
    
    %Old stuff
    if p <= 0 || size(DATA.CellList,2) < p || size(DATA.CellList,1) < e
        return;
    end
    cellid = DATA.CellList(e,p, cluster);
    if cellid > 0 && length(XC) >= cellid && length(XC{cellid}) >= e
        x = [x XC{e, cellid, cluster}];
        details.fraction = length(x)./details.nt;
    end


     
function go = CheckSpikeFile(DATA, C)

name = SpikeFileName(DATA,C);
d = dir(name);
tdiff = C.savetime(1)-datenum(d.date); 
if tdiff > 0.001 %at least 1 min old to avoid rounding in d.date (string)
    fprintf('%s %.2f older\n',name,tdiff);
    go = 1;
else
   if tdiff > 0.0001
    go = 0;
   end
    go = 0;
end

function Get2DMaxima(DATA)
        
mid = DATA.mahaltype;
X = squeeze(DATA.mahal(:,:,mid));
[a,b,G] = Gauss2D(1,[-3:1:3]);
Y = conv2(X,G,'same');
SetFigure(DATA,DATA.tag.smoothed);
imagesc(Y);
xm = diff(sign(diff(Y)));
[a,b] = find(xm < 0);
ym = diff(sign(diff(Y,1,2)),1,2);
[a,b] = find(xm(:,2:end-1) < 0 & ym(2:end-1,:) < 0);
hold on;
mxs = xm(:,2:end-1) + ym(2:end-1,:)
plot(b+1,a+1,'r+');
title(sprintf('%d maxima',length(b)));



function name = SpikeFileName(DATA, C)
            xs = '';
            p  = C.probe;
    if rem(C.exptno,1) > 0.001
        xs = 'a';
    end
    [a,b] = fileparts(DATA.name);
    [c,d] = fileparts(a);
    name = [DATA.name '/Spikes/' b '.p' num2str(p)  't' num2str(floor(C.exptno)) xs '.mat'];


function [handles, details] = PlotPopPoints(X,Y, varargin)

j = 1;
colorscheme = 'none';
shapelist = [];
include = ones(size(X));
while j <= length(varargin)
    if isstruct(varargin{j})
        DATA = varargin{j};
    elseif strncmpi(varargin{j},'exclude',6)
        j = j+1;
        if strcmp(varargin{j},'offpeak') %trigger chan not biggest in mean
            id = find(abs(DATA.peakdiff) > 0);
            include(id) = 0;
        end
    elseif strncmpi(varargin{j},'colors',6)
        j = j+1;
        colorscheme = varargin{j};
    elseif strncmpi(varargin{j},'includelist',10)
        j = j+1;
        include = varargin{j};
    elseif strncmpi(varargin{j},'shapes',6)
        j = j+1;
        shapelist = varargin{j};
    end
    j = j+1;
end
    
colors = mycolors;
types = 1;
if isnumeric(colorscheme)
    types = colorscheme;
    c = abs(colorscheme);
    fill = 1-sign(colorscheme);
elseif strmatch(colorscheme,'probe')
    colors = mycolors('ncol',24);
    c = repmat([1:size(X,2)],size(X,1),1);
    fill = zeros(size(c));
elseif strmatch(colorscheme,'expt')
    colors = mycolors('ncol',size(X,1));
    c = repmat([1:size(X,1)]',1,size(X,2));
    fill = zeros(size(c));
elseif strmatch(colorscheme,'cutspace')
    c = DATA.cutspace;
    fill = zeros(size(c));
else
    c = ones(size(X));
    fill = zeros(size(c));
end
typelist = unique(types);
hold off;
np = 0;
for j = 1:size(X,1)
    for k = 1:size(X,2)
        if include(j,k)
            if isempty(shapelist)
                shape = 'o';
            else
                shape = shapelist(j,k);
            end
        h = plot(X(j,k),Y(j,k),shape,'buttondownfcn',{@HitPopPoint, j,k},'color',colors{c(j,k)});
        np = np+1;
        if fill(j,k)
            set(h,'markerfacecolor',colors{c(j,k)});
        end
        if length(types) > 1
            type = find(types(j,k) == typelist);
            handles(type) = h;
        else
            handles(1) = h;
        end
        hold on;
        end
    end
end
details.np = np;

function PlotMahalImage(DATA, type, varargin)
colorbarpos = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'colorbarpos',10)
        j = j+1;
        colorbarpos = varargin{j};
    end
    j = j+1;
end

mid = DATA.mahaltype;
if isempty(type) || strcmp(type,'mahal') && mid < 5
    imdata = squeeze(DATA.mahal(:,:,mid));
elseif strncmp(type,'fitdprime',5) || mid == 5
    imdata = squeeze(DATA.GaussFitdp(:,:,1));
elseif strncmp(type,'autodiff',6)
    mid = DATA.mahaltype;
    imdata = squeeze(DATA.mahal(:,:,mid)-DATA.automahal(:,:,mid));
elseif strncmp(type,'BuildDates',6)
    imdata = squeeze(DATA.ctimes(:,:,1));
end

cmenu = AddContextMenu(DATA,'cellplot');
imagesc(imdata,'buttondownfcn',{@HitImage, 1},'uicontextmenu',cmenu);
set(gcf, 'KeyPressFcn',{@KeyPressed, 1},'Keyreleasefcn',{@KeyReleased, 3},'uicontextmenu',cmenu);
 
if isempty(colorbarpos)
h = colorbar('EastOutside');
else
    h = colorbar('Position',colorbarpos);
end
%set(h,'position',[0.9 0.05 0.05 0.9])
hold on;
if DATA.plot.mahalcmax > 0
caxis([0 DATA.plot.mahalcmax]);
end
for j = 1:size(DATA.peakpos,2)
    plot(cat(1,j, DATA.peakpos(:,j)),cat(2,0.5,1:size(DATA.peakpos,1)),'w-',...
    'buttondownfcn',{@HitImage, 1},'uicontextmenu',cmenu);
    hold on;
end
for j = 1:size(DATA.mahal,1)
    text(1,j,ExLabel(DATA,j),'color','k');
end
Clusters = getappdata(DATA.toplevel,'Clusters');
        for j = 1:length(Clusters)
            good = find(sum(CellToMat(Clusters{j},'space'),2) > 0)';
            for k = good
                if Clusters{j}{k}.auto == 0 && DATA.plot.markmanual
                    DrawBox(j,k,1,'color','r','uicontextmenu',cmenu);
                end
                if Clusters{j}{k}.auto == 1 && DATA.plot.markauto
                    DrawBox(j,k,1,'color','y','uicontextmenu',cmenu);
                end
                if isfield(Clusters{j}{k},'marked') && Clusters{j}{k}.marked == 4 && DATA.plot.markgoodmu
                    DrawBox(j,k,1,'color','y','uicontextmenu',cmenu);
                end
                if Clusters{j}{k}.quick == 1 && DATA.markcell.quick
                    DrawBox(j,k,1,'color','y','uicontextmenu',cmenu);
                end
                if DATA.markcell.quick && Clusters{j}{k}.manual == 2
                    DrawBox(j,k,1,'color','m','uicontextmenu',cmenu);
                end
                if DATA.plot.markcopy && isfield(Clusters{j}{k},'manual') && Clusters{j}{k}.manual == 2
                    DrawBox(j,k,1,'color','g','uicontextmenu',cmenu);
                end
                if isfield(DATA.CellDetails,'excludetrials') && length(DATA.CellDetails.excludetrials{j,k}) > 0
                    DrawBox(j,k,1,'color','m','uicontextmenu',cmenu);
                end
                if DATA.plot.markgood && isfield(Clusters{j}{k},'marked') && Clusters{j}{k}.marked == 2
                    DrawBox(j,k,1,'color','g','uicontextmenu',cmenu);
                end
            end
        end
 set(gca,'UserData',DATA.toplevel);
 set(gcf,'UserData',DATA.toplevel);
    
function RateMenu(src, ks, fcn)
DATA = GetDataFromFig(src);
Expts = getappdata(DATA.toplevel,'Expts');
D = get(gcf,'UserData');
f = gcf;

    if strcmp(fcn,'hzoomout')
        RateZoom(DATA, 3,DATA.currentcell);
    elseif strcmp(fcn,'vzoomout')
        RateZoom(DATA, 2,DATA.currentcell);
    elseif strcmp(fcn,'zoomin')
        RateZoom(DATA, 1,DATA.currentcell);
    elseif strcmp(fcn,'spool')
        SpoolSpikes(DATA, DATA.currentpoint);
    elseif strcmp(fcn,'xyseq')
        PlotXYSequence(DATA, DATA.currentpoint);
    end
    figure(f);
    
function RateSeqKeyPressed(src, ks, fcn)

    if strmatch(ks.Key,{'shift' 'control' 'alt'})
        return;
    end
if ~(ismember(ks.Character,'+-' ) | sum(strcmp(ks.Key,{'leftarrow' 'rightarrow' 'uparrow' 'downarrow'})))
    return;
end

DATA = GetDataFromFig(src);
Expts = getappdata(DATA.toplevel,'Expts');
D = get(gcf,'UserData');
if ~isfield(D,'currentexpt') || D.currentexpt > length(D.exptlist)
    D.currentexpt = 1;
end
tag = get(src,'tag');
bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
if strcmp(ks.Key,'rightarrow')
    if D.currentexpt < length(D.exptlist)
        D.currentexpt = D.currentexpt+1;
    end
elseif strcmp(ks.Key,'leftarrow')
    if D.currentexpt > 1
        D.currentexpt = D.currentexpt-1;
    end
elseif strcmp(ks.Key,'uparrow')
    if DATA.currentcell < max(DATA.CellList(:))
        DATA.currentcell = DATA.currentcell+1;
    end
elseif strcmp(ks.Key,'downarrow')
    if DATA.currentcell > 1
        DATA.currentcell = DATA.currentcell-1;
    end
end
eid = D.exptlist(D.currentexpt);
DATA.currentpoint(1) = eid;
set(gcf,'UserData',D); 
if isempty(ks.Character)
elseif ismember(ks.Character,'+') || sum(strcmp(ks.Key,{'leftarrow' 'rightarrow'})) 
    t(1) = Expts{eid}.Header.timeoffset+Expts{eid}.Trials(1).TrialStart/10000;
    t(2) = Expts{eid}.Header.timeoffset + Expts{eid}.Trials(end).TrialStart*1.01/10000;
    t = t - D.timeadjust(D.currentexpt);
    set(gca,'xlim',t);
    if strcmp(ks.Modifier,'shift')
        DATA = RateZoom(DATA, 1, DATA.currentcell);
    else
        DATA= RateZoom(DATA,2, DATA.currentcell);
        title(sprintf('Expt%d %s',eid,DATA.expnames{eid}));
    end
end
if  sum(strcmp(ks.Key,{'uparrow' 'downarrow'})) 
    D = get(gcf,'UserData');
    if length(D.cellids) > 1
    DATA = RateZoom(DATA, 1, DATA.currentcell);  
    else
        PlotMenu(DATA, [], 'cells', 'rateseqone');
    end
end

if ks.Character == '-'
    X = get(gcf,'UserData');
    set(gca,'ylim',X.yrange);
    if strcmp(ks.Modifier,'shift')
        set(gca,'xlim',X.xrange);
    end
end
set(DATA.toplevel,'UserData',DATA);

function DATA = RateZoom(DATA,inout, cell)
%DATA = RateZoom(DATA,inout, cell)
% inout = 1 one cell/one expt
% inout = 2 zoom out Y axis only
% inout = 3 zoom out all the way X axis onl%y
% inout = 5 one cell, all expts
    D = get(gcf,'UserData');
    eid = DATA.currentpoint(1);
    if ~isfield(D,'handles')
        return;
    end
    if cell > length(D.handles)
        th = title(sprintf('Expt%d %s cell %d No Spikes',eid,DATA.expnames{eid},cell));
        return;
    else
        h = D.handles(cell);
    end
    if isfield(D,'cllabel') && D.cllabel > 0 && ishandle(D.cllabel)
        delete(D.cllabel);
    end
    [a,b,cl] = FindCell(DATA, cell, eid);
    if isempty(cl)
        cl = 0;
    else
    DATA.currentpoint(2) = b;
    DATA.currentcellcluster = cl;
    end
    th = title(sprintf('Expt%d %s cell %d ',eid,DATA.expnames{eid},cell),'color','k');
    if inout == 1 || inout == 5%zoom in Y axis
        if inout == 5
            set(gca,'xlim',D.xrange);
        end
        [yl, a] = GetYRange(h, get(gca,'xlim'));
        if a.result > 0
           set(th,'color',a.color);
           set(gca,'ylim',yl);
        else
            th = title(sprintf('Expt%d %s cell %d No Spikes',eid,DATA.expnames{eid},cell));
            if diff(D.RateRange(cell,:)) > 0
                set(gca,'ylim',D.RateRange(cell,:));
            end
        end
    elseif ismember(inout,[2])  %zoom out Y only
        yl = minmax(D.RateRange(:));
        yls = [];
        for j = 1:length(D.handles)
            [yl, a] = GetYRange(D.handles(j), get(gca,'xlim'));
            if a.result > 0
                yls = [yls yl];
            end
        end
        if length(yls)
        set(gca,'ylim',minmax(yls));
        end
    elseif inout == 3
        set(gca,'xlim',D.xrange);
    elseif inout == 4
        set(gca,'ylim',D.yrange);
        set(gca,'xlim',D.xrange);
    end
    %do this after any zooming
    D.cllabel = AddClusterString(DATA, th,cl);
    set(gcf,'UserData',D);

function [yl, details] = GetYRange(h, xl)
    details.result = 0;
    if h <= 0 || ~ishandle(h)
        yl = get(gca,'ylim');
        details.color = 'b';
        return;
    end
    X = get(h,'xdata');
    Y = get(h,'ydata');
    id = find(X >= xl(1) & X <= xl(2));
    details.color = get(h,'color');
    if length(id) <= 1
        yl = get(gca,'ylim');
    else
        yl = minmax(Y(id));
        if diff(yl) == 0
            yl(1) = yl(1)+1;
        end
        details.result = 1;
    end
        
    
function ch = AddClusterString(DATA, h, cl)
    if isempty(cl)
        cl = 0;
    end
    if h > 0 && ishandle(h)
        x = get(h,'extent');
        ch = text(x(1)+x(3),x(2)+x(4)/2,sprintf('cl%d',cl),'color',DATA.colors{cl+1});
    else
        ch = [];
    end
    
    
  function XYKeyPressed(src, ks, fcn)

      if strcmp(ks.Key, 'control')
          setappdata(0,'control_is_down',1);
          return;
      elseif strcmp(ks.Key,{'alt'})
          setappdata(0,'alt_is_down',1);
          return;
      elseif sum(strcmp(ks.Key,{'shift' 'control' 'alt'}))
          return;
      end
      
DATA = GetDataFromFig(src);
tag = get(src,'tag');
bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
if sum(strcmp(ks.Key,{'uparrow' 'downarrow'})) & strcmp(ks.Modifier,'shift')
    if DATA.NewCut.saved == 0
        DATA = SaveCluster(DATA, DATA.currentpoint,1);
    end
    DATA = AddSelectedCells(DATA);
end


if strcmp(tag,DATA.tag.allxy) && DATA.allclustering
    if strcmp(ks.Key,'l')
        XYCluster(src, [],'lines');
    elseif strcmp(ks.Key,'1')
        XYCluster(src, [],'ellipse1');
    elseif strcmp(ks.Key,'2')
        XYCluster(src, [],'ellipse2');
    elseif strcmp(ks.Key,'3')
        XYCluster(src, [],'ellipse3');
    elseif strcmp(ks.Key,'4')
        XYCluster(src, [],'ellipse4');
    elseif strcmp(ks.Key,'f')
        DATA = FlipLineCrit(DATA);
    elseif strcmp(ks.Key,'q')
        XYCluster(src, [],'quicksave');
    elseif strcmp(ks.Key,'n')
        XYCluster(src, [],'nocut');
    end
else
    if strcmp(ks.Key,'f')
        DATA = FlipLineCrit(DATA);
    elseif strcmp(ks.Key,'d')
        DATA.plot.density = ~DATA.plot.density;
        ReplotXY(DATA,[],DATA.currentpoint(1),DATA.currentpoint(2),DATA.currentcluster);
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmp(ks.Key,'q')       
        XYCluster(src, [],'quicksave');
    elseif strcmp(ks.Key,'add') | ks.Character == '+' %characte can be empty
        DATA = AddSelectedCells(DATA);    
    end
end
if strcmp(ks.Key,'rightarrow')
    DATA = NextButton(src, ks, 'r');
elseif strcmp(ks.Key,'leftarrow')
    DATA = NextButton(src, ks, 'l');
elseif strcmp(ks.Key,'downarrow')
    DATA =  NextButton(src, ks, 'd');
elseif strcmp(ks.Key,'uparrow')
    DATA = NextButton(src, ks, 'u');
elseif strcmp(ks.Key,'s')
    SpoolSpikes(DATA, DATA.currentpoint);
end
set(DATA.toplevel,'UserData',DATA);

function DATA = FlipLineCrit(a,b)
    DATA = GetDataFromFig(a);
    
    if DATA.NewCut.probe > 0
        DATA.elmousept.pos = DATA.NewCut.pos
    else
    end
    C =  ClassifySpikes(DATA,DATA.elmousept,'flip');
    DATA = get(DATA.toplevel,'UserData');
    DATA = ConditionalPlotXY(DATA, C, 1);


function DATA = AddSelectedCells(DATA)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    cellid = DATA.currentcell;
    DATA.selectprobe(DATA.currentpoint(1),DATA.currentpoint(2)) = 1;
    [a,b] = find(DATA.selectprobe);
    for j  = 1:length(a)
        DATA = SetCellEntry(DATA, Clusters{a(j)}{b(j)}, a(j),b(j),DATA.cellcluster, cellid,'nosave');
    end
    SaveCellList(DATA);
    DATA = MarkCurrentCluster(DATA);
    PlotCellList(DATA,'showfig');
    set(DATA.toplevel,'UserData',DATA);

function KeyReleased(src, ks, fcn)

      if strcmp(ks.Key, 'control')
          setappdata(0,'control_is_down',0);
          return;
      elseif strmatch(ks.Key,{'alt'})
          setappdata(0,'alt_is_down',0);
          return;
      end
      
function KeyPressed(src, ks, fcn)

      if strcmp(ks.Key, 'control')
          setappdata(0,'control_is_down',1);
          return;
      elseif strcmp(ks.Key,{'alt'})
          setappdata(0,'alt_is_down',1);
          return;
      elseif sum(strcmp(ks.Key,{'shift' 'control' 'alt'}))
          return;
      end
DATA = GetDataFromFig(src);
tag = get(src,'tag');
figtag = get(gcf,'tag');

bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});

%Shift up/down saves any cut cluster, set the currecnt cell, then moves on
if sum(strcmp(ks.Key,{'uparrow' 'downarrow'})) & strcmp(ks.Modifier,'shift')
    if DATA.NewCut.saved == 0
        DATA = SaveCluster(DATA, DATA.currentpoint,1);
    end
    DATA = AddKey(DATA,figtag);
end

if strmatch(ks.Key,{'delete' 'backspace'}) 
    if fcn == 2 % delete box marking cell on template track plot
        DATA.usepeaks(DATA.currentpoint(1)) = 0;
        if ishandle(DATA.markcc)
            delete(DATA.markcc)
        end
        h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),2);
        set(h,'color','r');
        if strcmp(ks.Modifier,'shift')
            eid = find(DATA.usepeaks > 0);
            for j = 1:length(eid)
                h = DrawBox(eid(j),DATA.usepeaks(eid(j)),'color','r');
            end
            DATA.usepeaks = zeros(size(DATA.usepeaks));
        end
    elseif fcn == 3  %delete current location from CellList
        it = findobj(src,'Tag','CellNumberId');
        if length(it) == 1
            cellid = get(it,'value');
            DATA.selectprobe(DATA.currentpoint(1),DATA.currentpoint(2)) = 1;
            [a,b] = find(DATA.selectprobe);
            for j  = 1:length(a)
                id = find(ismember(DATA.CellList(a(j),b(j),:),cellid));
                if length(id) == 1
                    DATA.cellcluster = id;
                    DATA = DeleteCell(DATA, a(j), b(j), DATA.cellcluster);
                end
            end
            SaveCellList(DATA);
        end
        DATA = PlotCellList(DATA,'showfig');
    end
elseif strcmp(ks.Key,'add') | ks.Character == '+'
    DATA = AddKey(DATA,figtag);
    
elseif strmatch(ks.Key,'rightarrow')
    DATA = NextButton(src, ks, 'r');
elseif strmatch(ks.Key,'leftarrow')
    DATA = NextButton(src, ks, 'l');
elseif strcmp(ks.Key,'downarrow')
    DATA =  NextButton(src, ks, 'd');
elseif strmatch(ks.Key,'uparrow')
    DATA = NextButton(src, ks, 'u');
    n =1;
elseif strmatch(ks.Key,'f')
    DATA = FlipLineCrit(DATA);
elseif strmatch(ks.Key,'q')
    DATA = SaveCluster(DATA, DATA.currentpoint, 1);
elseif strmatch(ks.Key,'add')
elseif strmatch(ks.Key,'subtract')
elseif strmatch(ks.Key,'space')
end
set(DATA.toplevel,'UserData',DATA);

    if strcmp(tag,DATA.tag.celllist) && bt == 1
        figure(src);
    end

    
function DATA = AddKey(DATA, tag, varargin)
    if strcmp(tag,DATA.tag.celllist)
        DATA = AddSelectedCells(DATA);
        imtype = 3;
    else
        if DATA.usepeaks(DATA.currentpoint(1)) > 0
            h = DrawBox(DATA.currentpoint(1),DATA.usepeaks(DATA.currentpoint(1)),1);
            set(h,'color','r');
        end
        DATA.usepeaks(DATA.currentpoint(1)) = DATA.currentpoint(2);
        imtype = 1;
    end
    h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),imtype );
    set(h,'color','w');
    set(DATA.toplevel,'UserData',DATA);

function DATA = DeleteCell(DATA, e, p, cl)
    oldcell = DATA.CellList(e,p,cl);
    DATA.CellList(e,p,cl) = 0;
    DATA.CellChanges = cat(1,DATA.CellChanges,[e p 0 cl now oldcell]);
    SaveCellList(DATA);
          
    
function s = OldExLabel(DATA, j)
 s = '';
 Clusters = getappdata(DATA.toplevel,'Clusters');
 Expts = getappdata(DATA.toplevel,'Expts');
 k = floor(Clusters{j}{1}.exptno);
 if DATA.show.exptno
     s = [s num2str(k)];
 end
 if DATA.show.exptname && ~isempty(Expts)
     s = [s ':' Expts{j,1}.Header.expname];
 end
 if DATA.show.ed && ~isempty(Expts)
     s = [s ''  sprintf('%.2f', mean([Expts{j,1}.Trials.ed]))];
 end

 function s = ExLabel(DATA, j)
 s = '';

 if DATA.show.exptno
     s = [s num2str(DATA.exptid(j))];
 end
 if DATA.show.exptname
     s = [s ':' DATA.expnames{j}];
 end
 if DATA.show.ed 
     s = [s ''  sprintf('%.2f', DATA.electrodedepth(j))];
 end

function DATA = PlotAllClusters(mn,b, varargin)
    DATA = GetDataFromFig(mn);
    Expts  = getappdata(DATA.toplevel,'Expts');
    
    if ~ishandle(mn)
        plottype = DATA.plot.alltype;
        mn = findobj(DATA.toplevel,'Tag','plotalltype');
    else
        str = get(mn,'string');
        a = get(mn,'value');
        plottype = deblank(str(a,:));
        DATA.plot.alltype = plottype;
    end
    Clusters = getappdata(DATA.toplevel,'Clusters');
    SetFigure(DATA,DATA.tag.clusters);
    ts = [];
    durs = [];
    ns = [];
    DATA.markh = NaN;
    plotargs = PlotArgs(DATA);
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'includelist',10)
            plotargs = {plotargs{:} varargin{j} varargin{j+1}};
            j = j+1;
        end
        j = j+1;
    end
    
    if strcmp(plottype,'spkshapecorr')
        CompareShapes(DATA,1);
        return;
    elseif strcmp(plottype,'followcorr')
        DATA = TrackCluster(DATA,DATA.templateid(1),DATA.templateid(2));
        return;
    elseif strcmp(plottype,'spkshape')
        CompareShapes(DATA,'exptlines');
        return;
    elseif strcmp(plottype,'spkshapeim')
        DATA = CompareShapes(DATA,'exptimage');
        set(DATA.toplevel,'UserData',DATA);
        return;
    elseif strcmp(plottype,'spkshape')
        CompareShapes(DATA,2);
        return;
    elseif strcmp(plottype,'spkpeaks')
        subplot(1,1,1);
        hold off;
        if 1
        imagesc(DATA.peakpos,'buttondownfcn',{@HitImage, 1})
        set(gcf, 'KeyPressFcn',{@KeyPressed, 1});
        set(gcf, 'KeyReleaseFcn',{@KeyReleased, 1});
        set(gca,'UserData',DATA.toplevel);
        else
        for j = 1:size(DATA.peakpos,2)
            plot(1:size(DATA.peakpos,1),DATA.peakpos(:,j),'-');
            hold on;
        end
        end
        return;
    elseif strcmp(plottype,'spksize')
        SpikeDistances(DATA,1);
        return;
    elseif strcmp(plottype,'probexcorr')
        CompareProbes(DATA, 1);
        return;
    elseif strcmp(plottype,'probeco')
        CompareProbes(DATA, 2);
        return;
        return;
    elseif strcmp(plottype,'Exclusions')
        for j = length(Clusters):-1:1;
            for k = length(Clusters{j}):-1:1
                C = Clusters{j}{k};
                [xcl, details] = FindExcludedTrials(DATA,j,k,1,C);
                if length(xcl)
                    X(j,k) = -details.fraction;
                elseif isfield(C,'restricttimerange')
                    exptdur = Expts{j,k}.Trials(end).End(end) - Expts{j,k}.Trials(1).Start(1);
                    X(j,k) = diff(C.restricttimerange)./(exptdur./10000);
                else
                    X(j,k) = 0;
                end
            end
        end
        imagesc(X,'ButtondownFcn',{@HitImage, 1});
        return;
    elseif strcmp(plottype,'SpkRate')
        SetFigure(DATA,DATA.tag.clusters);
        for j = 1:size(DATA.nspks,1)
            m(j,:) = 10000 .* DATA.nspks(j,:)./sum([Expts{j,1}.Trials.dur]);
        end
        imagesc(m,'buttondownfcn',{@HitImage, 1});
        colorbar;
        return;
    elseif strcmp(plottype,'EventRate')
        SetFigure(DATA,DATA.tag.clusters);
        for j = 1:size(DATA.nevents,1)
            m(j,:) = 10000 .* DATA.nevents(j,:)./sum([Expts{j,1}.Trials.dur]);
        end
        imagesc(m,'buttondownfcn',{@HitImage, 1});
        colorbar;
        return;
    elseif strcmp(plottype,'TriggerFind')
        SetFigure(DATA,DATA.tag.clusters);
        m = max(DATA.mahal(:,:,[1 3]),[],3);
        hold off; 
        imagesc(m,'buttondownfcn',{@HitImage, 4});
        caxis([0 5]);
        hold on;
        [a,b] = find(m > 2 & squeeze(DATA.dropi(:,:,1)) < 1 & abs(DATA.peakdiff) < 1);
        for j = 1:length(a);
            DrawBox(a(j),b(j),4,'color',[0.5 0.5 0.5]);
        end
        vc = prctile(DATA.mspkvar(:),50);
        [a,b] = find(m < 2 & abs(DATA.peakdiff) < 1 & DATA.muampl > 0.9 & DATA.mspkvar > vc);
        for j = 1:length(a);
            DrawBox(a(j),b(j),4,'color',[0.0 0.0 0.0]);
        end
    elseif strcmp(plottype,'Tagged')
        imagesc(DATA.tagged,'ButtondownFcn',{@HitImage, 1});
    end
    
    if DATA.rebuild || (strcmp(plottype,'xcorr') && ~isfield(DATA,'synci'))
    for j = 1:length(DATA.strings)
        if length(Clusters{j}) > 1
            SetFigure(DATA,DATA.tag.clusters);
            clf
            if strcmp(plottype,'dips')
            fprintf('Plotting %s at %s\n',DATA.strings{j},datestr(now));
            a = PlotSpikeTimes(Clusters{j},plottype);
            n = size(a.dips,2);
            mahal(j,1:n) = a.dips(6,:);
            dips(j,1:n) = a.dips(1,:);
            ndmahal(j,1:n) = a.dips(7,:);
            DATA.qres{j} = a;
            elseif strncmp(plottype,'mahal',5) & DATA.rebuild
                fprintf('Plotting %s at %s\n',DATA.strings{j},datestr(now));
                a = PlotSpikeTimes(Clusters{j},'bestspace');
                n = size(a.dips,2);
                DATA.mahal(j,1:n,1) = a.dips(6,:);
                dips(j,1:n) = a.dips(1,:);
                DATA.mahal(j,1:n,2) = a.dips(7,:);
                DATA.qres{j} = a;  
            elseif strncmp(plottype,'xcorr',5)
                ofile = [DATA.name '/' strrep(DATA.strings{j},'ClusterTimes','Xcorrs')];
                if ~exist(ofile,'file') || DATA.rebuild == 1
                    xc = PlotSpikeTimes(Clusters{j},plottype);
                    save(ofile,'xc');
                    fprintf('Plotting %s (%d cells) at %s\n',DATA.strings{j},length(xc.cells),datestr(now));
                    ts = [ts xc.times];
                    durs = [durs xc.xctime];
                    ns = [ns prod(xc.ns,2)'];
                else
                    load(ofile);
                    ns = [ns prod(xc.ns,2)'];
                end
                DATA.synci{j} = xc;
            elseif strcmp(plottype,'probesd')
                xc = PlotSpikeTimes(Clusters{j},plottype);
                n = length(diag(xc.sds));
                sds(j,1:n) = diag(xc.sds);
                allsd(j,1:size(xc.sds,1),1:size(xc.sds,2)) = xc.sds;
                DATA.cres{j} = xc;
            elseif strcmp(plottype,'refitgm')
                for k = 1:length(Clusters{j})
                    fprintf('E%dP%d\n',floor(Clusters{j}{k}.exptno),k);
                    [a,b] = GMDip(Clusters{j}{k}.xy,0);
                    hold off; 
                    [nsp,x] = hist(Clusters{j}{k}.xy(:,1),500);
                    bar(x,nsp,1);
                    scale = trapz(x,nsp)
                    hold on;
                    plot(b.gxy(:,1),b.gxy(:,2).*scale,'g');
                    plot(b.gxy(:,1),b.gxy(:,3).*scale,'g');
                    plot(b.gxy(:,1),sum(b.gxy(:,2:3),2).*scale,'r');
                    title(sprintf('E%.0fP%d',Clusters{j}{k}.exptno,k));
                    drawnow;
                end
            elseif strcmp(plottype,'Evec') && DATA.rebuild
                for k = 1:length(Clusters{j})
                    C = Clusters{j}{k};
                    if isfield(Clusters{j}{k},'Evec') && isfield(Clusters{j}{k}.Evec,'Eval') 
                        DATA.eveci(j,k) = Clusters{j}{k}.Evec.Eval(1)./sum(Clusters{j}{k}.Evec.Eval);
                        mahals(j,k) = Clusters{j}{k}.mahal(1);
                        xc = corrcoef(C.MeanSpike.ms(:),C.MeanSpike.mu(:));
                        DATA.mucorrs(j,k) = xc(1,2);
                    end
                end
            elseif DATA.rebuild
                DATA.cres{j} = PlotSpikeTimes(Clusters{j},plottype);
            end
        end
    end
    end
    
    if strmatch(plottype,{'drop-mahal' 'mahalxc' 'mahaldp' 'dxpc' 'PcGms' 'mahalvar' 'BuildTimes' 'Fit-mahal' 'stability'},'exact')
        SetFigure(DATA, DATA.tag.popscatter,'front');
    else
        SetFigure(DATA, DATA.tag.all,'front');
    end
    hold off;
    if strcmp(plottype,'dips')
        imagesc(mahal,'buttondownfcn',{@HitImage, 1});
    elseif strcmp(plottype,'mydip')
        for j = 1:length(Clusters)
            for k = 1:length(Clusters{j})
                C = Clusters{j}{k};
                if isfield(C,'mydipsize')
                    x(j,k) = C.mydipsize(1);
                    y(j,k) = C.fitdprime(1);
                end
            end
        end
        colors = mycolors;
        [a,b] = PlotPopPoints( x,y, DATA, 'colors', 'expt', plotargs{:});
        fprintf('%d points\n',b.np);
    elseif strcmp(plottype,'Fit-mahal') || strcmp(plottype,'stability')
        colors = ones(size(DATA.GaussFitdp,1),size(DATA.GaussFitdp,2));
        shapelist(1:size(DATA.GaussFitdp,1),1:size(DATA.GaussFitdp,2)) = 'o';
        for j = 1:size(DATA.GaussFitdp,1)
            for k = 1:size(DATA.GaussFitdp,2)
                if Clusters{j}{k}.shape == 0
                    shapelist(j,k) = 's';
                end
            end
        end
        id = find(DATA.cutspace);
        id = find(sum(DATA.gmfitpos,3) > 1);
        colors(id) = 2;
        
        cid = find(sum(DATA.gmfitpos,3) > 1 & sum(DATA.CellList,3) == 0 & abs(DATA.peakdiff) == 0);
        colors(cid) = -2;
        did = find(sum(DATA.gmfitpos,3) < 2 & sum(DATA.CellList,3) > 0);
        colors(did) = -1;
        if strcmp(plottype,'Fit-mahal') 
            h = PlotPopPoints(DATA.GaussFitdp(:,:,1), DATA.GaussFitdp(:,:,2),DATA,'colorscheme',colors,'shapes',shapelist,plotargs{:});
            xlabel('dprime from GM');
            ylabel('dprime from Gauss fits');
        elseif  strcmp(plottype,'stability')
            h = PlotPopPoints(DATA.GaussFitdp(:,:,2), DATA.xysdindex(:,:),DATA,'colorscheme',colors,'shapes',shapelist,plotargs{:});
        end
        mylegend(h,{'not a cell','fit overlap but celll','fitdp < 3','fitdp > 3'},'location','southeast');
    elseif strcmp(plottype,'drop-mahal')
        colors = ones(size(DATA.mahal,1),size(DATA.mahal,2));
        cid = find( abs(DATA.peakdiff) > 0);
        colors(cid) = 2;
        did = find(sum(DATA.CellList,3) > 0);
        colors(did) = colors(did).*-1;
        PlotPopPoints( max(DATA.mahal(:,:,[1 3]),[],3),DATA.dropi(:,:,1), DATA,'colorscheme', colors, plotargs{:});
    elseif strcmp(plottype,'muamp-spkvar')
        PlotPopPoints(DATA.muampl, DATA.mspkvar(:,:,1),DATA,plotargs{:});
    elseif strcmp(plottype,'shapexc')
        PlotPopPoints(DATA.xcorrs.xc, DATA.xcorrs.effic,DATA,plotargs{:});
    elseif strcmp(plottype,'mahalxc')
        PlotPopPoints(DATA.mahal(:,:,1), DATA.mucorrs,DATA,plotargs{:});
    elseif strcmp(plottype,'mahaldp')
        PlotPopPoints(DATA.mahal(:,:,1), DATA.dpsum,DATA,plotargs{:});
    elseif strcmp(plottype,'npts')
        for j = 1:length(Clusters)
            for k = 1:length(Clusters{j})
                npts(j,k) = range(Clusters{j}{k}.spts);
            end
        end
        imagesc(npts);
    elseif strcmp(plottype,'templatescore')
        subplot('position',[0.05 0.05 0.4 0.9]);
        imagesc(DATA.tmplxcs,'buttondownfcn',{@HitImage, 1});
        title('Correlation');
        for j = 1:length(Clusters)
            text(1,j,sprintf('%d',Clusters{j}{1}.cluster{1}.exptno),'HorizontalAlignment','left','color','w');
        end
         set(gca,'UserData',DATA.toplevel);
        subplot('position',[0.45 0.05 0.4 0.9]);
        set(gca,'yticklabel',[]);
        imagesc(squeeze(DATA.mahal(:,:,1)),'buttondownfcn',{@HitImage, 1});
        set(gca,'yticklabel',[]);
        h = colorbar('Position',[0.85 0.05 0.1 0.9]);
        if DATA.plot.mahalcmax > 0
            caxis([0 DATA.plot.mahalcmax]);
        end
        SetFigure(DATA,DATA.tag.templatesrc);
        PlotMeanSpike(DATA.Templates{DATA.templatesrc},0,1);

    elseif strcmp(plottype,'finetemplatescores')
        subplot(1,1,1);
        tmpim =  cat(2,DATA.tmplscores{:});
        tmpim = smooth(tmpim,100);
        imagesc(tmpim);
        hold on;
        nc = DATA.templatesrc;
         set(gca,'UserData',DATA.toplevel);
         ns = 0;
         for j = 1:length(DATA.tmplscores)
             ns = ns+size(DATA.tmplscores{j},2);
             plot([ns ns],[0 24],'r');
             text(ns,0,sprintf('%d',Clusters{j}{1}.cluster{nc}.exptno),'horizontalalignment','right');
         end
    elseif strcmp(plottype,'dpxc')
        PlotPopPoints(DATA.dpsum, DATA.mucorrs,DATA,plotargs{:});
    elseif strcmp(plottype,'spkrate')
        
        PlotPopPoints(DATA.dpsum, DATA.mucorrs,DATA,plotargs{:});
    elseif strcmp(plottype,'PcGms')
        [a,b] = max(DATA.pcgms,[],3);
        PlotPopPoints(DATA.pcgms(:,:,1), max(DATA.pcgms(:,:,2:end),[],3),DATA, 'colors',b,plotargs{:});
        xlabel('Gm Distance Plain');
        ylabel('Distance (Transformed)');
    elseif strcmp(plottype,'BuildTimes')
        buildtimes = diff(DATA.ctimes(:,:,[1 3]),[],3).*24;
        subplot(1,1,1);
        colors = ones(size(buildtimes));
        id = find(DATA.thriter > 0)
        colors(id) = 2;
        id = find(DATA.thriter > 1)
        colors(id) = 3;
        
        PlotPopPoints(DATA.nevents, buildtimes,DATA, 'colorscheme',colors,plotargs{:});
        ylabel('Duration (hrs)');
        xlabel('N events');
        title(sprintf('Total %.1f hours, elapsed %.1f hours',sum(buildtimes(:)),max((DATA.ctimes(:,end,3))-min(DATA.ctimes(:,1,1))).*24));
        te = max(DATA.ctimes(:,end,3));
        ts = min(DATA.ctimes(:,1,1));
        [a,id] = sort(DATA.ctimes(:,1,1));
        for j = 1:length(id)-1
            gap(j) = DATA.ctimes(id(j+1),1,1)-DATA.ctimes(id(j),end,3);
        end
        fprintf('First started at %s, last finished at %s Total gap %.1f\n',datestr(ts),datestr(te),sum(gap).*24);
    elseif strcmp(plottype,'mahal1-n')
        subplot(1,1,1);
        colors = mycolors;
        hold off;
        clabels = {'ndPC' 'ndADC' 'ndVar'  'ndTemplate' 'none' 'C' 'ADC' 'varE' 'Template'};
            
        colors = mycolors('ncol',24);
        for j = 1:size(DATA.mahal,1)
            for k = 1:size(DATA.mahal,2)
                if DATA.peakdiff(j,k) == 0
                h = plot(DATA.mahal(j,k,1),DATA.mahal(j,k,3),'ro','buttondownfcn',{@HitPopPoint, j,k});
                set(h,'color',colors{DATA.cutspace(j,k)});
                hs(DATA.cutspace(j,k)) = h;
                hold on;
                end
            end
        end
        xlabel('mahal N-D');
        ylabel('mahal 1D');
        mylegend(hs,clabels);
    elseif strcmp(plottype,'mahaln')
        subplot(1,1,1);
        hold off;
        colors = mycolors('ncol',24);
        for j = 1:size(DATA.mahal,1)
            for k = 1:size(DATA.mahal,2)
                if DATA.peakdiff(j,k) == 0
                plot(DATA.mahal(j,k,1),DATA.mahal(j,k,2),'o','buttondownfcn',{@HitPopPoint, j,k});
                hold on;
                plot(DATA.mahal(j,k,1),DATA.mahal(j,k,3),'ro','buttondownfcn',{@HitPopPoint, j,k});
                end
            end
        end
        xlabel('mahal 2-D');
        ylabel('mahal 1D(red)/ND(blue)');
    elseif strcmp(plottype,'mahal1-2')
        subplot(1,1,1);
        hold off;
        colors = mycolors('ncol',24);
        for j = 1:size(DATA.mahal,1)
            for k = 1:size(DATA.mahal,2)
                if DATA.peakdiff(j,k) == 0
                plot(DATA.mahal(j,k,3),DATA.mahal(j,k,1),'o','buttondownfcn',{@HitPopPoint, j,k});
                hold on;
%                plot(DATA.mahal(j,k,3),DATA.mahal(j,k,1),'ro','buttondownfcn',{@HitPopPoint, j,k});
                end
            end
        end
        xlabel('mahal 1D');
        ylabel('mahal 2D/ND');
    elseif strcmp(plottype,'mahalvar')
        PlotPopPoints(DATA.mahal(:,:,1), DATA.mspkvar,plotargs{:});

    elseif strcmp(plottype,'mahal+var')
        subplot('position',[0.05 0.05 0.4 0.9]);
        imagesc(DATA.mspkvar,'buttondownfcn',{@HitImage, 1});
         set(gca,'UserData',DATA.toplevel);
        subplot('position',[0.45 0.05 0.4 0.9]);
        set(gca,'yticklabel',[]);
        PlotMahalImage(DATA,'', 'colorbarpos',[0.85 0.05 0.1 0.9]);
        set(gca,'yticklabel',[]);
    elseif strcmp(plottype,'man-auto')
        subplot(1,1,1);
        PlotMahalImage(DATA, 'autodiff');
    elseif sum(strcmp(plottype,{'mahal' 'fitdprime' 'fitdp' 'BuildDates'}))
        subplot(1,1,1);
        PlotMahalImage(DATA, plottype);

    elseif strcmp(plottype,'Evec')
        subplot(1,2,1);
        imagesc(DATA.eveci);
        hold off; 
        colors = mycolors('ncol',24);
        for j = 1:size(DATA.eveci,1)
            for k = 1:size(DATA.eveci,2)
                plot(mahals(j,k),DATA.eveci(j,k),'o','buttondownfcn',{@HitPopPoint, j,k},'color',colors{k});
                hold on;
            end
        end
        subplot(1,2,2);
        hold off; 
        for j = 1:size(DATA.eveci,1)
            for k = 1:size(DATA.eveci,2)
                plot(DATA.mucorrs(j,k),DATA.eveci(j,k),'o','buttondownfcn',{@HitPopPoint, j,k},'color',colors{k});
                hold on;
            end
        end
    elseif strcmp(plottype,'probesd')
        subplot(1,2,1);
        imagesc(sds);
        subplot(1,2,2);
        imagesc(squeeze(mean(allsd)));
    elseif strncmp(plottype,'xcorr',5)
        allsync = [];
        allexpt = [];
        allcells = [];
        nc = 0;
        
        for j = 1:length(DATA.synci)
            for k = DATA.synci{j}.cells
            id = find(DATA.synci{j}.cells < k);
            cells = [ones(size(id)).* k; DATA.synci{j}.cells(id)];
            allsync = cat(2,allsync,DATA.synci{j}.synci(k, DATA.synci{j}.cells(id),1));
            allexpt = cat(2,allexpt,j .* ones(size(id)));
            allcells = cat(2, allcells, cells);
            end
            n = 0;
            for k = 1:length(DATA.synci{j}.cells)
                for m = 1:k-1;
                    nc = nc+1;
                    n = n+1;
                    a = DATA.synci{j}.cells(k);
                    b = DATA.synci{j}.cells(m);
                    ccf = DATA.synci{j}.ccf{a,b};
                    midpt = ceil(length(ccf)./2);
                    asd = std(Clusters{j}{a}.MeanSpike.ms');
                    bsd = std(Clusters{j}{b}.MeanSpike.ms');
                    xc = corrcoef(Clusters{j}{a}.MeanSpike.ms(:),Clusters{j}{b}.MeanSpike.ms(:));
                    ampratio(nc) = (asd(a)./asd(b)).*(bsd(b)./bsd(a));
                    synci(nc) = mean(ccf([midpt+1 midpt-1]))./ccf(midpt);
                    xcs(nc) = xc(1,2);
                    sep(nc) = abs(b-a);
                    pairs(nc,:) = [a b];
                    nspk(nc) = prod(DATA.synci{j}.ns(n,:));
                    xctimes(nc) = DATA.synci{j}.xctime(n);
                end
            end
            ns = nspk;
        end
        hold off;
        if strcmp(plottype,'xcorra')
            for j = 1:length(ampratio)
                if sep(j) > 5
                    c = 'k';
                elseif sep(j) == 4
                    c = 'b';
                elseif sep(j) == 3
                    c = 'y';
                elseif sep(j) == 2
                    c = 'g';
                else
                    c = 'r';
                end
                plot(ampratio(j), synci(j),'o','color',c,'buttondownfcn',{@HitXcorr, j, allexpt(j), allcells(:,j)});
                hold on;
            end
            xlabel('Amp Ratio');
            ylabel('Synci Local');
        elseif strcmp(plottype,'xcorrt')
            subplot(1,1,1);
            hold off; 
            for j = 1:length(ns)
                plot(ns(j),xctimes(j),'o','color','b','buttondownfcn',{@HitXcorr, j, allexpt(j), allcells(:,j)});
                hold on;
            end
        elseif strcmp(plottype,'xcorrb') || strcmp(plottype,'xcorrd') 
            ic = unique(pairs(:));
            subplot(1,2,1);
            hold off;
            subplot(1,2,2);
            hold off;
            nc = 0;
            for j = 1:length(ic)
                for k = 1:j-1
                    d = abs(ic(j)-ic(k));
                    if d > 5
                        c = 'k';
                    elseif d == 4
                        c = 'b';
                    elseif d == 3
                        c = 'y';
                    elseif d == 2
                        c = 'g';
                    else
                        c = 'r';
                    end
                    id = find(pairs(:,1) == ic(j) & pairs(:,2) == ic(k));
                    if length(id)
                    subplot(1,2,1);
                    nc = nc+1;
                    mpairs(nc,:) = [ic(j) ic(k)];
                    x(nc) = mean(ampratio(id));
                    ccf = meanccf(DATA,id,ic(j),ic(k));
                    si = SyncIndices(ccf);
                    y(nc) = si(1);
                    z(nc) = mean(xcs(id));
                    if strcmp(plottype,'xcorrb')
                        plot(mean(ampratio(id)), mean(synci(id)),'o','color',c,'buttondownfcn',{@HitXcorr, id, 0, pairs(id(1),:)});
                        hold on;
                        subplot(1,2,2);
                        plot(mean(ampratio(id)), mean(xcs(id)),'o','color',c,'buttondownfcn',{@HitXcorr, id, 0, pairs(id(1),:)});
                        hold on;
                    elseif strcmp(plottype,'xcorrd')
                        plot(mean(xcs(id)), y(nc),'o','color',c,'buttondownfcn',{@HitXcorr, id, 0, pairs(id(1),:)});
                        hold on;
                        subplot(1,2,2);
                        plot(si(2), si(1),'o','color',c,'buttondownfcn',{@HitXcorr, id, 0, pairs(id(1),:)});
                        hold on;
                    end
                end
                end
            end
            subplot(1,2,1);
            if strcmp(plottype,'xcorrb')
                set(gca,'xscale','log');
                ylabel('syncii');
            else
                xlabel('shape xcorr');
                ylabel('synci');
            end
            set(gca,'UserData',DATA.toplevel);
            subplot(1,2,2);
            if strcmp(plottype,'xcorrb')
                set(gca,'xscale','log');
                xlabel('Amplitude ratio');
                ylabel('shape xcorr');
            else
                xlabel('allsync');
            end
            set(gca,'UserData',DATA.toplevel);
            id  = find(~isnan(x) & ~isnan(y) & ~isinf(y));
            X(:,1) = log(x(id));
            X(:,2) = y(id);
            X(:,3) = z(id);
            pairim = zeros(DATA.nprobes, DATA.nprobes);
            pairimy = pairim;
            pairimz = pairim;
            pairim = pairim + mean(log(x));
            [Eig, Ev] = eig(cov(X));
            pcs = X * Eig;
            SetFigure(DATA,DATA.tag.pairs);
            subplot(2,2,1);
            hold off;
            for j = 1:nc
                pairim(mpairs(j,1),mpairs(j,2)) = log(x(j));
                pairimy(mpairs(j,1),mpairs(j,2)) = log(1./y(j));
                pairimz(mpairs(j,1),mpairs(j,2)) = z(j);
            end
            
            imagesc(pairim,'buttondownfcn',{@HitImage, 2});
            set(gca,'UserData',DATA.toplevel);
            subplot(2,2,2);
            hold off;
            imagesc(pairimy,'buttondownfcn',{@HitImage, 2});
            set(gca,'UserData',DATA.toplevel);
            subplot(2,2,3);
            hold off;
            imagesc(pairimz,'buttondownfcn',{@HitImage, 2});
            set(gca,'UserData',DATA.toplevel);
            subplot(2,2,4);


           SetFigure(DATA, DATA.tag.all,'front');     
            E = AddEllipse(gcf,'wait','color','r','line');
            SetFigure(DATA,DATA.tag.pairs);
            plot(pcs(:,end),pcs(:,2),'.');
            if ~isempty(E)
            theta = atan((log(E.pos(3)./E.pos(1)))./(E.pos(4)-E.pos(2)));
            xy = xyrotate(log(x),y,theta);
            b = xyrotate(log(E.pos([1 3])),E.pos([2 4]),theta);
            id = find(xy(:,1) < mean(b([1 2])));
            hold on;
            plot(pcs(id,end),pcs(id,2),'r.');
            pairim = zeros(DATA.nprobes, DATA.nprobes);
            for j = 1:length(id)
                pairim(mpairs(id(j),1),mpairs(id(j),2)) = b(1)-xy(id(j),1);
            end
            SetFigure(DATA,DATA.tag.pairs);
            subplot(2,2,4);
            hold off;
            imagesc(pairim,'buttondownfcn',{@HitImage, 2});
            set(gca,'UserData',DATA.toplevel);
            end
        else
            for j = 1:length(allsync)
                plot(ns(j), allsync(j),'o','buttondownfcn',{@HitXcorr, j, allexpt(j), allcells(:,j)});
                hold on;
            end
            xlabel('Nspikes');
            ylabel('Synci (all)');
        end
        DATA.allexpt = allexpt;
        DATA.allsync = allsync;
        DATA.AllPairs = pairs;
        set(gca,'UserData',DATA.toplevel);
    end
    set(gca,'UserData',DATA.toplevel);
    set(DATA.toplevel,'UserData',DATA);

    
function DATA = ReadTemplateResults(DATA, nc)

for j = 1:length(Clusters)
    for k = 1:length(Clusters{j})
        DATA.tmplxcs(j,k) = Clusters{j}{k}.cluster{nc}.templatexc;
        DATA.tmplscores{j}(k,1:length(Clusters{j}{k}.cluster{nc}.xy)) = Clusters{j}{k}.cluster{nc}.xy(:,3);
        DATA.mahal(j,k,1) = Clusters{j}{k}.cluster{nc}.gmdprime;
        DATA.mahal(j,k,2) = Clusters{j}{k}.cluster{nc}.mahal(1);
    end
end




function bad = CheckFitDim(C)
    bad = 0;
   if isfield(C,'gmfit') && isfield(C.gmfit,'mu')     
       nd = size(C.gmfit.mu,2);
       if C.space(1) == 6 && (C.space(2) == 1 && nd ~= 4 || ...
                C.space(2) == 4 && nd ~= 6 || ...
                C.space(2) == 3&& nd ~= 4)
            if C.auto
                astr = '(auto)';
            else
                astr = [];
            end
                fprintf('E%.0fP%d Space ND%d, but Fit has %d dimensions %s \n',C.exptno,C.probe(1),C.space(2),nd,astr);
                bad = 1;
       end
   elseif ~isfield(C,'space')
       cprintf('errors','E%dP%d missing space\n',C.exptno,C.probe(1));
   end



function [DATA, Clusters] = ReadClusterResults(DATA, Clusters)
if DATA.datatype == 2
    DATA = ReadTemplateResults(DATA, DATA.templatesrc);
    return;
end

if DATA.useautoclusters
    AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
end

ctimes = CellToMat(Clusters{1},'ctime');
lastc.ctime = max(ctimes);
for j = 1:length(Clusters)
    %first make sure required fields are present
    for k = 1:length(Clusters{j})
        if ~isfield(Clusters{j}{k},'mahal') || length(Clusters{j}{k}.mahal) < 4 
            Clusters{j}{k}.mahal(4) =0;
        end
        if ~isfield(Clusters{j}{k},'quick')
            Clusters{j}{k}.quick =0;
        end
        if ~isfield(Clusters{j}{k},'MeanSpike')
            Clusters{j}{k}.MeanSpike = [];
        end
            
        if ~isfield(Clusters{j}{k},'Evec')
            Clusters{j}{k}.Evec.Eval = NaN;
        end
        if ~isfield(Clusters{j}{k},'probe')
            Clusters{j}{k}.probe = k;
        end
        if ~isfield(Clusters{j}{k},'bsetspace')
            Clusters{j}{k}.bestspace = [NaN 0];
        end

        if isfield(Clusters{j}{k},'gmfit') && isempty(Clusters{j}{k}.gmfit)
            Clusters{j}{k} = rmfield(Clusters{j}{k},'gmfit');
        end
        Clusters{j}{k}.exptid = j;
        for c = 1:length(Clusters{j}{k}.next)
            if isfield(Clusters{j}{k}.next{c},'space')
                Clusters{j}{k}.next{c}.exptid = j;
            end
        end
        if ~isfield(Clusters{j}{k},'space')
            Clusters{j}{k}.quick = NaN;
        end
        if DATA.useautoclusters
            AutoClusters{j}{k}.exptid = j;
        end
    end
    %then read values wanted in data
    DATA = ReadClustersIntoData(DATA,Clusters{j},j);
end


eid = CellToMat(Clusters,'exptid');
exptno = CellToMat(Clusters,'exptno');
exptno = exptno(:,1);
[a,b] = find(eid ==0);
for j = 1:length(a)
    DATA = AddError(DATA,'Cluster E%d(id%d)P%d exptid = 0!!!!!',exptno(a(j)),a(j),b(j));
    Clusters{a(j)}{b(j)}.exptid = a(j);
end
pid = CellToMat(Clusters,'probe');
[a,b] = find(pid ==0);
for j = 1:length(a)
    fprintf('Cluster E%dP%d probe = 0!!!!!\n',a(j),b(j));
    Clusters{a(j)}{b(j)}.probe = b(j);
end


setappdata(DATA.toplevel,'Clusters',Clusters);
if DATA.useautoclusters
    setappdata(DATA.toplevel,'AutoClusters',AutoClusters);
end
   
Expts = getappdata(DATA.toplevel,'Expts');
Expts = SetExptTimeOffset(Expts);
setappdata(DATA.toplevel,'Expts',Expts);

function DATA = AddError(DATA, varargin)
    s = sprintf(varargin{:});
    mycprintf('red',s);
    fprintf('\n');
    if ~isfield(DATA,'errs')
        DATA.errs = {};
    end
    DATA.errs = {DATA.errs{:} s};
    
    
function Expts = SetExptTimeOffset(Expts)
    if length(Expts)
        for j = 1:length(Expts)
            if isfield(Expts{j},'Header')
                ctimes(j) = Expts{j}.Header.CreationDate;
            end
        end
        estart = min(ctimes(ctimes > 0));
        for j = 1:size(Expts,1)
            if isfield(Expts{j},'Header')
                Expts{j}.Header.timeoffset = (ctimes(j)-estart).*24.*60.*60;
            end
        end
    end


function DATA = ReadClustersIntoData(DATA, Clusters, exid)
    
    j = exid;
    
    for k = 1:length(Clusters)
        C = Clusters{k};
    if isfield(C,'space')
        if C.space(1) == 6
            CheckFitDim(C);
        end

        if size(C.MeanSpike.ms,1) > 1
            A = std(C.MeanSpike.mu');
            DATA.muvar(j,k) = A(k);
            A = std(C.MeanSpike.ms');
            [a,b] = max(A);
            DATA.peakdiff(j,k) = b-k;
            DATA.mspkvar(j,k) = A(k);
            xi = [1:0.1:length(A)];
            y = InterpGauss1([1:length(A)],A,xi,1.5);
            [a,b] = max(y);
            DATA.peakpos(j,k) = xi(b);
        else
            DATA.mspkvar(j,k) = std(C.MeanSpike.ms);
            DATA.muvar(j,k) = std(C.MeanSpike.mu);
            DATA.peakpos(j,k) = k;
            DATA.peakdiff(j,k) = 0;
        end
        DATA.mahal(j,k,1) = C.mahal(1); %2-D
        if DATA.useautoclusters
            DATA.automahal(j,k,1) = AutoClusters{j}{k}.mahal(1); %2-D
            DATA.automahal(j,k,3) = AutoClusters{j}{k}.mahal(4); %1-D
            DATA.automahal(j,k,4) = AutoClusters{j}{k}.mahal(2);
            DATA.automahal(j,k,2) = AutoClusters{j}{k}.mahal(2);
        end
        if length(C.mahal) > 3
            DATA.mahal(j,k,3) = C.mahal(4); %1-D
        end
        if isfield(C,'clst')
            DATA.nclusters(j,k) = max(C.clst)-1;
        else
            DATA.nclusters(j,k) = 1;
        end

        if isfield(C,'pcgms')
            DATA.pcgms(j,k,1:length(C.pcgms)) = C.pcgms;
        end
        DATA.mahal(j,k,4) = C.mahal(2); %N-D
        if isfield(C,'bestspace')
            DATA.mahal(j,k,2) = C.bestspace(1); %N-D
            if isfield(C,'bestd')
                DATA.bestd(j,k,:) = C.bestd;
            end
        else
            DATA.mahal(j,k,4) = C.mahal(2); %N-D
        end
        if C.space(1) == 6
            DATA.cutspace(j,k) = C.space(2);
        else
            DATA.cutspace(j,k) = 5+C.space(1);
        end
        DATA.cutshape(j,k) = C.shape;
        if isfield(C,'Evec') && isfield(C.Evec,'Eval')
            DATA.eveci(j,k) = C.Evec.Eval(1)./sum(C.Evec.Eval);
        end
        xc = corrcoef(C.MeanSpike.ms(:),C.MeanSpike.mu(:));
        DATA.mucorrs(j,k) = xc(1,2);
        xc = (C.MeanSpike.ms(:)' * C.MeanSpike.mu(:)) ./ (C.MeanSpike.ms(:)'*C.MeanSpike.ms(:));
        DATA.muampl(j,k) = xc;
        if isfield(C,'dpsum')
        DATA.dpsum(j,k) = sum(C.dpsum);
        end
        if isfield(C,'gmfit')
            DATA.fitres(j,k) = C.gmfit.Converged;
        end
        if isfield(C,'gmfit1d') && isfield(C.gmfit1d,'Converged')
            DATA.fitres1d(j,k) = C.gmfit1d.Converged;
        end
        if C.space(1) == 6  %an N-D mahal fit;
            DATA.fitspace(j,k) = 10+C.space(2);
        else
            DATA.fitspace(j,k) = C.space(1);
        end
        DATA.pcspace(j,k) = C.dvdt + C.csd.*2;
        DATA.ctimes(j,k,2) = C.ctime;
        DATA.isauto(j, k) = C.auto;
        if isfield(C,'manual')
            DATA.isauto(j,k) = 3;
        end
        if isfield(C,'starttime') && C.auto ==1
            DATA.ctimes(j,k,1) = C.starttime;
        else
            DATA.ctimes(j,k,1) = C.ctime;
        end
        if isfield(C,'first')
            F = C.first;
            ndive = 1;
            while(isfield(F,'first'))
                F = F.first;
                ndive = ndive+1;
            end
            if DATA.ctimes(j,k,1)-F.starttime > 0.2
                fprintf('Very long firsttime Ex%d, P%d (%d) %.2f %s vs %s\n',...
                    j,k,C.auto,DATA.ctimes(j,k,1)-F.starttime,datestr(DATA.ctimes(j,k,1)),datestr(F.starttime));
            else
                DATA.ctimes(j,k,1) = F.starttime;
            end
            DATA.thriter(j,k) = ndive;
            %            DATA.ctimes(j,k,1) = F.starttime;
        else
            DATA.thriter(j,k) = 0;
        end
        if isfield(C,'savetime')
            DATA.ctimes(j,k,3) = C.savetime(end);
        else
            DATA.ctimes(j,k,3) = lastc.ctime;
        end
        if isfield(C,'fitdprime')
            DATA.GaussFitdp(j,k,2) = C.fitdprime(1);
            if isfield(C,'gmdprime')
                DATA.GaussFitdp(j,k,1) = C.gmdprime;
            else
                DATA.GaussFitdp(j,k,1) = C.mahal(1);
            end
            if length(C.fitdprime) >= 3
                DATA.gmfitpos(j,k,:) =C.fitdprime(2:3);
%now the value is signed, don't ned this check for overlap.                
                if sum(C.fitdprime(2:3)) < 2
%                    DATA.GaussFitdp(j,k,2) = 0;
                end
            end
        elseif isfield(DATA,'GaussFitdp') && j <= size(DATA.GaussFitdp,1) && k <= size(DATA.GaussFitdp,2)
% this iss only of old files where this was calculated in plotclusters,
% after allvpcs
%            Clusters{j}{k}.fitdprime(1) = DATA.GaussFitdp(j,k,2);
        end
        DATA.nevents(j,k) = C.nspks;
        if isfield(C,'clst')
        DATA.nspks(j,k) = sum(C.clst > 1);
        end
        DATA.dropi(j,k,1) = C.dropi(3);
        for c = 1:length(C.next)
            if isfield(C.next{c},'dropi')
                DATA.dropi(j,k,c+1) = C.next{c}.dropi(3);
            end
        end
        if ~isfield(C,'clusterprog') || isempty(C.clusterprog)
            DATA.clusterprog(j,k) = 0;
        else
            DATA.clusterprog(j,k) = strmatch(C.clusterprog(1:7),{'AllVPcs' 'PlotClusters'});
        end
    end

    lastc = C;
end

function CloseTag(tag)
it = findobj('Tag',tag,'type','Figure');
for j = 1:length(it)
    close(it(j));
end

function LoadAllSpikes(DATA, varargin)
    
    verbose = 1;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'quiet',5)
            verbose = 0;
        end
        j = j+1;
    end
    ts = now;
    [a,b] = fileparts(DATA.name);
    mkn = GetMonkeyName(DATA.name);
    if 0
    b = regexprep(b,'[GM]([0-9]*)(.*)','M$1');
    else
    b = regexprep(b,'[GM]([0-9]*)(.*)','$0');
    end
    [c,d] = fileparts(a);
    for k = 1:length(DATA.exptid)
        for p = 1:DATA.nprobes
            xs = '';
            if rem(DATA.exptid(k),1) > 0.001
                xs = 'a';
            end
            name = [DATA.name '/Spikes/' mnk b '.p' num2str(p)  't' num2str(floor(DATA.exptid(k))) xs '.mat'];
            if verbose > 0
            fprintf('Reading %s at %s\n',name,datestr(now));
            end
            AllSpikes{k,p} = ReadSpikeFile(name);
            AllSpikes{k,p}.probe = p;
            if isfield(AllSpikes{k,p},'Header')
            if ~isfield(AllSpikes{k,p}.Header,'ctime')
                AllSpikes{k,p}.Header.ctime = 0;
            end
            else
                AllSpikes{k,p}.Header.ctime = 0;
            end
        end
    end
    fprintf('Spike Load took %.2f\n',mytoc(ts));
    setappdata(DATA.toplevel,'AllSpikes',AllSpikes);

function LoadSelectedSpikes(DATA,eid,pid)
    ts = now;
    [a,b] = fileparts(DATA.name);
    mnk = GetMonkeyName(DATA.name);
    [c,d] = fileparts(a);
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    for j = 1:length(eid);
        e = eid(j);
        p = pid(j);
        xs = '';
        if rem(DATA.exptid(e),1) > 0.001
            xs = 'a';
        end
        name = [DATA.name '/Spikes/' mnk b '.p' num2str(p)  't' num2str(floor(DATA.exptid(e))) xs '.mat'];
        fprintf('Reading %s at %s\n',name,datestr(now));
        AllSpikes{e,p} = ReadSpikeFile(name);
        AllSpikes{e,p}.probe = p;
    end
    fprintf('Spike Load took %.2f\n',mytoc(ts));
    setappdata(DATA.toplevel,'AllSpikes',AllSpikes);

    
function [h, details] = QuickSpikes(DATA, pos, varargin)

    scalebar = 1;
    if isstruct(pos)
        Spks = pos;
        C = varargin{1};
        e = find(DATA.exptid == C.exptno) ;
        p = C.probe(1);
        pos = [e p];
        if size(C.MeanSpike.ms,1) == 1
            mprobe = 1;
        else
            mprobe = p;
        end
        if p > length(DATA.voffset)
            DATA.voffset(p) = 0;
        end
    else
        e = pos(1);
        p = pos(2);
        AllSpikes = CheckAllSpikes(DATA, e, p);
        Clusters = getappdata(DATA.toplevel,'Clusters');
        C = Clusters{e}{p};
        Spks = AllSpikes{e,p};
        maxv = CellToMat(AllSpikes(e,:),'VRange');
        DATA.voffset =  cumsum(cat(1,maxv(2:end,2), maxv(end,2))) - cumsum(maxv(:,1));
        DATA.voffset = DATA.voffset * 0.7;
        mprobe = pos(2);
    end
    
    showmean = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'scalebar',6)
            j = j+1;
            scalebar = varargin{j};
        elseif strncmpi(varargin{j},'showmean',8)
            showmean = 1;
        end
        j = j+1;
    end
    details = [];
    if ~isfield(Spks,'times') %will happen if this file is missing
        h = 0;
        return;
    end
    if isfield(Spks,'maxv') && strcmp(class(Spks(1).values),'int16')
        for j = 1:length(Spks)
        Spks(j).values = double(Spks(j).values) .* double(Spks(j).maxv./Spks(j).maxint);
        end
    end
    nevents = length(Spks(1).times);
    if nevents > 50000
        nspks = 1000;
    else
        nspks = 500;
    end
    istep = max([1 round(nevents/nspks)]);
    ispk = 1:istep:nevents;
    if length(C.clst) >= nevents && sum(C.clst(ispk) > 1) == 0 && sum(C.clst > 1) > 0 %no classified spikes
        xspk = find(C.clst > 1);
        ispk = [ispk xspk(1) xspk(end)];
    end
    details.ispk = ispk;
    DATA.usegmcid = 0;
    DATA.plotspk.bytrial = 0;
    yl = minmax(Spks(1).values(:));
    voff = DATA.voffset - DATA.voffset(p);
    if length(Spks) == 3
        yl(2) = voff(Spks(3).probe) + Spks(3).VRange(2);
        yl(1) = voff(Spks(2).probe) + Spks(2).VRange(1);
    end
    h = PlotSpikes(DATA, pos, ispk, Spks, C, 'fixy',yl);
    id = find(h == 0); %no spikes plotted
    for j = id
        ispk = find(C.clst == j);
        if length(ispk)
        ispk = ispk([1 end]);
        x = PlotSpikes(DATA, pos, ispk, Spks, C, 'fixy',yl,'holdon');
        h(j) = x(j);
        end
    end
    hold on;
    if showmean
        plot(C.MeanSpike.mu(mprobe,:),'-','color','k');
        for j = 2:length(h)
            if h(j) > 0 & ishandle(h)
                color = get(h(j),'color');
                color = 1-color;
                if j > 2
                    if isfield(C.next{j-2},'MeanSpike')
                    details.meanh(j) = plot(C.next{j-2}.MeanSpike.ms(mprobe,:),'-','color',color,'linewidth',2);
                    end
                else
                    details.meanh(j) = plot(C.MeanSpike.ms(mprobe,:),'-','color',color,'linewidth',2);
                end
            end
        end
        if length(h) < 2 
            details.meanh(2) = plot(C.MeanSpike.ms(mprobe,:),'-','color',DATA.colors{2},'linewidth',2);
        end
    end
    plot([1 1],[0 -scalebar],'k-','linewidth',2); %scale bar
    axdata = get(gca,'UserData');
    axdata.toplevel = DATA.toplevel;
    axdata.probe = p;
    axdata.eid = e;
    set(gca,'UserData',axdata);


function stopped = SpoolAllProbes(DATA, e, Spikes, Clusters, varargin)
    
    usetrials = [];
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'Trials',5)
            j = j+1;
            usetrials = varargin{j};
        end
        j = j+1;
    end
    DATA.usegmcid = 0;
    DATA.plotspk.bytrial = 1;
    DATA.voffset = [1:DATA.nprobes] .*2;
    Expts = getappdata(DATA.toplevel,'Expts');
    DATA.Expt = Expts{e,1};
    stopped = 0;
    
    chspk = [];
    for j = 1:length(chspk)
        if ~isempty(AllSpikes{e,chspk(j)})
            Spks(j+1) = AllSpikes{e,chspk(j)};
        end
    end
    for j = 1:length(chspk)
        Spks(j+1).probe = chspk(j);
    end

   F =  SetFigure(DATA, DATA.tag.allspikes);
   %subplot(1,1,1);
   stopbtn = findobj(F, 'tag','StopSpool');
   for k = 1:length(Spikes)
       p = Spikes{k}.probe;
       if isfield(Spikes{k},'VRange')
           ranges(k,:) = Spikes{k}.VRange;
       else
           ranges(k,:) = minmax(Spikes{k}.values(:));
       end
       if DATA.plotspk.showexcluded
           xcls{k} = [];
       else
       xcls{k} = FindExcludedTrials(DATA,e,p,1,Clusters{p});
       end
   end

   tic; 
   if length(usetrials)
       Trials = Expts{e,1}.Trials(usetrials);
   else
       Trials = Expts{e,1}.Trials;
   end
    xcl = [];
    if length(Spikes) == 24
        nr = 6;
        nc=4;
    elseif length(Spikes) == 96
        nr = 10;
        nc=10;
        plotpos = SetPlotPos(DATA, length(Spikes),nr,nc);
    else
        [nr,nc] = Nsubplots(length(Spikes));
    end

    if DATA.profiling
        profile on;
    end
    
    c = get(gcf,'Children');
    types = get(c,'type');
    axs = c(strcmp('axes',types));
    for j = 1:length(Trials)
        set(F,'name',sprintf('Ex%.0f Trial %d',DATA.exptid(e),Trials(j).Trial));
        for k = 1:length(Spikes)
            p = Spikes{k}.probe;
            chspk = [p-1 p+1];
            chspk = chspk(chspk > 0 & chspk <= DATA.nprobes);
            subplot(axs(p));
            %mysubplot(nr,nc,p,'tight');
            C = Clusters{p};
            if ~ismember(j, xcls{k})
                DATA.currenttrial = j;
                if DATA.usesavedcodes
                ispk = find(Spikes{k}.times > Trials(j).Start(1) & Spikes{k}.times < Trials(j).End(end));
                else
                ispk = find(C.times > Trials(j).Start(1)./10000 & C.times < Trials(j).End(end)./10000);
                end
                PlotSpikes(DATA, [e p], ispk, Spikes{k}, C, 'fixy',ranges(k,:),'empty');
                [a,b] = isacell(DATA,e,p);
                if a
                    xl = get(gca,'xlim');
                    yl = get(gca,'ylim');
                    if length(xcls{k})
                        str = sprintf('%d *',b);
                    else
                        str = sprintf('%d ',b);
                    end
                    set(gca,'xcolor','b','ycolor','b','linewidth',2);
                    text(mean(xl), yl(2), sprintf('Cell %s',str),'VerticalAlignment','top','color','b','fontweight','bold');
                end
                set(gca,'ButtonDownFcn',{@HitXYPlot, e,p});
            end
        end
        drawnow;
        stopped = get(stopbtn,'value');
        if stopped
            set(stopbtn,'value',0);
            return;
        end
    end

toc
    if DATA.profiling
        profile viewer;
    end


function PlotSyncSpikes(DATA, eid, probes, clnum, varargin)
j = 1;
if DATA.plot.synctmax > 0
    tlim = DATA.plot.synctmax .* 10;
else
    tlim = 20;
end
    
nstep = 10;
while j <= length(varargin)
    if strncmpi(varargin{j},'flip',4)
        probes = fliplr(probes);
    elseif strncmpi(varargin{j},'tmax',4)
        j = j+1;
        tlim = varargin{j}(1);
        if length(varargin{j}) > 1
            nstep = round(varargin{j}(2)/2);
        end
    end
    j = j+1;
end
    
 AllSpikes = CheckAllSpikes(DATA, eid, probes);
 Clusters = getappdata(DATA.toplevel,'Clusters');
 
 atid = find(Clusters{eid}{probes(1)}.clst == clnum(1)+1);
 btid = find(Clusters{eid}{probes(2)}.clst == clnum(2)+1);
 if isempty(AllSpikes)
     xct = -(tlim/10000):0.0002:(tlim/10000);
     [xc, b] = xcorrtimes(Clusters{eid}{probes(1)}.times(atid),Clusters{eid}{probes(2)}.times(btid),'times',xct);
     SetFigure(DATA,DATA.tag.spikes);
     plot(b.xpts, xc);
     return;
 end
 ta = AllSpikes{eid,probes(1)}.times(atid);
 tb = AllSpikes{eid,probes(2)}.times(btid);
 [a, aid, bid]= intersect(round(ta/10),round(tb/10));
% bid = find(ismember(round(ta/10),round(tb+10/10)));
 
 
aid = [];
bid = [];
dts = [];
for j = 1:length(ta)
    dt = ta(j)-tb;
%find all events witin tlime and show them    
    id = find(abs(dt) < tlim)';
    if length(id)
        aid = [aid j * ones(size(id))];
        bid = [bid id];
        dts = [dts dt(id)'];
    end
end
% id = union(aid,bid);
voff = AllSpikes{eid,probes(1)}.maxv .* AllSpikes{eid,probes(1)}.maxint/2;
x = [1:size(AllSpikes{eid, probes(1)}.values,2)]./40;
hold off;
 for j = 1:length(aid)
     %+v toff = 1 after 2, so subtract from 2
     toff = AllSpikes{eid,probes(1)}.times(atid(aid(j)))-AllSpikes{eid,probes(2)}.times(btid(bid(j)));
     toff = toff ./10;
     plot(x,AllSpikes{eid,probes(1)}.values(atid(aid(j)),:),'k');
     hold on;
     plot(x-toff,double(AllSpikes{eid,probes(2)}.values(btid(bid(j)),:))+voff,'b');
 end
 text(x(end),voff,sprintf('%d/%d',probes(2),clnum(2)),'horizontalalignment','left','color','b','fontsize',20);
 text(x(end),0,sprintf('%d/%d',probes(1),clnum(1)),'horizontalalignment','left','color','k','fontsize',20);
 yl = get(gca,'ylim');
 xl = get(gca,'xlim');
 if tlim <= 50
     [y,x] = hist(dts,-tlim:2:tlim);
 else
     [y,x] = hist(dts,-tlim:5:tlim);
 end
pmax = max(y)./length(atid);
y = y .* yl(2)./max(y);
x  = -x./10;
title(sprintf('%d events (of %d) within +- %.2fms pmax %.3f',length(aid),length(atid),tlim./10,pmax));
%x = (x-min(x)) .* xl(2)./max(x);
%x = x+xl(1);
plot(x,y,'ro-');

 function voffset = SetVOffset(DATA, AllSpikes, e)
     if isempty(AllSpikes)
         voffset = [1:DATA.nprobes] .* 2; %? could use meanspiek - get values at startup
     else
     maxv = CellToMat(AllSpikes(e,:),'VRange');
     if isempty(maxv)
         voffset = [1:DATA.nprobes] .* 2; %? could use meanspiek - get values at startup
     else
     voffset =  cumsum(cat(1,maxv(2:end,2), maxv(end,2))) - cumsum(maxv(:,1));
     voffset = voffset * 0.7;
     end
     end


function [stopped, h] = SpoolSpikes(DATA, pos, varargin)

    useids = [];
    args = {};
    setlist = 1;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'ids',3)
            j = j+1;
            useids = varargin{j};
        else
            args = {args{:} varargin{j}};
        end
        j = j+1;
    end
    
    stopped = 0;
    e = pos(1);
    p = pos(2);
    DATA.usegmcid = 0;
    DATA.plotspk.bytrial = 1;
    DATA.voffset = [1:DATA.nprobes] .*2;
    Expts = getappdata(DATA.toplevel,'Expts');
    DATA.Expt = Expts{e,1};
    Trials = Expts{e,1}.Trials;
    Clusters = getappdata(DATA.toplevel,'Clusters');
    C = Clusters{e}{p};
    chspk = [p-DATA.nprobespool:p+DATA.nprobespool];
    chspk = chspk(chspk > 0 & chspk <= DATA.nprobes & chspk ~= p);
    AllSpikes = CheckAllSpikes(DATA, e, [chspk p]);
    DATA.voffset = SetVOffset(DATA, AllSpikes, e);
    maxv = CellToMat(AllSpikes(e,:),'VRange');

    voff = (DATA.voffset - DATA.voffset(p));
    pa = min([chspk p]);
    pb = max([chspk p]);
    yl(1) = voff(pa)+maxv(pa,1);
    yl(2) = voff(pb)+maxv(pb,2);
    Spks(1) = AllSpikes{e,p};
    for j = 1:length(chspk)
        if ~isempty(AllSpikes{e,chspk(j)})
            Spks = [Spks AllSpikes{e,chspk(j)}];
        end
    end
    for j = 1:length(chspk)
        Spks(j+1).probe = chspk(j);
    end

    F = SetFigure(DATA, DATA.tag.spikes,'front');
    stopbtn = findobj(F,'Tag','StopSpool');
    tic;
    if length(useids)
        Trials = Trials(useids);
        firsttrial = 1;
    else
        mint = min(cat(1,Spks.times));
        useids = find([Trials.TrialStart] > mint-10000);
        firsttrial = useids(1);
        useids = 1:length(Trials);
    end
    xcl = FindExcludedTrials(DATA,e,p,1,C);
    for j = firsttrial:length(Trials)
        if ~ismember(j, xcl) || length(useids) == 1
        DATA.currenttrial = useids(j);
        ispk = find(C.times > Trials(j).Start(1)./10000 & C.times < Trials(j).End(end)./10000);
        h = PlotSpikes(DATA, pos, ispk, Spks, C, 'fixy',yl,args{:});
        drawnow;
        end
        stopped = get(stopbtn,'value');
        if stopped
            set(stopbtn,'value',0);
            stopped = 1;
            return;
        end
    end
    toc;
    clear AllSpikes;
    clear Clusters;
    if length(useids) ==1 || setlist
        SetTrialList(DATA, C, useids);
    end
    set(DATA.toplevel,'UserData',DATA);
    
function ploth = PlotSpikes(DATA, pos, spkid, Spks, C, varargin)

dvdt = 1;
ploth = [];
fixy = 0;
colors {1} = [0.5 0.5 0.5];
colors {2} = [1 0 0];
colors {3} = [0 1 0];
colors = DATA.colors;
scale = 1;
yl = [];
j = 1;
quicktest = 1;
showax = 1;
showtitle = 1;
holdoff = 1;
showall = 0;
e = pos(1);
p = pos(2);

while j <= length(varargin)
    if strncmpi(varargin{j},'colors',6)
        j = j+1;
        colors = varargin{j};
    elseif strncmpi(varargin{j},'fixy',3)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            yl = varargin{j};
        else
            fixy = 1;
        end
    elseif strncmpi(varargin{j},'empty',5)
        showax = 0;
        showtitle = 0;
    elseif strncmpi(varargin{j},'holdon',5)
        holdoff = 0;
    elseif strncmpi(varargin{j},'showall',7)
        showall = 1;
    end
    j = j+1;
end

if ~isfield(Spks,'values')
    fprintf('Missing Spikes for Expt %d\n',e);
    return;
end
for j = 1:length(Spks)
    if isfield(Spks,'maxv') && strcmp(class(Spks(j).values),'int16')
        Spks(j).values = double(Spks(j).values) .* double(Spks(j).maxv)./16000;
    end
end
if fixy
    yl = get(gca,'ylim');
end
chspk = pos(2);

if isempty(spkid)
    for c= chspk;
        if c > 0 & c <= DATA.nprobes
        for j = 1:2
            h = plot([0 30],[0 0],'color',colors{j});
        end
        hold on;
        end
    end
    if fixy
        set(gca,'ylim',yl);
    end
    if DATA.plotspk.bytrial && showtitle
        nt = DATA.currenttrial;
        T = DATA.Expt.Trials(nt);
        title(sprintf('Ex%d T%d %.2f-%.2f: No Spikes  ed%.2f',e,T.Trial,T.Start(1),T.End(end),T.ed));
    end
    return;
end

if spkid(end) > size(Spks(1).values,1)
    fprintf('Spike Events Mismatch E%dP%d\n',e,p);
    spkid = spkid(spkid <= size(Spks(1).values,1));
end

if length(spkid) == 0
    return;
end

if isfield(C,'clst') && length(C.clst) >= max(spkid)   
    clst = C.clst;
    if showall
        clst(clst < 1) = 1;
    end
nc = unique(clst(spkid));
nc = nc(nc > 0); %excludes excluded trials/times


twoclusters = 1;
else
    twoclusters = 0;
end
if isempty(spkid)
    ids{2} = DATA.clid;
    ids{1} = DATA.nid;
    id = ids{2};
   nid = ids{1};
elseif DATA.usegmcid && length(DATA.gmcid) >= max(spkid)
        nc = unique(DATA.gmcid);
        for j = 1:length(nc)
            id = find(DATA.gmcid(spkid) == nc(j));
            ids{j} = spkid(id);
        end
elseif DATA.usesavedcodes
    nc = unique(Spks(1).codes(spkid,1));
    for j = 1:length(nc)
        ids{nc(j)+1} = spkid(find(Spks(1).codes(spkid,1) == nc(j)));
    end
    id = ids{end};
elseif twoclusters
    ids = {[] []}; %min necessary
    for j = 1:length(nc)
    ids{nc(j)} = spkid(find(C.clst(spkid) == nc(j)));
    end
    id = ids{2};
   nid = ids{1};
else
    nc = unique(Spks(1).codes(spkid,1)+1);
    nc = nc(nc > 0);
    for j = 1:length(nc)
        ids{nc(j)} = spkid(find(Spks(1).codes(spkid,1) == nc(j)-1));
    end
    id = ids{2};
end
ispk = pos(2);
for j = 1:length(ids)
    V{j} = Spks(1).values(ids{j},:)';
end

nspks = length(id);
voff = DATA.voffset - DATA.voffset(ispk(1));

    l = size(V{1},1);
    if holdoff
        hold off;
    else
        hold on;
    end
    x = [1:l NaN];
if isfield(Spks(1),'Vrange')
    dv = max(Spks(1).Vrange)/3;
else
    dv = 0;
end
for c= chspk;
    if c > 0 & c <= DATA.nprobes
        for j = 1:length(V)
        nV = V{j} + voff(c);
        if length(nV) == 0
        elseif size(nV,1) > 1 
            ploth(j) = plot([0 30],[0 0],'color',colors{j});
            nV(l+1,:) = NaN;
            set(ploth(j),'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));
            hold on;
        else
            ploth(j) = plot(1:l,nV,'color',colors{j});
            hold on;
        end
        end
        h = text(size(V{1},1)-2,voff(c)+dv,sprintf('%d',c));
        set(h,'fontweight','bold','horizontalalignment','right');
    end
end
spkid = spkid(spkid <= length(C.clst));
for c = 2:length(Spks)
    if isfield(Spks(c),'Vrange')
        dv = max(Spks(c).Vrange)/3;
    else
        dv = abs(voff(Spks(c).probe)-voff(Spks(1).probe))/5;
    end
    %Spks.times is in 0.1ms ticks
    [ix, ida, idb] = intersect(round(Spks(c).times./10),round(Spks(1).times(spkid)./10));
    for j = 1:length(nc)
        id = find(C.clst(spkid(idb)) ==nc(j));
        V{j} = Spks(c).values(ida(id),:)';
        l = size(V{1},1);
        x = [1:l NaN];
    nV = V{j} + voff(Spks(c).probe);
    if size(nV,1) > 1 
        h = plot([0 30],[0 0],'color',colors{j});
        nV(l+1,:) = NaN;
        set(h,'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));
        hold on;
    end
    end
    h = text(size(V{1},1)-2,voff(Spks(c).probe)+dv,sprintf('%d',Spks(c).probe));
    set(h,'fontweight','bold','horizontalalignment','right');
end

if length(yl) == 2
    set(gca,'ylim',yl);
end
if showax == 0
    set(gca,'xtick',[],'ytick',[]);
end
hold off;

if showtitle
if DATA.plotspk.bytrial
    nt = max([DATA.currenttrial 0]);
    T = DATA.Expt.Trials(nt);
title(sprintf('E%d T%d %.2f-%.2f: %d/%d ed%.2f',e,T.Trial, T.Start(1)./10000,T.End(end)./10000, nspks,length(spkid),DATA.Expt.Trials(nt).ed));
else
title(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...
    spkid(1),spkid(end),Spks(1).times(spkid(1)),Spks(1).times(spkid(end)),nspks,length(spkid)));
end
end


function CheckRates(DATA, expname, cell)
    exid = find(strcmp(expname, DATA.expnames));
    nx = 0;
    for e = exid
        [p,clnum] = find(DATA.CellList(e,:,:) == cell);
        if ~isempty(p)
            nx = nx+1;
            Expts{nx} = CountExptSpikes(DATA, Expts{e},Clusters{e}{p},clnum);
        end
    end
    
    function R = CheckAllRateSequences(DATA)
        
        plottype = 2;
Expts = getappdata(DATA.toplevel,'Expts');
Clusters = getappdata(DATA.toplevel,'Clusters');
R.cellid = [];
if isempty(Clusters)
    return;
end
exname = unique(DATA.expnames);
Im = zeros(length(Expts),max(DATA.CellList(:)));
    cellid = unique(DATA.CellList);
    cellid = cellid(cellid > 0);
x = getappdata(DATA.toplevel,'RateCheckData');
if isfield(x,'plottype')
    R.plottype = x.plottype;
else
    R.plottype = 'cv';
end

if isempty(cellid)
    return;
end
for j = 1:length(exname)
    eid = find(strcmp(exname{j},DATA.expnames));
    for k = 1:length(cellid)
        c = cellid(k);
        E = {};
        for e = 1:length(eid)
            [p, cl] = find(squeeze(DATA.CellList(eid(e),:,:)) == c);
            if length(p) > 1
                fprintf('Cell %d Doubly Defined Expt %d P %s\n',c,eid(e),sprintf('%d ',p));
                p = p(1);
                cl - cl(1);
            end
            if ~isempty(p)
                E{e} = CountExptSpikes(DATA, Expts{eid(e)}, Clusters{eid(e)}{p},cl);
            end
        end
        if ~isempty(E)
        check = CheckExptRates(E,'print');
        R.ff(j,k) = max([check.errs.ff]);
        R.blkcv(j,k) = check.blkcv;
        R.blkff(j,k) = check.blkff;
        R.blkskew(j,k) = check.blkskew;
        R.warning(j,k) = sum(check.warning(2:end));
        R.slopes(j,k) = max(abs([check.errs.slope]));
        exid = [check.errs.exptno];
        Im(exid,c) = check.errs.ff;
        ffid = find([check.errs.ff] > DATA.checkrate.ff);
        Im(exid(ffid),c) = Im(exid(ffid),c)+1; 
        a = find(abs([check.errs.slope]) > DATA.checkrate.slope);
        Im(exid(a),c) = Im(exid(a),c)+2; 
        a = find(abs(log([check.diffs])) > DATA.checkrate.diff);
        if length(a) > 2 && sum(check.warning(2:end)) == 0
            a; %not useful? 
        end
        R.diffs(j,k) = max(abs(log([check.diffs])));
        Im(exid(a),c) = Im(exid(a),c)+4; 
        R.Im = Im;
        end
    end
end
R.cellid = cellid;
setappdata(DATA.toplevel,'RateCheckData',R)
      
SetFigure(DATA,DATA.tag.clusters);
PlotRateCheck(R);

function RatePlotMenu(a,b,fcn)
    DATA = GetDataFromFig(a);
    R = getappdata(DATA.toplevel,'RateCheckData');
    if sum(strcmp(fcn,{'cv' 'skew' 'image'}))
        R.plottype = fcn;
    end
    SetFigure(DATA,DATA.tag.clusters);
    PlotRateCheck(R);
    setappdata(DATA.toplevel,'RateCheckData',R);

function PlotRateCheck(R)
    hold off;
    if strcmp(R.plottype,'image');
        imagesc(R.Im,'buttondownfcn',{@HitImage, 'CellRates'});
    else
        sym = 'osx+<v>^*ph.osx+<v>^*ph.';
        colors = mycolors;
        if strcmp(R.plottype,'skew');
            x = R.blkff;
            y = R.blkskew;
        else
            x = R.blkff;
            y = R.blkcv;
        end
        hold off;
    for j = 1:size(R.blkcv,1);
        for k = 1:size(R.blkcv,2);
            h = plot(x(j,k),y(j,k),sym(j),'color',colors{k},'buttondownfcn',{@HitExptPoint,j,R.cellid(k)});
            if R.warning(j,k)
                set(h,'linewidth',2);
            end
            hold on;
        end
    end
    xlabel('Fano Factor (blockwise)');
    if strcmp(R.plottype,'skew');
        set(gca,'xscale','log','yscale','lin');
        ylabel('Skewness');
    else
        set(gca,'xscale','log','yscale','log');
        ylabel('SD/mean(blockwise)');
    end
end
fgf;


function HitExptPoint(a,b, ex, cell)
DATA = GetDataFromFig(a);
exname = unique(DATA.expnames);
fprintf('%s Cell%d\n',exname{ex},cell);
SetFigure(DATA,DATA.tag.rateseq);
PlotClusterRates(DATA, 'rateseqone','expt',exname{ex},'cell',cell);
%RateZoom(DATA,5,cell);

function res = CheckExptsRateSequence(Expts)
   
    if iscell(Expts)
        for j = 1:length(Expts)
            [res.errs(j), counts{j}] = CheckExptRateSequence(Expts{j});
            m(j) = mean(counts{j});
        end
        id = find(m > 0);
        mrate = mean(cat(2,counts{id}));
        medrate = median(m(m>0));
        sd = std(cat(2,counts{id}));
        res.vmr = sd./mrate;
        for j = 1:length(Expts)
            if m(j) > 0;
            res.diffs(j) = (res.errs(j).mean)./medrate;
            else
            res.diffs(j) = 1;
            end
        end
    elseif isstruct(Expts)
        res = CheckExptRateSequence(Expts);
    end
    
function [err, counts] = CheckExptRateSequence(Expt)
    counts = [Expt.Trials.count];
    if length(counts) > 100
        smw = 10;
    elseif length(counts) > 20
        smw = 5;
    elseif length(counts) == 0
        smw = 2;
    else
        smw = 2;
    end
    smc = smooth(counts, smw);
    err.ff = var(smc)./mean(smc);
    b = polyfit(1:length(counts),sqrt(counts),1);
    err.mean = mean(counts);
    err.slope = length(counts).*b(1)./mean(sqrt(counts));
    err.smw = smw;
    if isfield(Expt.Header,'exptno')
        err.exptno = Expt.Header.exptno;
    end
    if sum(counts) == 0
        fprintf('E%dCell%d no spikes\n',err.exptno,Expt.Header.cellnumber);
        err.slope = 0;
        
    end
    
    
    