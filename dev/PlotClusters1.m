function varargout = PlotClusters(name, varargin)
%PlotClusters(dir, ......)
%reads in Cluster files made by AllVPcs and plots up clusters for each
%expt/cell
prefsdir = '/bgc/group/matlab/preferences/PlotClusters';

TAGTOP = 'PlotClusters';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        TOPTAG = varargin{j};
    end
    j = j+1;
end
initcall =0 ;
if ishandle(name)  % callback
    gui = name;
    DATA = GetDataFromFig(name);
    name = varargin{2};
    if length(varargin) > 2
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
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
    end
    
elseif isdir(name)
    d = dir(name);
    strings = {};
    for j  = 1:length(d)
        if strfind(d(j).name,'ClusterTimes.mat') & isempty(strfind(d(j).name,'OldClusterTimes.mat'))
            strings = {strings{:} d(j).name};
        end
    end
    DATA.name = name;
    DATA.strings = strings;
    DATA.dates = zeros(size(DATA.strings));
    DATA.tag.top = TAGTOP;
    DATA = SetDefaults(DATA);
    DATA = InitInterface(DATA);
end

if strcmp(name,'getstate')
    varargout{1} = DATA;
    return;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'LoadAuto',8)
        DATA = LoadAll(DATA,'auto','force');
        set(DATA.lstui,'string',DATA.strings);
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
    elseif strncmpi(varargin{j},'Load',4)
        if ~isfield(DATA,'name') && ischar(name)
            fprintf('%s is not a directory\n',name);
            return;
        elseif isdir(DATA.name)
          if strncmpi(varargin{j},'LoadSpikes',8)
            DATA = LoadAll(DATA,[],'loadspikes');
          else
            DATA = LoadAll(DATA,[],'force');
          end
            set(DATA.lstui,'string',DATA.strings);
        else
            fprintf('%s is not a directory\n',DATA.name);
            return;
        end
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
    end
    j = j+1;
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
    CheckClusters(DATA,'exclusions');
    CheckClusters(DATA,'nclusters');
    CheckClusters(DATA,'fitspace');
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
    DATA = LoadExtra(DATA);
elseif strncmpi(name,'reloadcellist',9)
    DATA = LoadCellFile(DATA);
elseif strncmpi(name,'readres',7)
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
elseif strncmpi(name,'CallAllVPcs',7)
    vname = FileName(DATA,DATA.currentpoint(1),'FullV');
    if exist(vname,'file')
        AllVPcs(vname,'thcan',DATA.currentpoint(2),'nprobepc',1,'tryall','reapply');
    end

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
    GetFigure(DATA.tag.spikes,'front');
    PlotSyncSpikes(DATA, eid, fliplr(probes),[1 1], varargin{:});
    return;
    DATA.plotcells.showmahal = 2;
    DATA.plotcells.showfitdp = 2;
    PlotCellList(DATA, 'showfig');
elseif strncmpi(name,'plotxcorrprobes',12)
    SetFigure(DATA, DATA.tag.xcorr);
    PlotXcorrs(DATA, DATA.xcorrs, 0);
elseif strncmpi(name,'plotxcorr',6)
    SetFigure(DATA, DATA.tag.xcorr);
    PlotXcorrs(DATA, DATA.xcorrs, 1);
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
if initcall && isfield(DATA,'exptid') %Do pop plot if data loaded
    DATA = PlotAllClusters(DATA,[]);
elseif initcall
    DATA = LoadCellFile(DATA);
    DATA = LoadExtra(DATA);
    PlotCellList(DATA,'showfig');
end
set(DATA.toplevel,'UserData',DATA)

function out = doentry(DATA, Clusters, id)

type = DATA.plot.onetype;

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
     GetFigure('Xcorrs')
     imagesc(squeeze(cres.synci(:,:,1)));
     DATA.synic{id} = cres;
     xc = cres;
     ofile = [DATA.name '/' strrep(DATA.strings{id},'ClusterTimes','Xcorrs')];
     save(ofile,'xc');
end
out = DATA;

function DATA = SetDefaults(DATA)
    
    DATA.version = 1.2;
    DATA.colors = mycolors('spkcolors')

    DATA.clusteroffset = 0;
    DATA.currenttrial = 1;
    DATA.spikesloaded = 0;
DATA.usesavedcodes = 0;
DATA.Comments = [];
DATA.cellcluster = 1;  %current cluster for defining cells
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
DATA.tag.allexpt = [DATA.tag.top 'AllExpt'];
DATA.tag.xyplot = [DATA.tag.top 'XYplot'];
DATA.tag.xcorr = [DATA.tag.top 'xCorr'];
DATA.tag.hist = [DATA.tag.top 'Histogram'];
DATA.tag.misc = [DATA.tag.top 'Misc'];
DATA.tag.spikes = [DATA.tag.top 'Spikes'];
DATA.tag.allspikes = [DATA.tag.top 'AllSpikes'];
DATA.tag.comments = [DATA.tag.top 'Comments'];
DATA.tag.xyseq = [DATA.tag.top 'XYseq'];
DATA.rebuild = 0;
DATA.cellbackup = 0;
DATA.spoolspikes = 0;
DATA.showspkxy=0;
DATA.showspkmean=1;
DATA.plotmeantype='both';
DATA.refitgm = 0;
DATA.plothist = 0;
DATA.plotexpt = 0;
DATA.plotallxy = 0;
DATA.plotxyseq = 0;
DATA.plottrighist = 0;
DATA.plotspks = 0;
DATA.plotspk.bytrial = 1;
DATA.plotspk.showexcluded = 0;
DATA.plotspk.showcellmeans = 0;
DATA.plotcells.showmahal = 0;
DATA.plotcells.showfitdp = 0;
DATA.plot.gmcid = 0;
DATA.plot.xcmax = 200;
DATA.plot.synctmax = 2;

DATA.id = 1;
DATA.exclude.offpeak = 0;
DATA.exclude.onpeak = 0;
DATA.exclude.noncell = 0;
DATA.colorscheme = 'Plain';
DATA.datatype = 1;
DATA.templatesrc = 1;
DATA.show.exptno = 0;
DATA.show.exptname = 0;
DATA.show.ed = 0;
DATA.SpaceTypeLabels = {'PCs',  'ADC', 'Template', 'Var-E', 'undef', 'N-D'};
DATA.mahaltype = 1;
DATA.mahalcrit = 2;
DATA.xccrit = 0.9;
DATA.CellList = [];
DATA.CellDetails = [];
DATA.CellChanges = [];
DATA.steptype = 1;
DATA.nprobespool = 1;
DATA.currentpoint = [1 1];
DATA.NewCut.exptid = 0;
DATA.NewCut.probe = 0;
DATA.plotexpttype = 'means';

DATA.plot.onetype = 'dips';
DATA.plot.density = 0;
DATA.plot.mahalcmax = 5;
DATA.plot.alltype = 'mahal+var';
DATA.plot.showboundary = 0;
DATA.markcell.candidates = 1;
DATA.markcell.dropi = 0; %show cells with low drop index
DATA.markcell.ellipses = 0;
DATA.markcell.mahal = 0; %show cells with low mahal distance


function TagMenu(a, b, fcn)
    DATA = GetDataFromFig(a);
    reason = strmatch(fcn,{'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear'});
    if isempty(reason)
        reason = NaN;
    end
    if reason == 10
        DATA.selectprobe = zeros(size(DATA.selectprobe));
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
    SaveComments(DATA);
    set(DATA.toplevel,'UserData',DATA);

function GetComment(DATA,a,b)
    f = SetFigure(DATA,DATA.tag.comments);
    str = sprintf('Ex%.1f (%d), Probe%d',DATA.exptid(a),a,b);
    it = findobj(f,'tag','ProbeLabel');
    if length(it) == 1
        set(it,'string',str);
    end
  
function AddComment(a,b,str)
    DATA = GetDataFromFig(a);
    n = length(DATA.Comments) + 1;
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
    SaveComments(DATA);
    set(DATA.toplevel,'UserData',DATA);
    
function SaveComments(DATA)
    outname = [DATA.name '/Comments.mat'];
    Comments = DATA.Comments;
    Tagged = DATA.tagged;
    save(outname,'Tagged','Comments');

function PrintComments(DATA,e,p)
    if isempty(DATA.Comments)
        return;
    end
    
    for j = 1:length(p)
        id = find([DATA.Comments.ex] == e & [DATA.Comments.p] == p(j));
        for k = 1:length(id)
            fprintf('E%.0f P%d %s\n',DATA.Comments(id(k)).exptno,p(j),DATA.Comments(id(k)).comment);
        end
    end

            

    
function DATA = LoadComments(DATA)
    outname = [DATA.name '/Comments.mat'];
    if exist(outname,'file')
    load(outname);
    if exist('Comments','var')
        DATA.Comments = Comments;
    end
    if exist('Tagged','var')
        DATA.tagged = Tagged;
    end
    end
    
function PlotMenu(a, b, fcn, type, varargin)
    DATA = GetDataFromFig(a);
    onoff = {'off' 'on'};
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if strcmp(fcn,'expt')
        PlotAllCell(DATA, type);
        DATA.plotexpttype = type;
        m = get(a,'Parent');
        c = get(m,'Children');
        for j = 1:length(c)
            set(c(j),'Checked','off');
        end
        set(a,'checked','on');
        PlotExptCounts(DATA,DATA.currentpoint(1),DATA.currentpoint(2));
    elseif strcmp(fcn,'spikes')
        if strmatch(type,'Spool')
        end
    elseif strcmp(fcn,'xcorr')
        if strmatch(type,'zoom')
            DATA.plot.xcmax = varargin{1};
        elseif strmatch(type,'spkzoom')
            DATA.plot.synctmax = varargin{1};
        end
    elseif strcmp(fcn,'probes')
       if strmatch(type,{'spoolall' 'spoolcells' 'allspks' 'allcellspks' 'spooleverything' 'spooleverycell'})
            DATA = PlotAllProbe(DATA, type);
       elseif strmatch(type,{'AllMean' 'AllMeanIm' 'AllExptIm'})
           PlotAllProbeMean(DATA, type);
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
           AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
           eid = DATA.xcorrs(DATA.xcid).eid;
           probes = find(DATA.selectprobe(eid,:));
           GetFigure(DATA.tag.spikes);
           SpoolAllProbes(DATA, eid, AllSpikes(eid,DATA.xcorrs(DATA.xcid).cells), Clusters{eid});

       elseif strcmp(type,'AllXY')
           PlotAllProbeXY(DATA);
       else
            PlotAllCell(DATA, type);
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
        elseif strcmpi(type,'allxy')
            PlotAllCellXY(DATA);
        elseif strcmpi(type,'markdropi')
            if DATA.markcell.dropi > 0
                DATA.markcell.dropi = 0;
                set(a,'Checked','off');
            else
                DATA.markcell.dropi = 1.5;
                set(a,'Checked','on');
            end
            PlotCellList(DATA);
        elseif strcmpi(type,'markmahal')
            if DATA.markcell.mahal > 0
                DATA.markcell.mahal = 0;
                set(a,'Checked','off');
            else
                DATA.markcell.mahal = 2.5;
                set(a,'Checked','on');
            end
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
            PlotAllCellMean(DATA,'lineonly');
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
        end
    elseif strncmp(fcn,'means',4)
        DATA.plotmeantype = type;
        C = getappdata(DATA.toplevel,'Clusters');
        GetFigure(DATA.tag.spkmean);
        set(a,'Checked','on');
        PlotMeanSpike(C{e}{p}, p,0,type);
    end
    
    set(DATA.toplevel,'UserData',DATA);
    
    
    function E = CountSpikes(Expt, C, clnum,xcl)
        
        latency = 500;
        id = find(C.clst == clnum);
        t = (C.times(id) .*10000)-latency;
        for j = 1:length(Expt.Trials)
            T = Expt.Trials(j);
            Expt.Trials(j).count = sum(t > T.Start(1) & t <= T.End(end));
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
            Expts{id(j)} = CountSpikes(Expts{id(j)},Clusters{e(j)}{p(j)},c(j)+1,xid);
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
        
    function [true, cellid] = isacell(DATA, row, p)
        cellid = 0;
        
        true = sum(DATA.CellList(row,p,:),3) > 0;
        if true
            cellid = DATA.CellList(row,p,:);
            cellid = cellid(cellid>0);
        end
    
        function ApplyLayout(DATA)
    
    if ~exist(DATA.layoutfile)
        fprintf('Cant read %s\n',DATA.layoutfile);
        return;
    end
    load(DATA.layoutfile);
    setappdata(DATA.toplevel,'Figpos',Figpos);
    f = fields(Figpos);
    for j = 1:length(f)
        it = findobj('type','figure','tag',f{j});
        if length(it) == 1
                set(it,'Position',Figpos.(f{j}));
        end
    end


    function OptionMenu(a, b, fcn)
        
        DATA = GetDataFromFig(a);
        onoff = {'off' 'on'};
        if strcmp(fcn,'usesavedcodes')
            DATA.usesavedcodes = ~DATA.usesavedcodes;
            set(a,'Checked', onoff{1+DATA.usesavedcodes});
        elseif strcmp(fcn,'tofront')
            FiguresToFront(DATA.tag);
        elseif strcmp(fcn,'loadlayout')
            if ~isfield(DATA,'layoutfile') || isempty(DATA.layoutfile)
                DATA.layoutfile =  '/bgc/group/matlab/preferences/PlotClusters/Bruce.layout.mat';
            end
            [afile, pathname] = uigetfile(DATA.layoutfile);
            DATA.layoutfile = [pathname afile];
            ApplyLayout(DATA);
        elseif strcmp(fcn,'savelayout')
            f = fields(DATA.tag);
            for j = 1:length(f);
                it = findobj('type','figure','Tag', DATA.tag.(f{j}));
                if length(it) == 1
                    Figpos.(DATA.tag.(f{j})) = get(it,'Position');
                end
            end
            [outname, pathname] = uiputfile(DATA.layoutfile);
            if outname
                DATA.layoutfile = [pathname outname];
                save(DATA.layoutfile,'Figpos');
                fprintf('Layout saved to %s\n',DATA.layoutfile);
            end
        elseif fcn == 1
            DATA.plot.density = ~DATA.plot.density;
            set(a,'Checked',onoff{DATA.plot.density+1});
            if isfield(DATA,'currentpoint')
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
        elseif strcmp(fcn,'excludecl2')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),2);
        elseif strcmp(fcn,'usealltrialscl1')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),1, 'reset');
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
        elseif strcmp(fcn,'xcorrcells')
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
                    id = find(cellid(:,1) ==j && cellid(:,2) == k || cellid(:,2) ==j && cellid(:,1) == k);
                    if length(id) == 1 && DATA.xcorrs(id).calctime > P.savetime(1)
                        xc = DATA.xcorrs(id).xc;
                    else
                        pid = find(P.clst == cnum(j));
                        qid = find(Q.clst == cnum(j));
                        xc = xcorrtimes(P.times(pid),Q.times(qid));
                    end
                    if k == j
                        xc(201) = 0;
                    end
                    subplot(np,np,k+(j-1)*np);
                    bar(-200:200,xc);
                    axis('tight');
                    if k == j
                        title(sprintf('P%d',cid(k)));
                    end
                    if k == 1
                        ylabel(sprintf('P%d',cid(j)));
                    end
                end
            end
            toc
        elseif strcmp(fcn,'xcorrallcells') ||  strcmp(fcn, 'xcorrallprobes') || strcmp(fcn, 'recalcxcorrcell')

            if  strcmp(fcn, 'recalcxcorrcell')
                recalc = 1;
                DATA = rmfield(DATA,'xcorrs');
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
            if isfield(DATA,'xcorrs') && recalc == 0
                probes = cat(1,DATA.xcorrs.probes);
                cellid = cat(1,DATA.xcorrs.cells);
            else
                probes = [0 0];
                cellid = [0 0];
            end
            np = length(cid);
            tic;
            nxc = 1;
            DATA.xcorrval.times = [-0.2:0.001:0.2];
            for e = 1:length(DATA.exptid)
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
                        if recalc == 0 && length(id) == 1 && DATA.xcorrs(id).calctime > max([P.savetime(1) Q.savetime(1)])
                            xc = DATA.xcorrs(id).xc;
                        else
                            if length(id) >= 1
                                nxc = id(1);
                            else
                                nxc = length(DATA.xcorrs)+1;
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
                            DATA.xcorrs(nxc).calctime = now;
                        end
                    end
                end
            end
            toc
            SaveExtras(DATA);
            set(DATA.toplevel,'UserData',DATA);
            PlotXcorrs(DATA, DATA.xcorrs, bycell);
        elseif strcmp(fcn,'usealltrialscl2')
            DATA = ExcludeTrialsForCell(DATA, DATA.currentpoint(2),2, 'reset');
        elseif fcn == 8
            DATA.plot.gmcid = ~DATA.plot.gmcid;
            set(a,'Checked',onoff{DATA.plot.gmcid+1});
        elseif fcn == 9
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
                        pid = find(P.clst == 1+j);
                        Q = Clusters{ps(k)}.next;
                        qid = find(Q.clst == 2);
                        xc = xcorrtimes(P.times(pid),Q.times(qid));
                    end
                    if k == j
                        xc(201) = 0;
                    end
                    subplot(np,np,k+(j-1)*np);
                    bar(-200:200,xc);
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
            GetFigure(DATA.tag.allxy,'front');
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
            CellFinder(DATA);
        elseif strmatch(fcn,{'showexcluded' 'showcellmeans'})
            DATA.plotspk.(fcn) = ~DATA.plotspk.(fcn);
            set(a,'Checked', onoff{1+DATA.plotspk.(fcn)});
        end
        set(DATA.toplevel,'UserData',DATA);
        
function PlotExpts(DATA)
    GetFigure('Means');
    Expts = Expts{DATA.currentpoint(1)};
   ex =      PlotExpt(Expts{1});
   ts = now;
    for j = 1:length(Expts)
        PlotRates(Expts{j},ex);
%        PlotExpt(Expts{j},'hold');
        hold on;
    end
     mytoc(ts);
     
     
function PlotXcorrs(DATA, xcorrs, bycell)
    ClearPlot;
    if bycell
    cellids = cat(1,xcorrs.cells);
    cellids(isnan(cellids)) = 0;
    else
        cellids = cat(1,xcorrs.probes);
    end
    weights = prod(cat(1,xcorrs.n)');
    cells = unique(cellids);
    cells = cells(cells > 0);
    probes = cat(1,xcorrs.probes);
    np = length(cells);
    for j = 1:length(cells)
        ida = find(cellids(:,1) == cells(j));
        idb = find(cellids(:,2) == cells(j));
        cellpos(j) = (sum(probes(ida,1))+sum(probes(idb,2)))./(length(ida)+length(idb));
        for k = 1:j
            ida = find(cellids(:,1) == cells(j) & cellids(:,2) == cells(k));
            idb = find(cellids(:,2) == cells(j) & cellids(:,1) == cells(k));
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
            if length(id)
                if length(id) > 1
                    xc = WeightedSum(cat(1,xcorrs(id).xc),weights(id));
                else
                    xc = xcorrs(id).xc;
                end
                mysubplot(np,np,k+(j-1)*np);
                h = plot(-200:200,xc,'k-','linewidth',2);
                axis('tight');
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k)]});
                set(h,'buttondownfcn',{@HitXcorrAll,bycell,[cells(j) cells(k)]});
                if k == j
                    [a,b] = find(DATA.CellList == cells(k));
                    p = 1+mod(b-1,DATA.nprobes);
                    if j == 1
                        title(sprintf('Cell%d',cells(k)));
                    else
                        ii = find(icells == cells(k-1));
                        ij  = find(icells == cells(j));
                        if bycell
                            h = text(xl(1),yl(2),sprintf('C%d at %.1f (%.1f)',cells(k),mean(p),separation(ii,ij)));
                            set(h,'HorizontalAlignment','left','verticalalignment','bottom');
                        else
                            title(sprintf('P%d',k));
                        end
                    end
                end
                if k == 1
                    ylabel(sprintf('Cell%d',cells(j)));
                end
            end
        end
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
DATA.plotexpt = GetCheck('PlotExpt',DATA.toplevel);
DATA.refitgm = GetCheck('RefitGM',DATA.toplevel);
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

function id = FindSpikes(DATA,C, xcl)
    xid = [];
    e = DATA.currentpoint(1);
    p = C.probe(1);
    t = C.times.*10000;
    Expts = getappdata(DATA.toplevel,'Expts');
    e = floor(C.exptno); %exlusion list is by Exp # not row
    e = C.exptid; %but Expts list is as long as rows
    if max(xcl) > length(Expts{e, p}.Trials)
        fprintf('Error In Trial Exclusion list\n');
    else
        for j = 1:length(xcl)
            id = find(t > Expts{e, p}.Trials(xcl(j)).Start(1) & t< Expts{e, p}.Trials(xcl(j)).End(end));
            xid = [xid id];
        end
    end
    id = setdiff(1:length(C.times),xid);

    function plots = PlotClusterPoints(C, uid, varargin)
        plotgmcid = 0;
        clnum = 1;
        j = 1;
        colors = mycolors('spkcolors');
        while j <= length(varargin)
            if strncmpi(varargin{j},'colors',6)
                j = j+1;
                colors = varargin{j};
            elseif strncmpi(varargin{j},'plotgmcid',8)
                plotgmcid = 1;
            end
            j = j+1;
        end
        if isempty(uid)
            uid = 1:size(C.xy,1);
        end
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
            id = find(C.clst(uid) == j);
            plots(j) = plot(C.xy(uid(id),1),C.xy(uid(id),2),'.','markersize',1,'color',colors{j});
        end

function plots = PlotClusterXY(DATA, C, varargin)
    titlemode = 0;
    clnum = 1;
    xyargs = {};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plotgmcid',8)
            xyargs = {xyargs{:} varargin{j}};
        elseif strncmpi(varargin{j},'shorttitle',8)
            titlemode = 1;
        else
            xyargs = {xyargs{:} varargin{j}};
        end
        j = j+1;
    end
    lincolor = 'r';
    
    p = C.probe(1);
    e = C.exptid;

    xcl = FindExcludedTrials(DATA, e, p,1,C);
    if length(xcl)
        uid = FindSpikes(DATA, C, xcl);
    else
        uid = 1:size(C.xy,1);
    end
    if DATA.plot.density
        lincolor = 'w';
        plots = DensityPlot(C.xy(uid,1),C.xy(uid,2),'ynormal');
        hold on;
    else
        plots = PlotClusterPoints(C,uid,xyargs{:});
    end
    if isfield(C,'space')
        spstr = [DATA.SpaceTypeLabels{C.space(1)} sprintf(' %d',C.space(2:end))];
    else
        spstr = '';
    end
    if isfield(DATA,'GaussFitdp') && size(DATA.GaussFitdp,1) >= e
        spstr = [sprintf('(Fit %.1f)',DATA.GaussFitdp(e,p,2)) spstr];
    end
    if isfield(C,'gmfit2d') && DATA.plot.showboundary
            xy = GMBoundary(C); %,'plot to show image
            plot(xy(:,1),xy(:,2),'k');
    end
    axis('tight');
    if isfield(C,'bestspace') && C.space(1) == 6
        if titlemode == 1
           h= title(sprintf('E%.0fP%d %.2f  %s %.2f',C.exptno,p,C.mahal(1),spstr,C.mahal(4)));
        else
           h= title(sprintf('P%d Ex %.0f Gm %.2f  %s (%.2f,%.2f for space %.0f)%.0f',p,C.exptno,C.mahal(1),spstr,C.mahal(4),C.bestspace(1),C.bestspace(2),C.sign));
        end
    else
        if titlemode == 1
            h=title(sprintf('E%.0fP%d %.2f  %s %.2f',C.exptno,p,C.mahal(1),spstr,C.mahal(4)));
        else
            h =title(sprintf('P%d Ex %.0f Gm %.2f  %s (%.2f)%.0f',p,C.exptno,C.mahal(1),spstr,C.mahal(4),C.sign));
        end
    end
    if C.shape == 1
        line([C.crit(1) C.crit(1)],get(gca,'ylim'),'color',lincolor);
    elseif C.shape == 2
        line([C.crit(1) C.crit(1)],get(gca,'ylim'),'color',lincolor);
    elseif C.shape == 0
        C.down = 0;
        DrawEllipse(C,'color',lincolor);
    end
    for j = 1:length(C.next) 
        if ~isempty(C.next{j})
        C.next{j}.down = 0;
        if C.next{j}.shape == 0
            DrawEllipse(C.next{j},'color',DATA.colors{j+2});
        end
        end
    end
    if isacell(DATA,e,p)
        sz = get(h,'fontsize');
        sz = sz .* 1.4;
        set(h,'color','b','fontweight','bold','fontsize',sz);
    elseif C.auto == 1
        set(h,'color','r');
    elseif DATA.plot.density
        set(h,'color','w');
    end
    
    function r = CalcRadius(E)
                
        rx = E.xyr(3);
        ry = E.xyr(4);
        xys = xyrotate(E.xy(:,1)-E.xyr(1),E.xy(:,2)-E.xyr(2),E.angle);
        r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;

    function [x, nsp] = PlotHist(xy, varargin)
    E = [];
    j = 1;
    while j <= length(varargin)
        if isfield(varargin{j},'pos') || isfield(varargin{j},'crit')
            E = varargin{j};
        end
        j = j+1;
    end
    
    hold off;

    if E.shape == 0
        r = CalcRadius(E);
        [nsp,x] = hist(r,200);
        bar(x,nsp,1);
        scale = trapz(x,nsp);
    else
    [nsp,x] = hist(xy(:,1),200);
    bar(x,nsp,1);
    scale = trapz(x,nsp);
    end

        
    if E.shape == 0
        hold on;
        plot([1 1],[0 max(nsp)],'r');
    elseif isfield(E,'gmfit1d')
        [a,b] = GMDip(xy,E.gmfit1d);
        hold on;
        plot(b.gxy(:,1),b.gxy(:,2).*scale,'color','g');
        plot(b.gxy(:,1),b.gxy(:,3).*scale,'color','g');
        plot(b.gxy(:,1),sum(b.gxy(:,[2 3]),2).*scale,'r');
        plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
        plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
    else
        b.gxy(:,1) = linspace(min(E.xy(:,1)),max(E.xy(:,1)));
        gmsd = NaN;
    plot([E.crit(1) E.crit(1)],[0 max(nsp)],'r');
    end

 function PlotExptCounts(DATA, e, p, varargin)
       tag = DATA.tag.expt;

       j = 1;
       while j <= length(varargin)
           if strcmpi(varargin{j},'tag')
               j = j+1;
               tag = varargin{j};
           end
           j = j+1;
       end
       
     GetFigure(tag);
     Clusters = getappdata(DATA.toplevel,'Clusters');
     Expts = getappdata(DATA.toplevel,'Expts');
     clnum = 2;
     Expt = Expts{e,p};
     spkt = Clusters{e}{p}.times(Clusters{e}{p}.clst == clnum) .*10000;

     for j = 1:length(Expt.Trials)
         id = find(spkt > Expt.Trials(j).Start(1) & spkt < Expt.Trials(j).End(end));
         Expt.Trials(j).Spikes = round(spkt(id)'-Expt.Trials(j).Start(1));
     end
     if strncmpi(DATA.plotexpttype,'trialcounts',10)
         PlotExpt(Expt,'seqt');
     else
         PlotExpt(Expt,'shown');
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
      
Clusters = getappdata(DATA.toplevel,'Clusters');
Expts = getappdata(DATA.toplevel,'Expts');
    oldf = gcf;
Expt = [];
if size(Expts,1) >= ex
    Expt = Expts{ex,p};
    DATA.Expt  = Expt;
end
if DATA.datatype == 2
C = Clusters{ex}{p}.cluster{DATA.templatesrc};
else
C = Clusters{ex}{p};
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
end
DATA.markh = MarkCurrentCluster(DATA);
PrintComments(DATA,ex,p);
fprintf('Mahal ND %.2f(%d), 2D %.2f, 1D %.2f. Dropi %.2f (T%.2f). Made %s\n',C.bestspace(1),bestspace,C.mahal(1),C.mahal(4),C.dropi(3),C.Trigger,datestr(C.ctime))
if isfield(DATA,'GaussFitdp')
    fprintf('Fit %.1f %.1f\n',DATA.GaussFitdp(ex,p,1),DATA.GaussFitdp(ex,p,2));
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
    QuickSpikes(DATA, [ex p]);
    drawnow;
end
if DATA.showspkxy
    SetFigure(DATA,DATA.tag.xyplot,'front');
    hold off;
    PlotClusterXY(DATA,C);
    DATA.NewCut.exptid = 0;
    DATA.NewCut.probe = 0;
    if DATA.plottrighist
        AddTrigHist(C);
    end
%    drawnow;
end
if ~strcmp(DATA.plotexpttype,'none')
    PlotExptCounts(DATA, ex,p);
end
if DATA.plothist || DATA.refitgm
    GetFigure(DATA.tag.hist,'front');
    [x, nsp] = PlotHist(C.xy,C);
    hold on;
    scale = trapz(x,nsp);
    [dp, fits,c] = Fit2Gauss(C,200);
    if sum(c.fitpos) == 2
        qstr = [];
    else
        qstr = '*';
    end
    if C.shape == 0
        b.gxy(:,1) = linspace(min(x),max(x),100);
    elseif isfield(C,'gmfit1d') 
        [a,b] = GMDip(C.xy,C.gmfit1d);
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
    plot(b.gxy(:,1),fya,'m-','linewidth',2);
    plot(b.gxy(:,1),fyb,'c-','linewidth',2);
    end
    title(sprintf('%d/%d events M%.2f From Fit %.2f%s',C.ncut,size(C.xy,1),C.mahal(4),dp,qstr));
end
if DATA.refitgm
    if C.shape == 0
        C.r = CalcRadius(C);
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
    GetFigure(DATA.tag.spkmean,'front');
    PlotMeanSpike(C,p,0,'addtitle',str,DATA.plotmeantype);
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
    DATA.markcc = DrawBox(ex, p);
end
DATA.currentpoint = [ex, p];
set(DATA.toplevel,'UserData',DATA);
PlotCellList(DATA,'showfig');
drawnow;
figure(oldf);

function AddTrigHist(C, varargin)
    %superimposes a histogram of trigger values on current axes.
if ~isfield(C,'vhist');
    return;
end
    h = ishold;
hold on;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
x = linspace(xl(1),xl(2),length(C.vhist));
plot(x,yl(1)+C.vhist.*diff(yl)./max(C.vhist));
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
    GetFigure('Correlogram');
    plot(xc);

function h = PlotMeanSpike(C, p, cluster, varargin)
    addstr = [];
    plots = [ 1 1];
    colors = mycolors('spkcolors');
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'addtitile',5)
            j = j+1;
            addstr = varargin{j};
        elseif strncmpi(varargin{j},'lineonly',5)
            plots = [0 1];
        elseif strncmpi(varargin{j},'imageonly',7)
            plots = [1 0];
        elseif strncmpi(varargin{j},'exptim',6)
            plots = [2 0];
        elseif strncmpi(varargin{j},'meanlines',6)
            plots = [3 0];
        end
        j = j+1;
    end

    chspk = C.probe(1) + [-1:1];
    voff = [-1:1] .*2;
    
    if cluster == 0 %plot all
        nclusters = 1;
        for j = 1:length(C.next)
            if ~isempty(C.next{j})
                nclusters = nclusters+1;
            end
        end
        if plots(1) == 3 
        PlotMeanSpike(C,p,1,'meanlines');
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
        if plots(2)
        PlotMeanSpike(C, p, 1, 'lineonly');
        end
        nc = 1;
        for j = 1:length(C.next)
            if ~isempty(C.next{j})
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
        C.next{cluster-1}.exptno = C.exptno;
        C = C.next{cluster-1};
    end
    if plots(1) == 1
        hold off;
            h(1) = imagesc(C.MeanSpike.ms);
        if p <= 0 && isfield(C,'probe');
            p = C.probe(1);
        end
        line([0 5],[p p],'color','r');
        title(sprintf('P%d/%d Ex %.1f Gm %.2f (%.2f) %s',p,cluster,C.exptno,C.mahal(1),C.mahal(2),addstr));
        if sum(plots) > 1
            subplot(1,2,2);
        end
    elseif plots(1) == 3 %mean lines
        subplot(1,1,1);
        for j = 1:length(chspk)
             plot(C.MeanSpike.ms(chspk(j),:)+voff(j),'r');
             hold on;
             for k = 1:length(C.next)
                 if isfield(C.next{k},'MeanSpike')
                     plot(C.next{k}.MeanSpike.ms(chspk(j),:)+voff(j),'color',colors{k+2});
                 end
             end
        end
    end
    if plots(2)
        hold off;
        v = std(C.MeanSpike.ms');
        id = find(v > max(v)/2);

        for j = id
            h(2) = plot(C.MeanSpike.ms(j,:),'r');
            hold on;
            if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j
                plot(C.MeanSpike.dp(j,:),'g');
            end
            plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);
        end
    end

function HitImage(a,b, type)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');
ax = get(a,'Parent');
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
ex = round(xy(1,2));
p = round(xy(1,1));
bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})
zval = NaN;
for j = 1:length(l)
    a = get(l(j));
    if isfield(a,'CData')
        Z = get(l(j),'CData');
        zval = Z(ex,p);
    end
end
fprintf('Hit %.0f,%.0f %.3f\n',ex,p,zval);

offset = DATA.clusteroffset;
if ismember(type, [1 3]) % 3 = hit cell image - set cell#
    if DATA.datatype == 2
        C = Clusters{ex}{p}.cluster{DATA.templatesrc};
    else
        C = Clusters{ex}{p};
    end
    fprintf('P%d E%.0f cut on %s\n',p,C.exptno,datestr(C.ctime));
    if type == 3
       it = findobj('Tag','CellNumberId');
       if DATA.CellList(ex,p,offset+1) > 0
           DATA.currentcell = DATA.CellList(ex,p,offset+1);
           set(it,'value',DATA.currentcell);
       end
    end
    if DATA.plotspks
        PlotAllCell(DATA,'allspks');
    end
    if DATA.plotallxy && type == 3  %hit cell image plot
        PlotAllCellXY(DATA);
    end
elseif type == 2
    id = find(DATA.AllPairs(:,1) == ex & DATA.AllPairs(:,2) == p);
    if length(id)
    xc = meanccf(DATA, id, ex, p);
    dist(1) = meanmahal(DATA, id, ex);
    dist(2) = meanmahal(DATA, id, p);
    GetFigure('OneCell');
    plot(xc);
    title(sprintf('Expts %s mahal P%d %.2f, P%d %.2f',num2str(DATA.allexpt(id)),ex,dist(1),p,dist(2)));
    end
    set(gca,'ylim',[0 max(xc)]);
    GetFigure(DATA.tag.spkmean);
    subplot(1,2,1);
    ms = MeanSpike(DATA, DATA.allexpt(id), p);
    imagesc(ms);
    subplot(1,2,2);
    ms = MeanSpike(DATA, DATA.allexpt(id), ex);
    imagesc(ms);
    C = Clusters{DATA.allexpt(id(1))}{p};
    ex= DATA.allexpt(id(1));
end
DATA.currentpoint(1) = ex;
DATA.currentpoint(2) = p;
ShowData(DATA,ex,p,'oneprobe');



function HitXcorrAll(a,b, type, cells)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');
if type == 2
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
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
    text(xl(2),yl(2),sprintf('Cell %d->%d P%d->%d %.0fms',xc.cells(1),xc.cells(2),cells(2),cells(3),xpts(b).*1000),'horizontalalignment','right','verticalalignment','top');
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
 icells = unique([DATA.xcorrs.cells]);
 icells = icells(icells > 0);
 np = length(icells);
 pa = mean(cat(1,probes(ida,1), probes(idb,2))); 
 pb = mean(cat(1,probes(idb,1), probes(ida,2))); 
end
 mysubplot(2,2,2);
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
    mysubplot(np, np, k);
    h = plot(DATA.xcorrval.times,X.xc,'r-');
    set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@HitXcorrAll, 2, [X.eid X.probes X.clnum ids(j)]});
    set(h,'buttondownfcn',{@HitXcorrAll, 2, [X.eid X.probes X.clnum ids(j)]});
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
    mysubplot(np, np, k);
    delete(gca);
    end
end
 
 function HitXcorr(a,b, id, ex, cells)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');

if strcmp(id,'zoom')
    xl = get(gca,'xlim');
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})
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
GetFigure('OneCell');
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
GetFigure(DATA.tag.onecluster);
hold off; 
plot(C{pt}.xy(:,1),C{pt}.xy(:,2),'.', 'markersize', sz);
hold on;
yl = get(gca,'ylim');
plot([C{pt}.crit(1) C{pt}.crit(1)],yl,'r-');
[a,b] = smhist(C{pt}.xy(:,1),100);
plot(b,yl(1)+ a.*range(yl)./max(a),'r');
title(sprintf('Probe %d',pt));


function HitPopPoint(a,b, ex, p)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');

C = Clusters{ex}{p};
fprintf('Expt %d Probe %d\n',C.exptno,p);
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
 Clusters = FixClusters(Clusters);
     

function [Clusters, details] = LoadClusters(name)

exptno = -1;

load(name);
Clusters = Clusters;
for j = 1:length(Clusters)
    if isfield(Clusters{j},'exptno')
        exptno = Clusters{j}.exptno;
    end
    if isfield(Clusters{j},'clst')
        fprintf('%s P%d has clst\n',name,j);
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
load(dname);
if exist('FullVData','var')
    details.FullVData = FullVData;
end
for j = 1:length(ClusterDetails)
    if isfield(ClusterDetails{j},'xy')
    Clusters{j}.xy = ClusterDetails{j}.xy;
    if isfield(ClusterDetails{j},'clst')
        Clusters{j}.clst = ClusterDetails{j}.clst;
    else
        Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
        id = find(ismember(ClusterDetails{j}.t,Clusters{j}.times));
        Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
        Clusters{j}.clst(id) = 2;
    end
        Clusters{j}.times = ClusterDetails{j}.t;
    if exptno >= 0
    Clusters{j}.exptno = exptno;
    end
    end
    if isfield(ClusterDetails{j},'Evec')
        Clusters{j}.Evec = ClusterDetails{j}.Evec;
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
        if ~isfield(CellDetails,'exptids')
            for j = 1:size(Expts,1)      
                DATA.CellDetails.exptids(j) = Expts{j,1}.Header.exptno;
            end
        end
        if ~isfield(CellDetails,'trialids')
            for j = 1:size(Expts,1)
                DATA.CellDetails.trialids{j} = [Expts{j,1}.Trials.Trial];
            end
        else
            for j = 1:size(Expts,1)
                Tn = [Expts{j,1}.Trials.id];
                if length(setxor(DATA.CellDetails.trialids{j},Tn)) > 0 && j <= size(DATA.CellDetails.excludetrials,1)
                    fprintf('Expt Trial List Changed for Expt %d\n',Expts{j,1}.Header.exptno);
                    for k = 1:size(DATA.CellDetails.excludetrials,2)
                        if ~isempty(DATA.CellDetails.excludetrials{j,k,1})
                            [a, id] = ismember(Tn,DATA.CellDetails.excludetrials{j,k,1});
                        end
                    end
                end
            end
        end
        ConvertExclusion(DATA);
        if exist('CellListB','var')
            DATA.CellList(1:size(CellListB,1),1:size(CellListB,2),2) = CellListB;
        end
        if exist('CellChanges','var')
            DATA.CellChanges = CellChanges;
        else
            CellChanges = [];
        end
        if DATA.cellbackup == 0
        bakfile = strrep(cellfile,'.mat','back.mat');
        save(bakfile,'CellList','CellDetails','CellChanges');
        [a,b] = NetFilename(bakfile);
        if exist(b,'dir')
            fprintf('Backing up Cell List to %s\n',a);
            save(a,'CellList','CellDetails','CellChanges');
        else
            fprintf('Can''t Backup CellList to Network\n');
        end
        DATA.cellbackup = 1;
        set(DATA.toplevel,'UserData',DATA);
        end
    else
        DATA.CellList = zeros([length(DATA.exptid) DATA.nprobes 2])
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
AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
e = DATA.currentpoint(1);
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
    GetFigure('CheckTimes');
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

function [f, isnew] = SetFigure(DATA, tag, varargin)
        
    [f,isnew] = GetFigure(tag,varargin{:});
    if isnew
        if strfind(tag,'AllXY');
            set(f,'UserData', DATA.toplevel);
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
        elseif strfind(tag,'CellList');
            hm = uimenu(f,'Label','Plots','Tag','PlotMenu');
            uimenu(hm,'Label','AllXY','Callback',{@PlotMenu, 'cells' 'allxy'});
            uimenu(hm,'Label','All Means','Callback',{@PlotMenu, 'cells' 'allmean'});
            uimenu(hm,'Label','All Mean image','Callback',{@PlotMenu, 'cells' 'allmean'});
            uimenu(hm,'Label','All Trigger hist','Callback',{@PlotMenu, 'cells' 'trighist'});
            uimenu(hm,'Label','All Spks','Callback',{@PlotMenu, 'cells' 'allspks'});
            uimenu(hm,'Label','SpkRate','Callback',{@PlotMenu, 'cells' 'spkrate'});
            uimenu(hm,'Label','SpkRate+xy sequence','Callback',{@PlotMenu, 'cells' 'spkrate+xy'});
            uimenu(hm,'Label','All','Callback',{@PlotMenu, 'cells' 'allplots'});
            uimenu(hm,'Label','All - next cell','Callback',{@PlotMenu, 'cells' 'allnext'});
            hm = uimenu(f,'Label','Mark','Tag','MarkMenu');
            uimenu(hm,'Label','Dropi <1.5','Callback',{@PlotMenu, 'cells' 'markdropi'});
            uimenu(hm,'Label','Candidates','Callback',{@PlotMenu, 'cells' 'markcandidates'});
            uimenu(hm,'Label','mahal < 3','Callback',{@PlotMenu, 'cells' 'markmahal'});
            uimenu(hm,'Label','Ellipses','Callback',{@PlotMenu, 'cells' 'markellipse'});
            bp = [0.01 0.001 0.08 0.04];
            uicontrol(gcf,'style','pop','string',num2str([1:30]'), ...
        'Tag','CellNumberId','Callback', {@SetCellNumber, 'change'},...
        'units', 'norm', 'position',bp,'value',1);
            bp = [0.1 0.001 0.08 0.04];
            uicontrol(gcf,'style','pushbutton','string','Set', ...
                'Callback', {@SetCellNumber, 'set'}, 'Tag','SetCell',...
                'units', 'norm', 'position',bp,'value',1);
            bp(1) = bp(1)+bp(3);
            bp(3) = 0.16;
            uicontrol(gcf,'style','pop','string','Just This Square|Cluster 2 This Square|Cluster 3|Cluster 4 This Square|Use Template Peaks|Better than mahal|Spool All Expts|AllXY|AllMean|AllMeanIm|Exclude Selected Trials|Exclude Trials for C2', ...
        'Tag','ClusterModifier','callback',{@CellTrackMenu, 'setmenu'},...
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

    
    set(f,'UserData', DATA.toplevel);
        elseif strfind(tag,'AllSpikes')
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
        uimenu(hm,'Label','Synchronous','Callback',{@PlotMenu, 'probes', 'plotsync'});
        uimenu(hm,'Label','Synchronous Flip','Callback',{@PlotMenu, 'probes', 'plotsyncb'});
             set(f,'UserData',DATA.toplevel);
            
        elseif strfind(tag,'ClustersAll');
            set(f,'UserData', DATA.toplevel);
        elseif strfind(tag,'xCorr')
            set(f,'UserData', DATA.toplevel);
            hm = uimenu(f,'Label','Zoom','Tag','xcZoom');
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

            set(f,'UserData',DATA.toplevel);
            hm = uimenu(f,'Label','Spool','Tag','Spool');
            uimenu(hm,'Label','This Probe','Callback',{@PlotMenu, 'spikes', 'Spool'});
            uimenu(hm,'Label','Sync Spikes','Callback',{@PlotMenu, 'probes', 'spoolsync'});
        elseif strfind(tag,'XYplot');
                hm = uimenu(f,'Label','Cluster','Tag','Cluster');
                uimenu(hm,'Label','SaveChanges','Callback',{@XYCluster, 'save'});
                uimenu(hm,'Label','Cut Ellipses','Callback',{@XYCluster, 'ellipses'});
                uimenu(hm,'Label','Cut Lines','Callback',{@XYCluster, 'lines'});
                uimenu(hm,'Label','Ellipse for Cl2','Callback',{@XYCluster, 'ellipse2'});
                uimenu(hm,'Label','No Cutting','Callback',{@XYCluster, 'nocut'});
                uimenu(hm,'Label','Optimize angle','Callback',{@XYCluster, 'bestangle'});
                uimenu(hm,'Label','Optimize angle (quick)','Callback',{@XYCluster, 'quickangle'});
                uimenu(hm,'Label','replot','Callback',{@XYCluster, 'replot'});
                uimenu(hm,'Label','revert','Callback',{@XYCluster, 'revert'});
                uimenu(hm,'Label','replot with GM ids','Callback',{@XYCluster, 'replotgm'});
                uimenu(hm,'Label','Density','Callback',{@OptionMenu, 1});
                tmp.parentfig = DATA.toplevel;
                set(f,'UserData',DATA.toplevel);
       else
            set(f,'UserData', DATA.toplevel);
        end
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
     GetFigure(DATA.tag.clusters);
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
     GetFigure(DATA.tag.clusters);
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
     GetFigure(DATA.tag.clusters);
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

 function DATA = LoadExtra(DATA)
     Xbits = [];
     outname = [DATA.name '/PlotClusterExtra.mat'];
     if exist(outname,'file')
         load(outname);
         f = fields(Xbits);
         for j = 1:length(f)
             DATA.(f{j}) = Xbits.(f{j});
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
     if ~isempty(Xbits)
         save(outname,'Xbits')
     end

function [dp, fits, details] = Fit2Gauss(C, varargin)

    plottype = 0;
    fx = linspace(min(C.xy(:,1)),max(C.xy(:,1)),200);
    id = find(C.clst ==2);
    nid = find(C.clst ==1);
    if C.shape == 0 %ellipse
        rx = C.xyr(3);
        ry = C.xyr(4);
        xys = xyrotate(C.xy(:,1)-C.xyr(1),C.xy(:,2)-C.xyr(2),C.angle);
        r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
        fx = linspace(min(r),max(r),200);
        [y,x] = hist(r,fx);
        [a,b] = min(abs(fx-1));
        elseif isfield(C,'crit')
        [a,b] = min(abs(fx-C.crit));
    [y,x] = hist(C.xy(:,1),fx);
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
             GetFigure('FitGauss');
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


function XYCluster(src,b, type)
DATA = GetDataFromFig(src);
EClusters = getappdata(DATA.toplevel,'Clusters');
e = DATA.currentpoint(1);
p = DATA.currentpoint(2);
C = EClusters{e}{p};
ex = C.exptno;
    onoff = {'on','off'};
        subtype = strmatch(type,{'bestangle' 'quickangle'});
    if strcmp(type,'save')
       C = DATA.NewCut;
        cfile = [DATA.name '/Expt'  num2str(ex) 'ClusterTimes.mat'];
        load(cfile);
        afile = [DATA.name '/Expt'  num2str(ex) 'AutoClusterTimesDetails.mat'];
        a = load(afile);
        dfile = [DATA.name '/Expt'  num2str(ex) 'ClusterTimesDetails.mat'];
        load(dfile);
        for j = 1:length(a.ClusterDetails)
            if j > length(ClusterDetails) || isempty(ClusterDetails{j}) || ~isfield(ClusterDetails{j},'t')
                ClusterDetails{j} = a.ClusterDetails{j};
            end
        end
        ClusterDetails{p}.clst = C.clst;
        ClusterDetails{p}.xy = C.xy;
        if isempty(Clusters{p})
            Clusters{p} = rmfield(C,{'xy' 'clst' 'Evec'});
        end
        if C.cluster > 1
            n = C.cluster-1;
            Clusters{p}.next{n}.shape = C.next{n}.shape;
            Clusters{p}.next{n}.crit = C.next{n}.crit;
            Clusters{p}.next{n}.xyr = C.next{n}.xyr;
            Clusters{p}.next{n}.sign = 0;
            Clusters{p}.next{n}.angle = C.angle;
            Clusters{p}.next{n}.times = C.times(C.clst == 2);
        else
        Clusters{p}.shape = C.shape;
        Clusters{p}.crit = C.crit;
        Clusters{p}.xyr = C.xyr;
        Clusters{p}.angle = C.angle;
        Clusters{p}.sign = C.sign;
        Clusters{p}.times = C.times(C.clst == 1);
        end
        Clusters{p}.auto = 2;
        Clusters{p}.clusterprog = sprintf('PlotClusters V %.2f',DATA.version);
        testing = 0;
        if testing
            ofile = strrep(cfile,'ClusterTimes','NewClusterTimes');
        else
            ofile = cfile;
        end
        save(ofile,'Clusters');
        if testing
            ofile = strrep(dfile,'ClusterTimes','NewClusterTimes');
        else
            ofile = dfile;
        end
        save(ofile,'ClusterDetails');
        
        fprintf('Save Cluster for probe %d Expt %d\n',p,ex);
    elseif strcmp(type,'nocut') 
        DATA.ellmousept.shape = -1;
    elseif strmatch(type,{'bestangle' 'quickangle'})
        if subtype == 2
            [a,b,c] = BestAngleGM(C.xy, [], [],'quick');
        else
            [a,b,c] = BestAngleGM(C.xy, [], []);
        end
        C.gmfit1d = c.gmfit;
        C.crit = c.crit;
        GetFigure(DATA.tag.xyplot);
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
            GetFigure(DATA.tag.hist)
            PlotHist(c.xy,C);
            title(sprintf('Best Angle %.1f deg mahal %.2f',a *180/pi,b));
            DATA.tag.misc = 'Misc';
            GetFigure(DATA.tag.misc);
            hold off;
            plot(c.angles,c.d);
            if subtype == 1
            hold on;
            plot(c.angles,abs(c.dprime),'r');
            end
        end
    elseif strncmp(type,'ellipse',7) || strcmp(type,'lines')
        F = get(get(src,'parent'),'parent');
        DATA.currentcluster = 1;
        DATA.elmousept.h= -1;
        shapes = [0 1 0];
        DATA.elmousept.shape= shapes(strmatch(type,{'ellipses' 'lines' 'ellipse2'}));
        DATA.elmousept.down = 0;
        DATA.elmousept.done = 0;
        DATA.elmousept.steps = 0;
        DATA.elmousept.angle = 0;
        DATA.elmousept.cluster = 1;
        if strcmp(type,'ellipse2')
            DATA.elmousept.cluster = 2;
            DATA.elmousept.plotargs = {'color' 'g'};
        else
            DATA.elmousept.plotargs = {'color' 'r'};
        end
        DATA.elmousept.color = [1 0 0];
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
        GetFigure(DATA.tag.xyplot); hold off;
        PlotClusterXY(DATA,X{p});
        set(DATA.toplevel,'UserData',DATA);
    end
    
     function in = InGraph(pt, ax)
        xl = get(ax,'Xlim');
        yl = get(ax,'Ylim');
      in = pt(1,1) > xl(1) & pt(1,1) < xl(2) & pt(1,2) > yl(1) & pt(1,2) < yl(2);
    
    function XYButtonPressed(src, data)
DATA = GetDataFromFig(src);

DATA.ts = now;
start = get(gca,'CurrentPoint');
if InGraph(start,gca)
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
%    DATA.elmousept.mode = mode;
    distance = DistanceToEllipse(DATA.elmousept,start(1,1:2));
    DATA.elmousept.steps = 0;
    if ishandle(DATA.elmousept.h)
        delete(DATA.elmousept.h);
    end
    set(gca,'xlimmode','manual','ylimmode','manual');
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
        DATA.elmousept.down = 1;
        DATA.elmousept.done = 0;
        DATA.elmousept.angle = 0;
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
pi = DATA.currentpoint(2);
ei = DATA.currentpoint(1);

if mode == 1
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
elseif mode == 2  %button was pressed inside ellipse, just move it
DATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)]; 
end
xyr = DATA.elmousept.xyr;
set(DATA.toplevel,'UserData',DATA);
%touch inside ellispe to make cut. If drawing a line, make cut when
%released
if mode == 2  || DATA.elmousept.shape  ==1 
%PCCluster(DATA,DATA.elmousept,1);
       C =  ClassifySpikes(DATA,DATA.elmousept,1);
        if (DATA.plothist | DATA.refitgm) && DATA.elmousept.shape == 1
            [a,b]  = GMDip(C.xy(:,1),0);
            [c,d,e] = GMfit(C.xy,2,1,'idlist',C.clst);
            title(sprintf('E%dP%d mahal %.2f,%.2f(1/2)\n',ei,pi,b.mahal(b.best),d));
            if DATA.plothist
            GetFigure(DATA.tag.hist,'front');
            PlotHist(C.xy, DATA.elmousept);
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

function C = ClassifySpikes(DATA, E, mode)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if isfield(DATA,'NewCut') && DATA.NewCut.probe == p && DATA.NewCut.exptid == e
        C = DATA.NewCut;
    else
        C = Clusters{e}{p};
    end
    C.exptid = e;
    C.cluster = E.cluster;
    if E.cluster > 1
    C.next{E.cluster-1}.shape = E.shape;
    else
    C.shape = E.shape;
    end
    if E.shape == 0 % ellipse
        xy = xyrotate(C.xy(:,1)-E.xyr(1),C.xy(:,2)-E.xyr(2),E.angle);
        r = (xy(:,1)).^2./E.xyr(3).^2 + (xy(:,2)).^2./E.xyr(4).^2;
        id = find(r < 1);
        C.clst(id) = E.cluster+1;
        id = find(r >=1 & C.clst' == E.cluster+1);
        C.clst(id) = 1;
        PlotClusterPoints(C,[]);
        if E.cluster > 1
            n = E.cluster-1;
            C.next{n}.xyr = E.xyr;
            C.next{n}.angle = E.angle;
            C.next{n}.crit = 0;
            DrawEllipse(E,'color',DATA.colors{n+2});
        else
            C.xyr = E.xyr;
            C.angle = E.angle;
            DrawEllipse(E,'color','r');
        end
    elseif E.shape == 1
        if E.mode == 3 %invert sign
            C.sign = C.sign * -1;
        else
            E.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        end
        C.addangle = E.angle;
        xy = xyrotate(C.xy(:,1),C.xy(:,2),E.angle);
        crit = xyrotate(E.pos([1 3]),E.pos([2 4]),E.angle);
        if C.sign < 0
        id = find(xy(:,1) < mean(crit(:,1)));
        nid = find(xy(:,1) > mean(crit(:,1)));
        else
        id = find(xy(:,1) > mean(crit(:,1)));
        nid = find(xy(:,1) < mean(crit(:,1)));
        end
        C.angle = C.angle+E.angle;
        C.clst(id) = 2;
        C.clst(nid) = 1;
        C.crit = mean(crit(:,1));
        C.xy = xy;
        E.pos(1) = mean(crit(:,1));
        E.pos(3) = mean(crit(:,1));
        E.pos(4) = max(xy(:,2));
        E.pos(2) = min(xy(:,2));
        E.pos(5) = max(xy(:,1));
        E.pos(7) = min(xy(:,1));
        E.pos(6) = mean(E.pos([4 2]));
        E.pos(8) = E.pos(6);
%        C.angle = 0;
        PlotClusterPoints(C,[]);
        fprintf('%.2f ',E.pos);
        fprintf('\n');
        DrawEllipse(E);
        DATA.elmousept.pos = E.pos;
        if C.space(1) == 6
        end
    end
    C.id = find(C.clst) ==2;
    Clusters{e}{p}.clst = C.clst;
    Clusters{e}{p}.xy = C.xy;
    Clusters{e}{p}.crit = C.crit;
    setappdata(DATA.toplevel,'Clusters',Clusters);
    DATA.NewCut = C;
%need to record new cluster, and new E.pos
    set(DATA.toplevel,'UserData',DATA);
    
function ScrollWheel(src, evnt)
DATA = GetDataFromFig(src);
DATA.elmousept.angle = DATA.elmousept.angle+0.02*evnt.VerticalScrollCount;
fprintf('Angle %.2f\n',DATA.elmousept.angle);
DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
set(DATA.toplevel,'UserData',DATA);

 function XYButtonDragged(src, data)
DATA = GetDataFromFig(src);
if isfield(DATA,'elmousept') &&  DATA.elmousept.down >0
%    fprintf('D%d,%d %.2f %.2f\n',DATA.elmousept.down,DATA.elmousept.steps,DATA.elmousept.pos(1),DATA.elmousept.pos(3));
    DATA.elmousept.steps = DATA.elmousept.steps +1;
if  DATA.elmousept.down == 1
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
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
end


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
set(DATA.toplevel,'UserData',DATA);
SaveCellList(DATA);
if length(F) == 1
    SetTrialList(DATA);
else
    SpoolSpikes(DATA.spoolfig,'excludelist', trials);
end

function [e,p] = cPoint(DATA)
    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    
function SetTrialList(DATA)
    F = findobj('type','figure','Tag',DATA.tag.spikes);
    Clusters = getappdata(DATA.toplevel,'Clusters');
    C = Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)};
    tn = [DATA.Expt.Trials.Trial]';
    xcl = FindExcludedTrials(DATA, DATA.currentpoint(1),DATA.currentpoint(2), 1, C);
    tn(xcl) = tn(xcl).*-1;
    strial = DATA.currenttrial;
    if strial > length(tn)
        strial = 1;
    end
    
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
    [nr,nc] = Nsubplots(length(eid));
    SetFigure(DATA,DATA.tag.allxy);
    for j = 1:length(eid)
        mysubplot(nr,nc,j);
        hold off; 
        PlotMeanSpike(Clusters{eid(j)}{cid(j)},0,clid(j),type);
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
        text(xl(2),yl(1),sprintf('E%dP%d/%d',eid(j),cid(j),clid(j)),'horizontalalignment','right','verticalalignment','bottom');
    end

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
    plottype = 1;
    T = DATA.Expt.Trials;
    Expt = DATA.Expt;
    setfig = 1;

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'expt',4)
            j = j+1;
            Expt = varargin{j};
            T = Expt.Trials;
        elseif strncmpi(varargin{j},'noplot',6)
            plottype = 0;
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
    if plottype > 0  && setfig
    SetFigure(DATA,DATA.tag.xyseq);
    hold off;
    end
%    need to use time or trial bins for this test, as rate changes 
    if isfield(Expt.Header,'timeoffset')
        toff = Expt.Header.timeoffset;
    else
        toff = 0;
    end

    for j = nc
        xcl = FindExcludedTrials(DATA,e,p,j-1, C);
        id = find(C.clst == j);
        if C.shape == 0
            r = CalcRadius(C);
        else
            r = C.xy(:,1);
        end
        if plottype == 1
            plot(C.times(id)+toff,r(id),'.','color',DATA.colors{j});
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
if plottype > 0
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
end

function HitXYseq(a,b)
    DATA = GetDataFromFig(a);
    SetTrialList(DATA);
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
    if DATA.spikesloaded
        SpoolSpikes(DATA,DATA.currentpoint,'ids',id(end))
    end


function [eid, cid, clid] = FindCell(DATA, cellid)
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);
    [eid, id] = sort(eid);
    cid = cid(id);
    clid = clid(id);

    
function PlotAllCell(DATA, type)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');
    colors = mycolors('spkcolors');
    
    
    id = find(DATA.CellList == cellid);
    [eid, cid, clid] = ind2sub(size(DATA.CellList),id);

    
    [eid, id] = sort(eid);
    cid = cid(id);
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
                set(h(np),'color','r');
                set(h(2),'color',colors{np});
            end
           h = title(sprintf('E%.1fP%d 1D%.1f 2D%.1f',C.exptno,cid(j),C.mahal(4),C.mahal(1)));
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

function SetCellFromLine(a,b, cluster, cell)
    axdata = get(gca,'UserData');
    DATA = get(axdata.toplevel,'UserData');
   
    DATA.CellList(DATA.currentpoint(1),axdata.probe,cluster) = cell;
    PlotCellList(DATA,'showfig');

function SetCellFromSubplot(a,b, cell)
    axdata = get(gca,'UserData');
    DATA = get(axdata.toplevel,'UserData');
   
    DATA.CellList(DATA.currentpoint(1),axdata.probe,DATA.cellcluster) = cell;
    PlotCellList(DATA,'showfig');

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

function PlotAllProbeMean(DATA, type)
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
    cmenu = AddCellContextMenu(DATA);
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
        set(gca,'UIContextMenu',cmenu, 'UserData', axdata);
    end

function cmenu = AddCellContextMenu(DATA)
    cmenu = uicontextmenu;
    for j = 1:20
        uimenu(cmenu,'label',sprintf('Cell %d',j),'Callback',{@SetCellFromSubplot, j});
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
    Clusters = getappdata(DATA.toplevel,'Clusters');
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    Expts = getappdata(DATA.toplevel,'Expts');
    
    eid = DATA.currentpoint(1);
    DATA.Expt = Expts{eid,1};
    [nr,nc] = Nsubplots(length(Clusters{eid}));
     if strcmp(type, 'spoolall')
         DATA.plotspk.allprobes =1;
            SpoolAllProbes(DATA, eid, AllSpikes(eid,:),Clusters{eid});
            return;
     elseif strcmp(type, 'spooleverything')
         DATA.plotspk.allprobes =1;
         for j = 1:length(DATA.exptid)
            stopped = SpoolAllProbes(DATA, j, AllSpikes(j,:),Clusters{j});
            if stopped
                return;
            end
         end
         return;
     elseif strcmp(type, 'spooleverycell')
         for j = 1:length(DATA.exptid)
             cid = isacell(DATA, j, 1:DATA.nprobes);
            stopped = SpoolAllProbes(DATA, j, AllSpikes(j,cid),Clusters{j});
            if stopped
                return;
            end
         end
         return;
     elseif strcmp(type, 'spoolcells')
         cid = isacell(DATA, eid, 1:DATA.nprobes);
            SpoolAllProbes(DATA, eid, AllSpikes(eid,cid),Clusters{eid});
            return;
     end

     if strmatch(type,{'allspks' 'allcellspks'})
         F = SetFigure(DATA,DATA.tag.allspikes);
         subplot(1,1,1);
     else
         SetFigure(DATA,DATA.tag.allxy);
     end

     cmenu = uicontextmenu;
    for j = 1:20
        uimenu(cmenu,'label',sprintf('Cell %d',j),'Callback',{@SetCellFromSubplot, j});
    end
    axdata.toplevel = DATA.toplevel;

     for j = 1:length(Clusters{eid})
        mysubplot(nr,nc,j);
        hold off; 
        if strmatch(type,{'allspks' 'allcellspks'})
            C = Clusters{eid}{j};
            [a,b] = isacell(DATA,eid,j);
            if a > 0 || strcmp(type,'allspks')
            QuickSpikes(DATA, AllSpikes{eid,j},C,'showmean');
            if a
                h = title(sprintf('Cell %s E%.0fP%d 1D%.1f 2D%.1f',sprintf('%d ',b),C.exptno,j,C.mahal(4),C.mahal(1)));
                set(h,'color','b','fontweight','bold');
            else
                h = title(sprintf('E%.0fP%d 1D%.1f 2D%.1f',C.exptno,j,C.mahal(4),C.mahal(1)));
            end
            end
        else
            PlotMeanSpike(Clusters{eid}{j},0,1,type);
        end
        axdata.probe = j;
        set(gca,'Xtick',[],'Ytick',[]);
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid, j});
        set(gca,'UserData',axdata,'UIContextMenu',cmenu);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
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
    
function DATA = PlotAllProbeXY(DATA,varargin)
    oneprobe = 0;
    j =1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'oneprobe',7)
            oneprobe = 1;
        end
        j = j+1;
    end
    
    Clusters = getappdata(DATA.toplevel,'Clusters');
    eid = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    if oneprobe
        np = length(Clusters);
    else
        np = length(Clusters{eid});
    end
    [nr,nc] = Nsubplots(np);
    SetFigure(DATA,DATA.tag.allxy);
    for j = 1:np
        mysubplot(nr,nc,j);
        hold off; 
        if oneprobe
            PlotClusterXY(DATA,Clusters{j}{p},'shorttitle');
            missed = MissedCell(DATA,[j p]);
        else
            PlotClusterXY(DATA,Clusters{eid}{j},'shorttitle');
            missed = MissedCell(DATA,[eid j]);
        end
        set(gca,'Xtick',[],'Ytick',[]);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        if missed 
            set(h,'color','g');
        end
        if size(DATA.CellList,1) >= eid && DATA.CellList(eid,j,1) > 0
            t = text(xl(2),yl(1),sprintf('Cell %d',DATA.CellList(eid,j,1)),'VerticalAlignment','bottom','HorizontalAlignment','right');
            if DATA.plot.density
                set(t,'color','w');
            end
        end
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid, j});
        c = get(gca,'Children');
        for k = 1:length(c)
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid, j});
        end
            
    end
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    
function missed = MissedCell(DATA, pos)    
    e = pos(1);
    p = pos(2);
    fitd = 0;
    if isfield(DATA,'GaussFitdp')
        fitd = DATA.GaussFitdp(e,p,2);
        if sum(DATA.gmfitpos(e,p,:)) < 2
            fitd = 0;
        end
    end
    d = max(DATA.mahal(e,p,[1 3]));
    if isacell(DATA,e,p) == 0 && (d > 2 && fitd > 1) && DATA.peakdiff(e,p) == 0 && DATA.dropi(e,p) > 1
        missed = 1;
    else
        missed = 0;
    end

function HitExptPlot(src, b, type, e)
    DATA = GetDataFromFig(src);    

    DATA.currentpoint(1) = e;
    DATA = ClearSelections(DATA);
    if DATA.plotspks
        PlotAllProbe(DATA,'allspks');
    end
    set(DATA.toplevel,'UserData',DATA);

 function HitXYPlot(src, b, e,p)
    DATA = GetDataFromFig(src);    
    cmenu = get(gca,'UIContextMenu');
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    if bt == 2  && ~isempty(cmenu)
        return;
    end
    if ~isfield(DATA,'selectprobe')
        DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    end
    DATA.selectprobe(e,p) = ~DATA.selectprobe(e,p);
    if DATA.selectprobe(e,p)
        set(gca,'color',[0 0 0]);
        DATA = ShowData(DATA,e,p,'oneprobe');
    else
        set(gca,'color',[1 1 1]);
    end
    DATA.currentpoint = [e p];
    GetFigure(DATA.tag.all);
%    DATA.markh = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),'color','w');
    set(DATA.toplevel,'UserData',DATA);
    
function PlotAllCellXY(DATA)
    Clusters = getappdata(DATA.toplevel,'Clusters');
    it = findobj('Tag','CellNumberId');
    cellid = get(it,'value');

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
        plots = PlotClusterXY(DATA,Clusters{eid(j)}{cid(j)},'shorttitle');
        if clid(j) ~= 1
            np = clid(j)+1;
            set(plots(np),'color','r');
            set(plots(2),'color',colors{np});
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
           
        set(h,'fontweight','normal','fontsize',12,'color',color);
        set(gca,'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        c = get(gca,'Children');
        for k = 1:length(c)
            set(c(k),'ButtonDownFcn',{@HitXYPlot, eid(j), cid(j)});
        end
    end
    
    function h= DrawLine(E,varargin)

x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
if length(E.pos) > 7
    xp = E.pos([5 7]);
    yp = E.pos([6 8]);
else
xp = [mean(x) + diff(y)/4 mean(x)-diff(y)/4];
yp = [mean(y) - diff(x)/4 mean(y)+diff(x)/4];
end
if ishandle(E.h(1)) 
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


function h= DrawEllipse(E,varargin)

if E.shape(1) == 1
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

    GetFigure(DATA.tag.celllist);
    spkcolors = mycolors('spkcolors');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    for j = 1:length(eid)
        GetFigure(DATA.tag.celllist);
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

function SetCellNumber(a,b, fcn)
DATA = GetDataFromFig(a);
Clusters = getappdata(DATA.toplevel,'Clusters');
modifier = 0;
cellid = 0;
if strcmp(fcn,'set')
    it = findobj(get(a,'parent'),'Tag','ClusterModifier');
    if ~isempty(it)
        modifier = get(it,'value');
        if modifier == 4
            DATA.CellDetails.note(DATA.currentpoint(1),DATA.currentpoint(2)) = get(it,'value');
            DATA.CellDetails.notestr = get(it,'string');
        end
        set(it,'value',1);
    end


    e = DATA.currentpoint(1);
    p = DATA.currentpoint(2);
    it = findobj(get(a,'parent'),'Tag','CellNumberId');
    if ~isempty(it)
        cellid = get(it,'value');
        DATA.currentcell = cellid;
        if modifier == 5
            id = find(DATA.CellList == cellid);
            DATA.CellList(id) = 0;
            DATA.CellDetails.Templates(cellid,:) = DATA.templateid;
            id = find(DATA.usepeaks > 0);
            for j = id;
                DATA.CellList(j,DATA.usepeaks(j),DATA.cellcluster) = cellid;
                DATA.CellChanges = cat(1,DATA.CellChanges,[j DATA.usepeaks(j) cellid 1 now]);
            end
            %        id = find(DATA.usepeaks == 0);
            %       deletes = sum(DATA.CellList(id,:) == cellid);
        elseif ismember(modifier, [1 2 3 4])
            id = find(DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),:) == cellid);
            DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),id) = 0; %remove this cell from other clusters this expt
            DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),modifier) = cellid;
            DATA.CellChanges = cat(1,DATA.CellChanges,[DATA.currentpoint(1),DATA.currentpoint(2) cellid modifier now]);
            if modifier == 1
                C = Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)};
            else
                C = Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)}.next{modifier-1};
            end
            if isfield(C,'excludetrialids');
                % find trial #s that match the ids excluded
                DATA.CellDetails.excludetrials{e,p,modifier} = union(DATA.CellDetails.excludetrials{e,p,modifier}, C.excludetrialids);
            end
        end
        SaveCellList(DATA);
        set(DATA.toplevel,'UserData',DATA);
    end
    PlotCellList(DATA,'showfig');
elseif strcmp(fcn,'delete')
    it = findobj(get(a,'parent'),'tag','CellCluster');
    cl = get(it,'value');
    DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),cl) = 0;
    SaveCellList(DATA);
    PlotCellList(DATA, 'showfig');
    set(DATA.toplevel,'UserData',DATA);
elseif strcmp(fcn,'change')
    if DATA.plotallxy;
        PlotAllCellXY(DATA);
        drawnow;
    end
    if DATA.plotspks
        PlotMenu(DATA,[], 'cells', 'allspks');
        drawnow;
    end
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);
    DATA.currentcell = get(a,'value');
    set(DATA.toplevel,'UserData',DATA);
end

function SaveCellList(DATA)
    cellfile = [DATA.name '/CellList.mat'];
    CellList = DATA.CellList;
    CellDetails = DATA.CellDetails;
    CellChanges = DATA.CellChanges;
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
colors = mycolors;
subplot(1,1,1);
set(f,'UserData',DATA.toplevel);
hold off;
if plotmahal
im = zeros(size(DATA.CellList));
id = find(sum(DATA.CellList,3) > 0);
X = squeeze(DATA.mahal(:,:,1));
im(id) = X(id);
end
if ~isfield(DATA,'nclusters')
    for j = 1:size(DATA.CellList,1)
    for k = 1:size(DATA.CellList,2)
        id = find(DATA.CellList(j,k,:) > 0);
        DATA.nclusters(j,k) = length(id);
    end
    end
end
nc = max(DATA.nclusters); %max # clusters in any one expt, for each probe
offset = DATA.clusteroffset;
if plotmahal
    imagesc(im,'buttondownfcn',{@HitImage, 1});
    caxis([0 5]);
else
    for j = 1:size(DATA.CellList,2)*2
        CellIm(:,j) = DATA.CellList(:,ceil(j/2),1+offset);
    end
    id = find(nc >1);
    for j = 1:length(id)
        tid = find(DATA.nclusters(:,id(j)) > 1);
        CellIm(tid,id(j)*2) = DATA.CellList(tid,id(j),offset+2);
    end
        
imagesc([1 size(DATA.CellList,2)],[1 size(DATA.CellList,1)],CellIm,'buttondownfcn',{@HitImage, 3});
end
hold on;
cells = unique(DATA.CellList(:,:,1+offset));
cells = cells(cells > 0);

for j = 1:length(cells)
    [x,y] = find(DATA.CellList(:,:,1+offset) == cells(j));
    for k = 1:length(y)
        if DATA.nclusters(x(k),y(k)) > 1
            plot([y(k)-0.25 y(k)-0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
        else
            plot([y(k) y(k)], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
        end
    end
    [a,b] = min(x);
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
    set(h,'color',colors{j});
end
cells = unique(DATA.CellList(:,:,2+offset));
cells = cells(cells > 0);
for j = 1:length(cells)
    [x,y] = find(DATA.CellList(:,:,2+offset) == cells(j));
    for k = 1:length(y)
            plot([y(k)+0.25 y(k)+0.25], [x(k)-0.5 x(k)+0.5],'-','color',colors{j},'linewidth',2,'buttondownfcn',{@HitImage, 3});
    end
    [a,b] = min(x);
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)));
    set(h,'color',colors{j});
end
if isfield(DATA,'currentpoint')
    h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2),'color','w');
end

iscellim = sum(DATA.CellList,3) > 0;
if DATA.markcell.candidates
    [x,y] = find(isnan(DATA.CellList(:,:,1)));
    for j = 1:length(x)
        DrawBox(x(j),y(j),'color','g');
    end
end
if DATA.markcell.dropi > 0
    [x,y] = find(DATA.dropi < DATA.markcell.dropi & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),'color','m');
    end
end
if DATA.markcell.mahal > 0
    X = max(DATA.mahal(:,:,[1 3]),[],3);
    [x,y] = find(X < DATA.markcell.mahal & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),'color','r');
    end
end
if DATA.markcell.ellipses > 0
    [x,y] = find(DATA.cutshape ==0 & iscellim);
    for j = 1:length(x)
        DrawBox(x(j),y(j),'color','m');
    end
end

if DATA.plotcells.showmahal > 0
    m = (max(DATA.mahal(:,:,[1 3]),[],3) > DATA.plotcells.showmahal);
    dp = (squeeze(DATA.GaussFitdp(:,:,2)) > DATA.plotcells.showfitdp & sum(DATA.gmfitpos,3) == 2);
    peakok = (abs(DATA.peakdiff) < 1);
    dropok = (DATA.dropi > 1);
    [aid, bid] = find(m & peakok & dropok);
    for j = 1:length(aid)
        DrawBox(aid(j),bid(j),'color','r');
    end
    [cid, did] = find(dp & peakok & dropok);
    for j = 1:length(cid)
        DrawBox(cid(j),did(j),'color','g');
    end
    [aid, bid] = find(dp & m & peakok & dropok);
    for j = 1:length(aid)
        DrawBox(aid(j),bid(j),'color','w');
    end
end
set(gcf, 'KeyPressFcn',{@KeyPressed,3});

function DATA = InitInterface(DATA)

    scrsz = get(0,'Screensize');
    cntrl_box = figure('Position', [10 scrsz(4)-310 360 200],...
        'NumberTitle', 'off', 'Tag',DATA.tag.top,'Name',DATA.name,'menubar','none');
    DATA.toplevel = cntrl_box;
    listh = 0.3;
    lst = uicontrol(gcf, 'Style','listbox','String', DATA.strings,...
            'Callback', {@PlotClusters, 'setentry'},...
            'units','norm', 'Position',[0.01 0 0.98 listh]);
    DATA.lstui = lst;
    nr = 12;
    nc = 4;
    bp = [0.01 listh+0.1./nr 1./(nc+0.1) 1./nr];
    bp(1) = 0.01;
    bp(2) = bp(2)+ 0.5./nr;
    bp(3) = 0.1;
    uicontrol(gcf,'style','pushbutton','string','<', ...
        'Callback', {@NextButton, 'l'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)*2+0.01;
    uicontrol(gcf,'style','pushbutton','string','>', ...
        'Callback', {@NextButton, 'r'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = listh+ 0.1./nr;
    bp(1) = bp(3)+0.01;
    uicontrol(gcf,'style','pushbutton','string','\/', ...
        'Callback', {@NextButton, 'd'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = bp(2)+1./nr;
    uicontrol(gcf,'style','pushbutton','string','^', ...
        'Callback', {@NextButton, 'u'}, 'Tag','NextButton',...
        'units', 'norm', 'position',bp,'value',1);
    bp(2) = bp(2)+1./nr;
    bp(1) = 0.01;
    bp(3) = 1./nc;
    uicontrol(gcf,'style','pushbutton','string','spool', ...
        'Callback', {@SpoolCurrentSpikes}, 'Tag','SpoolButton',...
        'units', 'norm', 'position',bp,'value',1);


    bp(2) = listh+ 0.1./nr;
    bp(1) = 0.3;
    bp(3)=1./nc;
    uicontrol(gcf,'style','pop','string','quality|times|xcorr|ptscatter|mahaldip|bestspace|bestdip|bmcmahal|cov(shape)|xc(shape)', ...
        'Callback', {@PlotClusters, 'setplot'}, 'Tag','plottype',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string','mahal|dips|mahaln|mahal1-n|mahal+var|drop-mahal|Fit-mahal|muamp-spkvar|xcorr|probexcorr|probeco|spkshapecorr|spkshape|spksize|Evec|BuildTimes|PcGms|Exclusions|Tagged|TriggerFind|CandidateCells|SpkRate|EventRate', ...
        'Callback', {@PlotAllClusters}, 'Tag','plotalltype',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','pop','string','plain|probe|cutspace', ...
        'Callback', {@Update}, 'Tag','ColorScheme',...
        'units', 'norm', 'position',bp,'value',1);
    bp(1) = 0.01;
    bp(1) = 0.3;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','spool', ...
        'Callback', {@Update}, 'Tag','SpoolSpikes',...
        'units', 'norm', 'position',bp,'value',DATA.spoolspikes);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','spk xy', ...
        'Callback', {@Update}, 'Tag','SpikeXY',...
        'units', 'norm', 'position',bp,'value',DATA.showspkxy);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','spkmean', ...
        'Callback', {@Update}, 'Tag','SpikeMean',...
        'units', 'norm', 'position',bp,'value',DATA.showspkmean);
    bp(1) = 0.01;
    bp(1) = 0.3;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','Histogram', ...
        'Callback', {@Update}, 'Tag','PlotHistogram',...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','refit GM', ...
        'Callback', {@Update}, 'Tag','RefitGM',...
        'units', 'norm', 'position',bp,'value',DATA.refitgm);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','PlotExpt', ...
        'Callback', {@Update}, 'Tag','PlotExpt',...
        'units', 'norm', 'position',bp,'value',DATA.plotexpt);
    bp(1) = 0.01;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','checkbox','string','AllXY', ...
        'Callback', {@Update}, 'Tag','PlotAllXY',...
        'units', 'norm', 'position',bp,'value',DATA.plotallxy);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','xyseq', ...
        'Callback', {@Update}, 'Tag','PlotXYseq',...
        'units', 'norm', 'position',bp,'value',DATA.plotxyseq);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','TrigHist', ...
        'Callback', {@Update}, 'Tag','PlotTrigHist',...
        'units', 'norm', 'position',bp,'value',DATA.plottrighist);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','Quickspks', ...
        'Callback', {@Update}, 'Tag','PlotSpks',...
        'units', 'norm', 'position',bp,'value',DATA.plotspks);

    
    
    bp(1) = 0.01;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','text','string','Exclude:', ...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','Peak Mismatch', ...
        'Callback', {@Update}, 'Tag','PeakMismatch',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.offpeak);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','Peak Match', ...
        'Callback', {@Update}, 'Tag','PeakMatch',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.onpeak);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','NonCell', ...
        'Callback', {@Update}, 'Tag','NonCell',...
        'units', 'norm', 'position',bp,'value',DATA.exclude.noncell);
    
    bp(1) = 0.01;
    bp(3) = 0.2;
    bp(2) = bp(2)+ 1./nr;
    uicontrol(gcf,'style','text','string','Label', ...
        'units', 'norm', 'position',bp,'value',DATA.plothist);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','Ex#', ...
        'Callback', {@Update}, 'Tag','LabelExptno',...
        'units', 'norm', 'position',bp,'value',DATA.show.exptno);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','Name', ...
        'Callback', {@Update}, 'Tag','LabelExptName',...
        'units', 'norm', 'position',bp,'value',DATA.show.exptname);
    bp(1) = bp(1)+bp(3)+0.01;
    uicontrol(gcf,'style','checkbox','string','ed', ...
        'Callback', {@Update}, 'Tag','LabelEd',...
        'units', 'norm', 'position',bp,'value',DATA.show.ed);
    hm = uimenu(cntrl_box,'Label','File','Tag','XClusterMenu');
    uimenu(hm,'Label','Load All','Callback',{@LoadAll, 'force'});
    uimenu(hm,'Label','ReLoad New','Callback',{@LoadAll, 'ifnew'});
    uimenu(hm,'Label','Use Current as Template','Callback',{@PlotClusters, 'followcorr'});
    uimenu(hm,'Label','Reexamine raw V','Callback',{@PlotClusters, 'CallAllVPcs'});
    uimenu(hm,'Label','Close','Callback',{@PlotClusters, 'close'});
    hm = uimenu(cntrl_box,'Label','Options','Tag','OptionMenu');
    uimenu(hm,'Label','Density','Callback',{@OptionMenu, 1});
    uimenu(hm,'Label','2D GM Clustering','Callback',{@OptionMenu, 8});
    uimenu(hm,'Label','Plot Cell','Callback',{@OptionMenu, 2});
    uimenu(hm,'Label','Plot AllXY This Expt','Callback',{@OptionMenu, 4});
    uimenu(hm,'Label','Plot AllXY This Probe','Callback',{@OptionMenu, 11});
    uimenu(hm,'Label','Plot MeanSpike (line) This Expt','Callback',{@OptionMenu, 5});
    uimenu(hm,'Label','Plot MeanSpike(IM) This Expt','Callback',{@OptionMenu, 6});
    uimenu(hm,'Label','Plot Expt All Probes','Callback',{@OptionMenu, 3});
    uimenu(hm,'Label','Exclude Selected Trials','Callback',{@OptionMenu, 7});
    uimenu(hm,'Label','Exclude Trials Cl2','Callback',{@OptionMenu, 'excludec12'});
    uimenu(hm,'Label','clear exclusions Cl1','Callback',{@OptionMenu, 'usealltrialscl1'});
    uimenu(hm,'Label','clear exclusions Cl2','Callback',{@OptionMenu, 'usealltrialscl2'});
    uimenu(hm,'Label','xcorr selected','Callback',{@OptionMenu, 9});
    uimenu(hm,'Label','xcorr all','Callback',{@OptionMenu, 10});
    uimenu(hm,'Label','xcorr cells','Callback',{@OptionMenu, 'xcorrcells'});
    uimenu(hm,'Label','xcorr cell all Expts','Callback',{@OptionMenu, 'xcorrallcells'});
    uimenu(hm,'Label','xcorr all Expts/Probes','Callback',{@OptionMenu, 'xcorrallprobes'});
    uimenu(hm,'Label','Find Candidate Connections','Callback',{@OptionMenu, 'findconnect'});
    uimenu(hm,'Label','Smoothed Image','Callback',{@OptionMenu, 12});
    uimenu(hm,'Label','CellFinder','Callback',{@OptionMenu, 13});
    sm = uimenu(hm,'Label','Layout');
    uimenu(sm,'Label','Save','Callback',{@OptionMenu, 'savelayout'});
    uimenu(sm,'Label','Load','Callback',{@OptionMenu, 'loadlayout'});
    uimenu(sm,'Label','Windows to Front','Callback',{@OptionMenu, 'tofront'});
    sm = uimenu(hm,'Label','Spooling','Tag','SpoolOptions');
    uimenu(sm,'Label','Used Saved Cluster Codes','Callback',{@OptionMenu, 'usesavedcodes'});
    uimenu(sm,'Label','Show Excluded Trials','Callback',{@OptionMenu, 'showexcluded'});
    uimenu(sm,'Label','Superimpose Cell Meanss','Callback',{@OptionMenu, 'showcellmeans'});

    
    hm = uimenu(cntrl_box,'Label','Plots','Tag','PlotMenu');
    sm = uimenu(hm,'Label','Expts','Tag','ExptPlotMenu');
    uimenu(sm,'Label','Means','Callback',{@PlotMenu, 'expt', 'means'});
    uimenu(sm,'Label','Trial Counts','Callback',{@PlotMenu, 'expt', 'trialcounts'});
    uimenu(sm,'Label','None','Callback',{@PlotMenu, 'expt', 'none'});
    sm = uimenu(hm,'Label','Cells','Tag','CellPlotMenu');
    uimenu(sm,'Label','AllXY','Callback',{@PlotMenu, 'cells', 'AllXY'});
    uimenu(sm,'Label','AllMean','Callback',{@PlotMenu, 'cells', 'AllMean'}); 
    uimenu(sm,'Label','AllMeanIm','Callback',{@PlotMenu, 'cells', 'AllMeanIm'}); 
    uimenu(sm,'Label','AllTrigHist','Callback',{@PlotMenu, 'cells', 'TrigHist'}); 
    uimenu(sm,'Label','Spool All - This expt','Callback',{@PlotMenu, 'cells', 'spoolall'}); 
    uimenu(sm,'Label','Spool All - All expts','Callback',{@PlotMenu, 'cells', 'spooleverycell'}); 
    uimenu(sm,'Label','Spool selected, all expts','Callback',{@PlotMenu, 'cells', 'spoolone'}); 
    uimenu(sm,'Label','AllSpks Cell','Callback',{@PlotMenu, 'cells', 'allspks'});
    uimenu(sm,'Label','All Cells Spks Expt','Callback',{@PlotMenu, 'probes', 'cellspks'});
    uimenu(sm,'Label','Spool All Cells','Callback',{@PlotMenu, 'cells', 'spoolall'});
    sm = uimenu(hm,'Label','Probes','Tag','ProbePlotMenu');
    uimenu(sm,'Label','Spool All Probes','Callback',{@PlotMenu, 'probes', 'spoolall'});
    uimenu(sm,'Label','Spool All Probes All Expts','Callback',{@PlotMenu, 'probes', 'spooleverything'});
    uimenu(sm,'Label','Spool Probes With Cells','Callback',{@PlotMenu, 'probes', 'spoolcells'});
    uimenu(sm,'Label','All Spks This Expt','Callback',{@PlotMenu, 'probes', 'allspks'});
    uimenu(sm,'Label','All Cell Spks This Expt','Callback',{@PlotMenu, 'probes', 'allcellspks'});
    uimenu(sm,'Label','All XY This Expt','Callback',{@PlotMenu, 'probes', 'AllXY'});
    uimenu(sm,'Label',' Mean Spikes','Callback',{@PlotMenu, 'probes', 'AllMean'});
    uimenu(sm,'Label','Mean Spikes Image','Callback',{@PlotMenu, 'probes', 'AllMeanIm'});
    uimenu(sm,'Label','Mean Spike All Probes and Expts','Callback',{@PlotMenu, 'probes', 'AllExptMeanIm'});
    uimenu(sm,'Label','Spike dprime All Probes and Expts','Callback',{@PlotMenu, 'probes', 'AllExptdpIm'});
    sm = uimenu(hm,'Label','MeanSpike','Tag','ProbePlotMenu');
    uimenu(sm,'Label','Image','Callback',{@PlotMenu, 'mean', 'imageonly'});
    uimenu(sm,'Label','Image+lines','Callback',{@PlotMenu, 'mean', 'DpLines'});
    uimenu(sm,'Label','Mean Lines','Callback',{@PlotMenu, 'mean', 'meanlines'});
    uimenu(sm,'Label','Lines (with dprime)','Callback',{@PlotMenu, 'mean', 'lineonly'});
    uimenu(sm,'Label','Image','Callback',{@PlotMenu, 'mean', 'MeanIm'});
    sm = uimenu(cntrl_box,'Label','Trials','Tag','TrialMenu');
    uimenu(sm,'Label','Exclude Selected Trials','Callback',{@OptionMenu, 7});
    uimenu(sm,'Label','Exclude Trials Cl2','Callback',{@OptionMenu, 'excludec12'});
    uimenu(sm,'Label','clear exclusions Cl1','Callback',{@OptionMenu, 'usealltrialscl1'});
    uimenu(sm,'Label','clear exclusions Cl2','Callback',{@OptionMenu, 'usealltrialscl2'});
    uimenu(sm,'Label','Show Excluded Trials','Callback',{@OptionMenu, 'showexcluded'});
    sm = uimenu(cntrl_box,'Label','Tag','Tag','TagMenu');
    uimenu(sm,'Label','Possible Cell','Callback',{@TagMenu, '?cell'});
    uimenu(sm,'Label','AllVPcs - Threshold','Callback',{@TagMenu, 'threshold'});
    uimenu(sm,'Label','AllVPcs  > 1 cell','Callback',{@TagMenu, 'morecells'});
    uimenu(sm,'Label','AllVPCs - improve','Callback',{@TagMenu, 'improve'});
    uimenu(sm,'Label','AllVPCs - error','Callback',{@TagMenu, 'error'});
    %    hm = uimenu(sm,'Label','Comment','Tag','TagMenu');
    uimenu(sm,'Label','Add Comment Manually','Callback',{@TagMenu, 'comment'});
    uimenu(sm,'Label','Comment: Poor Stability','Callback',{@TagMenu, 'poor stability'});
    uimenu(sm,'Label','Comment: Poor Isolation','Callback',{@TagMenu, 'poor isolation'});
    uimenu(sm,'Label','Comment: Dropping Spikes','Callback',{@TagMenu, 'dropping spikes'});
    uimenu(sm,'Label','Comment: Same As another probe','Callback',{@TagMenu, 'repeat probe'});
    uimenu(sm,'Label','Clear Selected','Callback',{@TagMenu, 'clear'});
    DATA.tagstrings = {'?cell' 'morecells' 'threshold' 'improve' 'error', 'comment' 'poor stability' 'poor isolation' 'dropping spikes' 'clear'};

    DATA.fig.top = DATA.toplevel;
    DATA.fig.clusters = GetFigure(DATA.tag.clusters);
    tmp.parentfig = DATA.toplevel;
    set(DATA.fig.clusters,'Name','Clusters','UserData',tmp);
    set(cntrl_box,'UserData',DATA);
    
function DATA = LoadXcorrFiles(DATA)    
    for j = 1:length(DATA.strings)
    end
    
function DATA = LoadAll(a,b, type)
    useman = 1;
    loadexpts = 1;
    AllExpts = {};
    DATA = GetDataFromFig(a);
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

    %need to think about cases where there was just the autocut and a 
    %manual cut is created
    if strcmp(type,'ifnew') %Already set
        d = dir(DATA.name);
        for j  = 1:length(d)
            if strfind(d(j).name,'ClusterTimes.mat') & isempty(strfind(d(j).name,'OldClusterTimes.mat'))
               id = strmatch(d(j).name,DATA.strings);
               if isempty(id) %new file
                   str = strrep(d(j).name, 'ClusterTimes', 'AutoClusterTimes');
                   id = strmatch(str,DATA.strings);
                   if length(id) == 1
                       DATA.strings{id} = d(j).name;
                       AutoClusters{id} = Clusters{id};
                       [C, details] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                       for k = 1:length(C)
                           if ~isempty(C{k})
                               Clusters{id}{k} = C{k};
                           end
                       end
                   end
               end
            end
        end
        cid = 1:length(DATA.strings);
        loadexpts = 0;
        mid = [];
    elseif useman == 0
        DATA.strings = DATA.strings(find(AutoExpt));
        cid = 1:length(DATA.strings);
        mid = [];
    elseif length(ManExpt)
        %cid find expts where AutoExpt is 0 (= ManExpt set), or AutoExpt is
        %not in ManExpt;
        cid = find(~ismember(AutoExpt,unique(ManExpt(ManExpt > 0))));
        cid = setdiff(cid, find(NewExpt > 0));
        
        mid = find(ManExpt > 0);
    else
        cid = 1:length(AutoExpt);
        mid = [];
    end
    msgmode = 0;
    reloaded = [];
    for j = 1:length(cid)
        k = cid(j);
        if strfind(DATA.strings{k},'Times.mat') && isempty(strfind(DATA.strings{k},'NewClusterTimes'))
            d = dir([DATA.name '/' DATA.strings{k}]);
            if d.datenum > DATA.dates(k) || strcmp(type,'force')
                DATA.dates(j) = d.datenum;
                if ismember(k,mid)
                    autocl = 0;
                    [AutoClusters{j}, autodetails] = LoadClusters([DATA.name '/' strrep(DATA.strings{k},'Cluster','AutoCluster')]);
                else %if no manual cluster, use the AutoCluster as main....
                    AutoClusters{j} = [];
                    autocl = 1;
                end

      
                fprintf('Loading %s\n',DATA.strings{k});
                if strcmp(type,'ifnew') %Only Change modified Clusters, so hat don't lose bits from AutoClusters
                    [C, details] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                    for k = 1:length(C)
                        if isfield(C{k},'xy') & C{k}.ctime > Clusters{j}{k}.ctime
                            Clusters{j}{k} = C{k};
                            reloaded(j,k)=1;
                        end
                    end
                else
                    [Clusters{j}, details] = LoadClusters([DATA.name '/' DATA.strings{k}]);
                end
                np = max([np length(Clusters{j})]);
                exptno(j) = details.exptno;
                if isfield(details,'FullVData')
                    FullVData{j} = details.FullVData;
                else
                    fprintf('No FullV Data for %s\n',DATA.strings{k});
                end
                for k = 1:length(Clusters{j})
                    if ~isfield(Clusters{j}{k},'mahal')
                        s = sprintf('Cluster %d Ex %d missing',k,j);
                        if msgmode == 0
                            s = questdlg(sprintf('Cluster %d Ex %d missing',k,j),'ClusterError','OK','Ignore Others','Ignore Others');
                            if strcmp(s,'Ignore Others')
                                msgmode = 1;
                            end
                        else
                            errordlg(sprintf('Cluster %d Ex %d missing',k,j),'ClusterError','modal');
                            fprintf('Cluster %d Ex %d missing\n',k,j);
                        end
                        Clusters{j}{k} = AutoClusters{j}{k};
                        Clusters{j}{k}.auto = 1;
                    elseif length(AutoClusters{j}) >= k && AutoClusters{j}{k}.savetime(1) > Clusters{j}{k}.savetime(1) && Clusters{j}{k}.auto == 1
                        Clusters{j}{k} = AutoClusters{j}{k};
                    end
                    if autocl
                        Clusters{j}{k}.auto = 1;
                    end
                    if ~isfield(Clusters{j}{k},'clst') && ~isfield(AutoClusters{j}{k},'clst')
                        fprintf('Missing clst for %d,%d\n',exptno(j),k);
                    end
                    if ~isfield(Clusters{j}{k},'xy') && isfield(AutoClusters{j}{k},'xy')
                        Clusters{j}{k}.xy = AutoClusters{j}{k}.xy;
                        if ~isfield(Clusters{j}{k},'clst') 
                            Clusters{j}{k}.clst = AutoClusters{j}{k}.clst;
                            Clusters{j}{k}.times = AutoClusters{j}{k}.times;
                        end
                        %next test should not be true, but is sometimes. Eg
                        %211 Ex13 P24.  clst should not be in Clusters.
                        if length(Clusters{j}{k}.times) < length(Clusters{j}{k}.clst) && Clusters{j}{k}.auto == 1
                            Clusters{j}{k}.clst = AutoClusters{j}{k}.clst;
                            Clusters{j}{k}.times = AutoClusters{j}{k}.times;
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
                end
                smrname = regexprep(DATA.name,'lem/M([0-9]*)(.*)','$0/lemM$1');
                exfile = [smrname '.' sprintf('%d',floor(exptno(j))) 'idx.mat'];
                if exist(exfile,'file') && loadexpts > 0
                    fprintf('Loading %s\n',exfile);
                    [Trials, Expts] = APlaySpkFile(exfile,'noerrs');
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
                    for k = 1:length(Clusters{j})
                        C = Clusters{j}{k};
                        if isfield(C,'clst')
                            spkt = C.times(C.clst>1) .*10000;
                        else
                            spkt = C.times .*10000;
                        end
                        for t = 1:length(Expt.Trials)
                            id = find(spkt > Expt.Trials(t).Start(1) & spkt < Expt.Trials(t).End(end));
                            Expt.Trials(t).Spikes = round(spkt(id)'-Expt.Trials(t).Start(1));
                        end
                        Expt.Header.expname = Expt2Name(Expt);
                        AllExpts{j,k} = Expt;
                    end
                end
            else
                DATA.electrodedepth(j) = NaN;
                exptno(j) = Clusters{j}{1}.exptno;
            end
            end
    end
    for j = 1:length(Clusters)
        Clusters{j} = FixClusters(Clusters{j});
    end
    DATA = LoadComments(DATA);
   if strcmp(type,'ifnew')     && length(reloaded) && DATA.spikesloaded
      [a,b] = find(reloaded);
      LoadSelectedSpikes(DATA,a,b);
   end
    [exptno, id] = sort(exptno);
    DATA.strings = DATA.strings(cid(id));
    DATA.dates = DATA.dates(id);
    DATA.exptid = exptno;
    DATA.trialids = DATA.trialids(id);
    if length(FullVData) == length(id)
        FullVData = FullVData(id);
    end
    Clusters = Clusters(id);
    if loadexpts && length(AllExpts)
        Expts = AllExpts(id,:);
        for j = 1:size(Expts,1)
            DATA.expnames{j} = Expts{j,1}.Header.expname;
        end
        setappdata(DATA.toplevel,'Expts',Expts);
    DATA.electrodedepth = DATA.electrodedepth(id);
    end
   if ~strcmp(type,'ifnew')   
    DATA.nprobes = np;
    DATA.voffset = [1:DATA.nprobes] .*2;

   end

    set(DATA.toplevel,'UserData',DATA);
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    set(DATA.toplevel,'UserData',DATA);
    setappdata(DATA.toplevel,'Clusters',Clusters);
    setappdata(DATA.toplevel,'FullVData',FullVData);
    DATA = LoadCellFile(DATA);
    DATA = LoadExtra(DATA);
    fprintf('Done\n');

    if strcmp(type,'loadspikes')
        fprintf('Loading Spikes...');
        LoadAllSpikes(DATA,'quiet');
        DATA.spikesloaded = 1;
        set(DATA.toplevel,'UserData',DATA);
        fprintf('\n');
    end

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
    
function CheckExpts(DATA, type)
    p = 1;
    Expts = getappdata(DATA.toplevel,'Expts');
 for j = 1:length(Expts)
     if strcmp(type,'errs')
         for k = 1:length(Expts{j,p}.Header.errs)
             fprintf('%d: %s\n',j,Expts{j,p}.Header.errs{k});
         end
     elseif strcmp(type,'method')
         fprintf('Ex %d Method %d\n',j,Expts{j,1}.Header.ReadMethod);
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
GetFigure(DATA.tag.clusters);
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
         end
     end
 end

function DATA = ConvertExclusion(DATA)
    
    
    Expts = getappdata(DATA.toplevel,'Expts');
    if size(DATA.CellDetails.excludetrials,1) < length(Expts)
        DATA.CellDetails.excludetrials{length(Expts),DATA.nprobes,6} = [];
    end
    if isfield(DATA.CellDetails,'version') & DATA.CellDetails.version > 1.1
        CheckExclusion(DATA);
        return;
    end
    for j = 1:size(DATA.CellDetails.excludetrials,1)
        for k = 1:size(DATA.CellDetails.excludetrials,2)
            for c = 1:size(DATA.CellDetails.excludetrials,3)
                X = DATA.CellDetails.excludetrials{j,k,c};
                if length(X)
                    if max(X) > length(Expts{k,j}.Trials)
                        fprintf('Too Many Excluded E%dP%dC%d: %d vs %d\n',j,k,c,length(X),length(Expts{k,j}.Trials));
                        X = X(X < length(Expts{k,j}.Trials));
                        ids = [Expts{k,j}.Trials(X).id];
                        DATA.CellDetails.excludetrials{j,k,c} = ids;
                    else
                        ids = [Expts{k,j}.Trials(X).id];
                        DATA.CellDetails.excludetrials{j,k,c} = ids;
                    end
                end
            end
        end
    end
    DATA.CellDetails.version = DATA.version;
%    set(DATA.toplevel,'UserData',DATA);
    

function DATA = CheckExclusion(DATA)
    
    nerr = 0;
    Expts = getappdata(DATA.toplevel,'Expts');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    for j = 1:size(DATA.CellDetails.excludetrials,1)
        for k = 1:size(DATA.CellDetails.excludetrials,2)
            for c = 1:size(DATA.CellDetails.excludetrials,3)
                X = DATA.CellDetails.excludetrials{j,k,c};
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

function name = FileName(DATA, ex, type)
    
   smrname = regexprep(DATA.name,'lem/M([0-9]*)','$0/lemM$1');
   if strmatch(type,'Expt')
       name = [smrname '.' num2str(ex) 'idx.mat'];
   elseif strmatch(type,'FullV')
       name = [DATA.name '/Expt' num2str(ex) 'FullV.mat'];
   end

function PlotShapes(DATA)    
spk = 1;
ids = [];
DATA.nprobes = 24;
GetFigure(DATA.tag.clusters);
clf;
for j = 1:DATA.nprobes
    subplot(4,6,j);
    Shape(DATA,j);
    drawnow;
end
    
function CompareShapes(DATA, type)    
spk = 1;
ids = [];
DATA.nprobes = 24;
GetFigure(DATA.tag.clusters);
clf;
for j = 1:DATA.nprobes
    subplot(4,6,j);
    if type == 1
    CompareShape(DATA,j);
    elseif type == 2
    PlotShape(DATA,j);
    end
    drawnow;
end


function SpikeDistances(DATA, type)
spk = 1;
ids = [];
GetFigure('Clusters','front');
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
GetFigure('Clusters','front');
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
GetFigure('Clusters','front');
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
GetFigure('Eigenvectors','front')
hold off; 
colors = mycolors;
for j = 0:5
    plot(E(:,end-j),'color',colors{j+1});
    hold on;
end

function PlotShape(DATA, spk)
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

function CellFinder(DATA)
        
mid = DATA.mahaltype;
%First make a smoothed map of mahal distances, and find local maxima
X = squeeze(DATA.mahal(:,:,mid));
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
GetFigure(DATA.tag.templatesrc);
PlotMeanSpike(Clusters{ex}{spk},spk,1);
GetFigure(DATA.tag.clusters);

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
set(gcf,'KeyPressFcn',{@KeyPressed, 2});
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
    elseif mahald(j,peaks(j)) > DATA.mahalcrit && xc(j,peaks(j)) > DATA.xccrit
        plot(square(1,:)+peaks(j),square(2,:)+j,'w-','buttondownfcn',{@HitImage, 1},'linewidth',2);
        DATA.usepeaks(j) = peaks(j);
    else
        DATA.usepeaks(j) = 0;
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


function h = DrawBox(ex, p, varargin)
square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];
h = plot(square(1,:)+p,square(2,:)+ex,'k-','buttondownfcn',{@HitImage, 1},'linewidth',2);
if length(varargin)
    set(h,varargin{:});
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
    

function CompareShape(DATA, spk)
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


function h = MarkCurrentCluster(DATA)
if strmatch(DATA.plot.alltype,{'template' 'mahal+var'})
    SetFigure(DATA, DATA.tag.all);
    square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];
    hold on;
    if ishandle(DATA.markh)
        delete(DATA.markh);
    end
    h = plot(square(1,:)+DATA.currentpoint(2),square(2,:)+DATA.currentpoint(1),'k-','linewidth',2);
else 
    h = [];
end

function DATA = ClearSelections(DATA, force)
%wipe the list of selected probes
%Unless force is set, this may depend on GUI settings/state
    DATA.selectprobe = zeros(length(DATA.exptid),DATA.nprobes);



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
    DATA = ClearSelections(DATA);
elseif dir == 'l' %step 1 left
    if DATA.currentpoint(2) > 1
    DATA.currentpoint(2) = DATA.currentpoint(2)-1;
    else
        DATA.currentpoint(2) = DATA.nprobes;
        if DATA.currentpoint(1) > 1
        DATA.currentpoint(1) = DATA.currentpoint(1)-1;
        end
    end
    DATA = ClearSelections(DATA);
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
            DATA.currentpoint(1) = DATA.nprobes;
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
    DATA = ClearSelections(DATA);
end
ShowData(DATA,DATA.currentpoint(1),DATA.currentpoint(2),'oneprobe');
if DATA.plotallxy
    PlotAllProbeXY(DATA);
    drawnow;
end
if ismember(dir,'ud')  %moving expts up and down with arrow keys, spoole/show all probes
    e = DATA.currentpoint(1);
    PrintComments(DATA,e,1:DATA.nprobes);
    if DATA.spoolspikes
        PlotAllProbe(DATA,'spoolall');
    elseif DATA.plotspks
        PlotAllProbe(DATA,'allspks');
    end
    DATA = ClearSelections(DATA);
    str = sprintf('Ex %.1f: %s',DATA.exptid(DATA.currentpoint(1)),DATA.expnames{DATA.currentpoint(1)});
    fprintf('%s\n',str);
    SetFigure(DATA, DATA.tag.all);
    title(str);
    Expts = getappdata(DATA.toplevel,'Expts');
    DATA.Expt = Expts{e,1};
    SetTrialList(DATA);
end
stopped = 0;
set(DATA.toplevel,'UserData',DATA);


function DATA = SpoolCurrentSpikes(mn,b)
    DATA = GetDataFromFig(mn);
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');

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
        PlotClusterXY(DATA,Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)});
        drawnow;
        SetTrialList(DATA);
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
    details = [];
    x = [];
    if ~isfield(DATA,'trialids') || cluster < 1
        return;
    end
    if isfield(C,'excludetrialids')
        x = find(ismember(DATA.trialids{e},C.excludetrialids));
    end
    details.nt = length(DATA.trialids{e});
    details.fraction = 0;
    
    XC = DATA.CellDetails.excludetrials;
        
    if size(XC,1) < e || size(XC,2) < p || size(XC,3) < cluster
        return;
    end
    y = find(ismember(DATA.trialids{e},XC{e,p, cluster}));
    x = [x y];
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
GetFigure('smoothed');
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


function handles = PlotPopPoints(X,Y, varargin)

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
elseif strmatch(colorscheme,'cutspace')
    c = DATA.cutspace;
    fill = zeros(size(c));
else
    c = ones(size(X));
    fill = zeros(size(c));
end
typelist = unique(types);
hold off;
for j = 1:size(X,1)
    for k = 1:size(X,2)
        if include(j,k)
            if isempty(shapelist)
                shape = 'o';
            else
                shape = shapelist(j,k);
            end
        h = plot(X(j,k),Y(j,k),shape,'buttondownfcn',{@HitPopPoint, j,k},'color',colors{c(j,k)});
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

function PlotMahalImage(DATA,varargin)
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
imagesc(squeeze(DATA.mahal(:,:,mid)),'buttondownfcn',{@HitImage, 1});
set(gcf, 'KeyPressFcn',{@KeyPressed, 1});
 
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
    plot(cat(1,j, DATA.peakpos(:,j)),cat(2,0.5,1:size(DATA.peakpos,1)),'w-','buttondownfcn',{@HitImage, 1});
    hold on;
end
for j = 1:size(DATA.mahal,1)
    text(1,j,ExLabel(DATA,j),'color','w');
end
 set(gca,'UserData',DATA.toplevel);
 set(gcf,'UserData',DATA.toplevel);
    
 
  function KeyPressed(src, ks, fcn)

DATA = GetDataFromFig(src);
tag = get(src,'tag');

if strmatch(ks.Key,'delete') 
    if fcn == 2 % delete box marking cell on template track plot
        DATA.usepeaks(DATA.currentpoint(1)) = 0;
        if ishandle(DATA.markcc)
            delete(DATA.markcc)
        end
        h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2));
        set(h,'color','r');
        set(DATA.toplevel,'UserData',DATA);
    elseif fcn == 3  %delete current location from CellList
        it = findobj(src,'Tag','CellNumberId');
        if length(it) == 1
            cellid = get(it,'value');
            id = find(ismember(DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),:),cellid));
            if length(id) == 1
                DATA.cellcluster = id;
            end
        end

        DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),DATA.cellcluster) = 0;
        PlotCellList(DATA,'showfig');
        SaveCellList(DATA);
        set(DATA.toplevel,'UserData',DATA);
    end
elseif strmatch(ks.Key,'add')
    if strcmp(tag,DATA.tag.celllist)
        id = find(DATA.CellList(DATA.currentpoint(1),:,DATA.cellcluster) == DATA.currentcell);
        DATA.CellList(DATA.currentpoint(1),id) = 0;
        DATA.CellList(DATA.currentpoint(1),DATA.currentpoint(2),DATA.cellcluster) = DATA.currentcell;
        DATA.markh = MarkCurrentCluster(DATA);
        PlotCellList(DATA,'showfig');
    else
        if DATA.usepeaks(DATA.currentpoint(1)) > 0
            h = DrawBox(DATA.currentpoint(1),DATA.usepeaks(DATA.currentpoint(1)));
            set(h,'color','r');
        end
        DATA.usepeaks(DATA.currentpoint(1)) = DATA.currentpoint(2);
    end
    h = DrawBox(DATA.currentpoint(1),DATA.currentpoint(2));
    set(h,'color','w');
    set(DATA.toplevel,'UserData',DATA);
    
elseif strmatch(ks.Key,'rightarrow')
    DATA = NextButton(src, ks, 'r');
elseif strmatch(ks.Key,'leftarrow')
    DATA = NextButton(src, ks, 'l');
elseif strmatch(ks.Key,'downarrow')
    DATA =  NextButton(src, ks, 'd');
elseif strmatch(ks.Key,'uparrow')
    DATA = NextButton(src, ks, 'u');
    n =1;
elseif strmatch(ks.Key,'add')
elseif strmatch(ks.Key,'subtract')
elseif strmatch(ks.Key,'space')
end
    if strcmp(tag,DATA.tag.celllist)
        figure(src);
    end

    
function s = ExLabel(DATA, j)
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
    GetFigure(DATA.tag.clusters);
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
        CompareShapes(DATA,2);
        return;
    elseif strcmp(plottype,'spkpeaks')
        subplot(1,1,1);
        hold off;
        if 1
        imagesc(DATA.peakpos,'buttondownfcn',{@HitImage, 1})
        set(gcf, 'KeyPressFcn',{@KeyPressed, 1});
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
        GetFigure(DATA.tag.clusters);
        for j = 1:size(DATA.nspks,1)
            m(j,:) = 10000 .* DATA.nspks(j,:)./sum([Expts{j,1}.Trials.dur]);
        end
        imagesc(m,'buttondownfcn',{@HitImage, 1});
        colorbar;
        return;
    elseif strcmp(plottype,'EventRate')
        GetFigure(DATA.tag.clusters);
        for j = 1:size(DATA.nevents,1)
            m(j,:) = 10000 .* DATA.nevents(j,:)./sum([Expts{j,1}.Trials.dur]);
        end
        imagesc(m,'buttondownfcn',{@HitImage, 1});
        colorbar;
        return;
    elseif strcmp(plottype,'TriggerFind')
        GetFigure(DATA.tag.clusters);
        m = max(DATA.mahal(:,:,[1 3]),[],3);
        hold off; 
        imagesc(m,'buttondownfcn',{@HitImage, 4});
        caxis([0 5]);
        hold on;
        [a,b] = find(m > 2 & DATA.dropi < 1 & abs(DATA.peakdiff) < 1);
        for j = 1:length(a);
            DrawBox(a(j),b(j),'color',[0.5 0.5 0.5]);
        end
        vc = prctile(DATA.mspkvar(:),50);
        [a,b] = find(m < 2 & abs(DATA.peakdiff) < 1 & DATA.muampl > 0.9 & DATA.mspkvar > vc);
        for j = 1:length(a);
            DrawBox(a(j),b(j),'color',[0.0 0.0 0.0]);
        end
    elseif strcmp(plottype,'Tagged')
        imagesc(DATA.tagged,'ButtondownFcn',{@HitImage, 1});
    end
    
    if DATA.rebuild || (strcmp(plottype,'xcorr') && ~isfield(DATA,'synci'))
    for j = 1:length(DATA.strings)
        if length(Clusters{j}) > 1
            GetFigure(DATA.tag.clusters);
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
                    GetFigure('Monitor');
                    subplot(1,1,1);
                    plot(ts);
                    drawnow;
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
    
    if strmatch(plottype,{'drop-mahal' 'mahalxc' 'mahaldp' 'dxpc' 'PcGms' 'mahalvar' 'BuildTimes' 'Fit-mahal' 'stability'})
        SetFigure(DATA, DATA.tag.popscatter,'front');
    else
        SetFigure(DATA, DATA.tag.all,'front');
    end
    hold off;
    if strcmp(plottype,'dips')
        imagesc(mahal,'buttondownfcn',{@HitImage, 1});
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
        PlotPopPoints( max(DATA.mahal(:,:,[1 3]),[],3),DATA.dropi, DATA,'colorscheme', colors, plotargs{:});
    elseif strcmp(plottype,'muamp-spkvar')
        PlotPopPoints(DATA.muampl, DATA.mspkvar(:,:,1),DATA,plotargs{:});
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
        GetFigure(DATA.tag.templatesrc);
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
        PlotMahalImage(DATA,'colorbarpos',[0.85 0.05 0.1 0.9]);
        set(gca,'yticklabel',[]);
    elseif strcmp(plottype,'mahal')
        subplot(1,1,1);
        PlotMahalImage(DATA);
        for j = 1:length(Clusters)
            for k = 1:length(Clusters{j})
                if Clusters{j}{k}.auto == 0
                    DrawBox(j,k,'color','r');
                end
            end
        end
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
            GetFigure('Pairs');
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
            GetFigure('Pairs');
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
            GetFigure('Pairs');
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
   end



function [DATA, Clusters] = ReadClusterResults(DATA, Clusters)
if DATA.datatype == 2
    DATA = ReadTemplateResults(DATA, DATA.templatesrc);
    return;
end

lastc.ctime = Clusters{1}{1}.ctime;
for j = 1:length(Clusters)
    for k = 1:length(Clusters{j})
        if length(Clusters{j}{k}.mahal) < 4
            Clusters{j}{k}.mahal(4) =0;
        end
        if ~isfield(Clusters{j}{k},'Evec')
            Clusters{j}{k}.Evec.Eval = NaN;
        end
        if ~isfield(Clusters{j}{k},'probe')
            Clusters{j}{k}.probe = k;
        end
        if isfield(Clusters{j}{k},'gmfit') && isempty(Clusters{j}{k}.gmfit)
            Clusters{j}{k} = rmfield(Clusters{j}{k},'gmfit');
        end
        Clusters{j}{k}.exptid = j;
        C = Clusters{j}{k};
        if C.space(1) == 6
            CheckFitDim(C);
        end
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
        DATA.mahal(j,k,1) = C.mahal(1); %2-D
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
            Clusters{j}{k}.bestspace = [NaN 0];
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
        DATA.dpsum(j,k) = sum(C.dpsum);
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
        DATA.nevents(j,k) = C.nspks;
        DATA.nspks(j,k) = sum(C.clst > 1);
        DATA.dropi(j,k) = C.dropi(3);
        lastc = C;
    end
end
setappdata(DATA.toplevel,'Clusters',Clusters);
   
Expts = getappdata(DATA.toplevel,'Expts');
if length(Expts)
    for j = 1:size(Expts,1)
        for k = 1:size(Expts,2)
            ctimes(j,k) = Expts{j,k}.Header.CreationDate;
        end
    end
    estart = min(ctimes(ctimes > 0));
    for j = 1:size(Expts,1)
        for k = 1:size(Expts,2)
            Expts{j,k}.Header.timeoffset = (Expts{j,k}.Header.CreationDate-estart).*24.*60.*60;
        end
    end
    setappdata(DATA.toplevel,'Expts',Expts);
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
    b = regexprep(b,'M([0-9]*)(.*)','M$1');
    [c,d] = fileparts(a);
    for k = 1:length(DATA.exptid)
        for p = 1:DATA.nprobes
            xs = '';
            if rem(DATA.exptid(k),1) > 0.001
                xs = 'a';
            end
            name = [DATA.name '/Spikes/' b '.p' num2str(p)  't' num2str(floor(DATA.exptid(k))) xs '.mat'];
            if verbose > 0
            fprintf('Reading %s at %s\n',name,datestr(now));
            end
            AllSpikes{k,p} = ReadSpikeFile(name);
            AllSpikes{k,p}.probe = p;
        end
    end
    fprintf('Spike Load took %.2f\n',mytoc(ts));
    setappdata(DATA.toplevel,'AllSpikes',AllSpikes);

function LoadSelectedSpikes(DATA,eid,pid)
    ts = now;
    [a,b] = fileparts(DATA.name);
    [c,d] = fileparts(a);
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    for j = 1:length(eid);
        e = eid(j);
        p = pid(j);
        xs = '';
        if rem(DATA.exptid(e),1) > 0.001
            xs = 'a';
        end
        name = [DATA.name '/Spikes/' b '.p' num2str(p)  't' num2str(floor(DATA.exptid(e))) xs '.mat'];
        fprintf('Reading %s at %s\n',name,datestr(now));
        AllSpikes{e,p} = ReadSpikeFile(name);
        AllSpikes{e,p}.probe = p;
    end
    fprintf('Spike Load took %.2f\n',mytoc(ts));
    setappdata(DATA.toplevel,'AllSpikes',AllSpikes);

function Spikes = ReadSpikeFile(spkfile)
if ~exist(spkfile)
    fprintf('Cannot read %s\n',spkfile);
    Spikes = [];
    return;
end
load(spkfile);
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end
    
function h = QuickSpikes(DATA, pos, varargin)

    if isstruct(pos)
        Spks = pos;
        C = varargin{1};
        e = find(DATA.exptid == C.exptno) ;
        p = C.probe(1);
        pos = [e p];
    else
        e = pos(1);
        p = pos(2);
        AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
        Clusters = getappdata(DATA.toplevel,'Clusters');
        C = Clusters{e}{p};
        Spks = AllSpikes{e,p};
    end
    
    showmean = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'showmean',8)
            showmean = 1;
        end
        j = j+1;
    end
    nevents = length(Spks.times);
    if nevents > 50000
        nspks = 1000;
    else
        nspks = 500;
    end
    istep = max([1 round(nevents/nspks)]);
    ispk = 1:istep:nevents;
    DATA.usegmcid = 0;
    DATA.plotspk.bytrial = 0;
    DATA.voffset = zeros(1,DATA.nprobes);
    yl = minmax(Spks.values(:));
    h = PlotSpikes(DATA, pos, ispk, Spks, C, 'fixy',yl);
    if showmean
        hold on;
        plot(C.MeanSpike.mu(pos(2),:),'-','color','k');
        for j = 2:length(h)
            if h(j) > 0 & ishandle(h)
                color = get(h(j),'color');
                color = 1-color;
                plot(C.MeanSpike.ms(pos(2),:),'-','color',color);
            end
        end
    end
    for j = 2:length(h)
        cmenu = uicontextmenu;
        for k = 1:20
            uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, j-1,  k});
        end
        set(h(j),'UIContextMenu',cmenu);
    end

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
   subplot(1,1,1);
   stopbtn = findobj(F, 'tag','StopSpool');
   for k = 1:length(Spikes)
       p = Spikes{k}.probe;
       ranges(k,:) = minmax(Spikes{k}.values(:));
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
    for j = 1:length(Trials)
        set(F,'name',sprintf('Ex%.0f Trial %d',DATA.exptid(e),Trials(j).Trial));
        for k = 1:length(Spikes)
            p = Spikes{k}.probe;
            chspk = [p-1 p+1];
            chspk = chspk(chspk > 0 & chspk <= DATA.nprobes);
          mysubplot(4,6,p);
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
    
 AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
 Clusters = getappdata(DATA.toplevel,'Clusters');
 
 atid = find(Clusters{eid}{probes(1)}.clst == clnum(1)+1);
 btid = find(Clusters{eid}{probes(2)}.clst == clnum(2)+1);
 if isempty(AllSpikes)
     xct = -(tlim/10000):0.0002:(tlim/10000);
     [xc, b] = xcorrtimes(Clusters{eid}{probes(1)}.times(atid),Clusters{eid}{probes(2)}.times(btid),'times',xct);
     GetFigure(DATA.tag.spikes);
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
    id = find(abs(dt) < tlim)';
    if length(id)
        aid = [aid j * ones(size(id))];
        bid = [bid id];
        dts = [dts dt(id)'];
    end
end
% id = union(aid,bid);
voff = 4;
x = [1:size(AllSpikes{probes(1)}.values,2)]./40;
hold off;
 for j = 1:length(aid)
     %+v toff = 1 after 2, so subtract from 2
     toff = AllSpikes{eid,probes(1)}.times(atid(aid(j)))-AllSpikes{eid,probes(2)}.times(btid(bid(j)));
     toff = toff ./10;
     plot(x,AllSpikes{eid,probes(1)}.values(atid(aid(j)),:));
     hold on;
     plot(x-toff,AllSpikes{eid,probes(2)}.values(btid(bid(j)),:)+voff);
 end
 text(x(end),voff,sprintf('%d',probes(2)),'horizontalalignment','right','color','r');
 text(x(end),0,sprintf('%d',probes(1)),'horizontalalignment','right','color','r');
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


function stopped = SpoolSpikes(DATA, pos, varargin)

    useids = [];
    args = {};
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
    DATA.Expt = Expts{e,p};
    Trials = Expts{e,p}.Trials;
    AllSpikes = getappdata(DATA.toplevel,'AllSpikes');
    Clusters = getappdata(DATA.toplevel,'Clusters');
    C = Clusters{e}{p};
    Spks = AllSpikes{e,p};
    chspk = [p-1 p+1];
    chspk = chspk(chspk > 0 & chspk <= DATA.nprobes);
    for j = 1:length(chspk)
        if ~isempty(AllSpikes{e,chspk(j)})
            Spks(j+1) = AllSpikes{e,chspk(j)};
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
    else
        useids = 1:length(Trials);
    end
    xcl = FindExcludedTrials(DATA,e,p,1,C);
    for j = 1:length(Trials)
        if ~ismember(j, xcl) || length(useids) == 1
        DATA.currenttrial = useids(j);
        ispk = find(C.times > Trials(j).Start(1)./10000 & C.times < Trials(j).End(end)./10000);
        PlotSpikes(DATA, pos, ispk, Spks, C, 'fixy',[-4 4],args{:});
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
    if length(useids) ==1
        SetTrialList(DATA);
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
    elseif strncmpi(varargin{j},'showall',7)
        showall = 1;
    end
    j = j+1;
end

if ~isfield(Spks,'values')
    fprintf('Missing Spikes for Expt %d\n',e);
    return;
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
    for j = 1:length(nc)
        ids{nc(j)} = spkid(find(Spks(1).codes(spkid,1) == nc(j)-1));
    end
    id = ids{2};
end
ispk = pos(2);
for j = 1:length(ids)
    V{j} = Spks(1).values(ids{j},:)';
end

voff = DATA.voffset - DATA.voffset(ispk(1));

    l = size(V{1},1);
    hold off;
    x = [1:l NaN];
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
        h = text(size(V{1},1)-2,voff(c)+0.5,sprintf('%d',c));
        set(h,'fontweight','bold','horizontalalignment','right');
    end
end
for c = 2:length(Spks)
    [ix, ida, idb] = intersect(round(Spks(c).times.*1000),round(Spks(1).times(spkid).*1000));
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
    h = text(size(V{1},1)-2,voff(Spks(c).probe)+0.5,sprintf('%d',Spks(c).probe));
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
title(sprintf('E%d T%d %.2f-%.2f: %d/%d ed%.2f',e,T.Trial, T.Start(1)./10000,T.End(end)./10000, length(id),length(spkid),DATA.Expt.Trials(nt).ed));
else
title(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...
    spkid(1),spkid(end),Spks(1).times(spkid(1)),Spks(1).times(spkid(end)),length(id),length(spkid)));
end
end

