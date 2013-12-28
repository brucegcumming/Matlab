function res = AllVPcs(V, varargin)
%AllVPcs(V, ...)  takes an MxN matrix (electrode by voltage) of continuous
%voltages, extracts segments triggered off one row, and plots PCS.
%  V can also be a filename of a FullV file.  Will only load if new. 
%
%AllVPcs(V,'tchan',c,..    uses Channel c to find trigger points.
%AllVPcs(V,'tchan',c,'reclassify')  uses saved clustering and exaclty
%recaptiulates
%AllVPcs(V,'tchan',c,'reapply')  uses saved clustering parameters, but
%appies them to new data (e.g. changes in trigger, new probe).
%AllVPcs(V,'tchan',c,'reapply', Clusters{p}) uses cluster given
%
%AllVPcs(V, 'tchan',c,'trigdt')  trigggers off dVdt
%AllVPcs(V, 'tchan',c, 'triggertemplate', Clusters{p}) convolves trgger
%    channel with template in cluster, triggers on peaks


%AllVPcs(V, 'tchan',c,

%to build the V file, use MakeProbeIndex to build a list of files
% then u
%            FullV =
% %
% PlotSpikeC(files,5,'probes',probelist,'sumv','submean','makev');
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
DATA.addmean = 0; %New default Jan 2012. Add mean back just in avg.

ispk = 0;
DATA.logfid = -1;
%Save memory if load FullV as int and leave them  as ints.  If this causes
%trouble, setting convertarg to 'convert' forces conversion to doulbe at
%load time
convertarg = {'noconvert'};

nprobepc = 1; % number of probes to include in in pc calculation
DATA.clplot = 0; %1 for density plot
plottype = 1; %plots made on first pass Could probably be 0 now...
plotv = 0;  %plot spikes superimposed on full voltage trace.
clusterprops = [];
spts = [];
gettrigtimes = 0;
smoothsd = 0; %smoothing just for trigger criterion
smoothv = 0; %smoothing of all traces for subsequent processing
subtractadj = 0; %subtract off addjacent channeld
tryall = [1 0 0 0]; %default set of CSD/dvdt combinations to use
dvdts = [0 1 0 0];
csds = [0 0 1 2];
rejectbydiff = 0;
calcdistancematrix = 0;
savedistancematrix = 0;
DATA.idstr = '';
DATA.dvdt = 0;
DATA.plotdvdt = 0;
DATA.tag.Expt = 'ExptFig';
DATA.autocutmode = 'mahal';
DATA.autorefine = 0;
DATA.csd = 0;
DATA.dvdy = 0;
DATA.ncomponents =2;
DATA.name = [];  %directory for current cell
DATA.fullvname = []; %full V filename, current expt
DATA.plotspk.probes = [];
DATA.plotspk.bytrial = 1;
DATA.plotspk.allprobes = 0;
DATA.plotspk.showmean = 0;
DATA.plotspk.showfullv = 0;
DATA.plotspk.nfullvprobes = 1;
DATA.plotspk.includeprepost = 0;
DATA.savetrigger = 0;

DATA.plottype = 1;
DATA.hanldes = [];
DATA.ptsz = [1 1];
DATA.comparecell = [];
DATA.usealltrials = 0;
%1 PCS  2-ADCvals  3-template
%10-ICS
DATA.tmpnorm = 1;
DATA.plot.xcorr = 0;
DATA.version = 1.2;
DATA.plot.expt = 0;
DATA.plot.expttype = 'means';
DATA.plot.exptfit = 0;
DATA.plot.scaledensity = 0;
DATA.hidecluster = 0;
DATA.vsmps = [20 6 15 11 30 20];
DATA.currentcluster = 1;
DATA.autofit.maxthriter = 0; %max number of iterations lowering threshold
DATA.usebmi = 0; %now do evertyhing on GM fit. Calculating old indices wastes a lot of time
DATA.lastcut = [];
DATA.SpaceTypes = {'Pcs' 'VarE'  'RawV' 'Template'};
DATA.elmousept.down = 0;
DATA.elmousept.shape = -1;
DATA.elmousept.handles = [];
DATA.vsep = 4;
DATA.user = GetUserName();
DATA.readlayout = 0;
DATA.loadedClusterDetails = 0;
prefsdir = '/bgc/group/matlab/preferences/AllVPcs';
DATA.layoutfile = [];
DATA.configfile = [];
DATA.plot.isi = 0;
DATA.plot.xyseq = 0;
DATA.plot.xcorrprobes = 0;
DATA.subtractmeanV = 0;
DATA.isicheck = [20 3];
DATA.exptno = 0;
DATA.version = 1.1;
DATA.trigdt = 0;
DATA.maintitle=-1;
DATA.clustericon = -1;
DATA.clustersubdir = [];
DATA.quickcutmode.fit1cut = 0;
DATA.quickcutmode.fit2gauss = 1;
DATA.quickcutmode.fitallcut = 0;
DATA.quickcutmode.quickest = 0;
DATA.quickcutmode.plotspikes = 1;
DATA.quickcutmode.dropi = 0;
DATA.quickcutmode.mean = 0;
DATA.quickcutmode.triggerhist = 1;
DATA.fullvswitchmode.summary = 0;
DATA.iteratefit = 0;
DATA.interactive = 0;
DATA.expttype = 'Default';
DATA.auto.advanceprobe = 0;
DATA.auto.advanceexpt = 0;
DATA.auto.saveref = 0;
DATA.auto.checkxcorr = 0; %check xcorr for cl1/2 if both
DATA.auto.uselastcluster = 0;
DATA.auto.replotcluster = 0;
DATA.auto.checkcluster = 1;
DATA.checkclusters = 1;
DATA.check.dropi = [0 2.5];
DATA.probelabel = [];
DATA.setspkrate = 0;
DATA.refinemode = 'cautious';
DATA.LastClusters = {};
DATA.Comments = [];
DATA.ncelltotry = [5 3];
DATA.verbose  =0;
DATA.errs = {};
DATA.errstates = {};
DATA.saveallspikes = 1;
forcetrigger = 0;

matcheventcounts = 0;
fixerrs = 0;
plotsummary = 0;


DATA.StdTemplate(1,:) = [    -0.1054   -0.1231   -0.1500   -0.2199   -0.4136   -0.7840 ...
      -1.2472   -1.6826   -1.8786   -1.6787   -1.2104   -0.6746   -0.2085    0.1502 ...
          0.4361    0.6637    0.8064    0.8326    0.7798    0.7036    0.6378    0.5895 ...
             0.5435    0.5011    0.4573    0.4182    0.3855    0.3550    0.3323...
             0.3 0.27 0.24 0.21 0.19 0.17 0.15 0.14 0.13 0.12 0.11]; 
DATA.StdTemplate(2,:) = [ -0.036 0.058 0.351 0.804 1.108 0.712 -0.372 -1.346 -1.804 -1.445 0.000 0.700 1.000 0.800 0.650 0.540 0.450 0.400 0.350 0.300 0.250 0.218 0.173 0.139 0.108 0.081 0.064 0.049 0.042 0.035 0.031 0.026 0.022 0.022 0.018 0.014 0.010 0.011 0.012 0.010 ];
DATA.StdTemplate(3,:) = [ -0.036 0.058 0.351 0.804 1.108 0.712 -0.372 -1.346 -1.804 -1.445 0.000 0.700 1.000 0.800 0.650 0.540 0.450 0.400 0.350 0.300 0.250 0.218 0.173 0.139 0.108 0.081 0.064 0.049 0.042 0.035 0.031 0.026 0.022 0.022 0.018 0.014 0.010 0.011 0.012 0.010 ];
DATA.usestdtemplates = 0;
DATA.TemplateUsed = [];
DATA.Template = [];
DATA.RefClusters = {};
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
unsafetosave = 0;
newdata = 0;
DATA.nolog = 0;
saveautocut = 0;
readclusterfromlog = 0;
DATA.trigtimes = {};
DATA.artifacttimes = [];
DATA.savespikes = 0;
DATA.watcharg = {};
DATA.watchplots = 1;
DATA.profligate = 0;
DATA.keepsmooth = 0;
DATA.gui.fontsize(1) = 14;


DATA.usegmcid = 0;
DATA.restricttimerange = [];
DATA.excludetrialids = [];

DATA.colors{1} = [0.5 0.5 0.5];
DATA.colors {2} = [1 0 0];
DATA.colors {3} = [0 1 0];
DATA.colors {4} = [0 0 1];
DATA.colors {5} = [1 0 1];
DATA.colors {6} = [1 1 0];
DATA.colors {7} = [0 1 1];
DATA.colors {8} = [0 1 0];
DATA.colors {9} = [0 1 0];
DATA.preperiod = 0.05;
DATA.postperiod = 0.1;
DATA.usetrials = 1;
DATA.tag.vhist = 'Vhist';
DATA.tag.spikes = 'Spikes';
DATA.tag.allspikes = 'AllSpikes';
DATA.tag.allxy = 'AllXY';
DATA.tag.hist = 'Hist';
DATA.tag.top = 'PCs';
DATA.tag.tmplscore = 'TemplateScores';
DATA.tag.vare = 'VarE';
DATA.tag.meanspike = 'MeanSpike';
DATA.tag.covar = 'Covar';
DATA.tag.dips = 'Dips';
DATA.tag.xcorr = 'Xcorrs';
DATA.tag.oldxy = 'previousXY';
DATA.tag.comparexy = 'CompareXY';
DATA.probeswitchmode = 'reclassify';
DATA.probe = 0; % if no fields, tests cause errors
DATA.tag.allxy = 'AllXY';
DATA.tag.celllist = 'AllCellList';
DATA.tag.fullv = 'FullV';
DATA.tag.comments = 'AllVcomment';
DATA.profiling = 0;
DATA.defaultconfig = [prefsdir '/' gethostname '.' GetUserName '.config'];
DATA.defaultlayout = [prefsdir '/' gethostname '.' GetUserName '.layout.mat'];

DATA.plotrv = 0;

useguidata = 0;
useoldlst = 0;
userefcluster = 0;
checklast = 1;
reclassifyall = 0;
refineall = 0;
saveclusters = 0;
spkrate = 50;
autocutone = 0;
forcecluster = 0;
forceclusterexpt = 0;
maxspksallowed = 600000;
maxspikeset = NaN;
maxrate = NaN;
recluster = 0;
refinecluster = 0;
vt = [];
plotdprimemax = 0;  %old way to find boundaries. Really no good.
bmcrit = 0.21;
verbose = 1;
DATA.cstarttime = now;
oldxyloaded = 0;
loadfromspikes = 0;

spoolspikes = 0;
forcedrive = 'C:/bgc/data';
forcedrive = [];
forcename = [];
fullVname = [];
errs = {};
nerr = 0;
DATA.showerrs = 0;
DATA.showdipvals = 0;

ttn = 1;
tt(ttn).time = now;
tt(ttn).str = 'Start';
ttn = ttn+1;
if length(varargin) & strcmp(varargin{end},'autocutall')
    autocutall = 1;
else
    autocutall = 0;
end


callstring = [];
j = 1;
while j <= length(varargin)  %some varags must be parsed first
    if isnumeric(varargin{j});
        if length(varargin{j}) > 2
            callstring = [callstring ' ' num2str(varargin{j}(1)) ':' num2str(varargin{j}(end))];
        else
            callstring = [callstring ' ' num2str(varargin{j})];
        end
    elseif ischar(varargin{j})
        callstring = [callstring ' ' varargin{j}];
    elseif isstruct(varargin{j})
        if isfield(varargin{j},'probe')
            callstring = [callstring 'Clusters{' num2str(varargin{j}.probe) '}'];
        end
    end
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
    elseif strncmp(varargin{j},'Expt',4)
        j = j+1;
        DATA.Expt = varargin{j};
    elseif strncmpi(varargin{j},'exptname',7)
        j = j+1;
        DATA.exptname = varargin{j};
    elseif strncmpi(varargin{j},'Gridonline',8)
        DATA.DataType = 'GridData 96';
        DATA.expttype = 'gridonline';
    elseif strncmpi(varargin{j},'GridData',5)
        DATA.DataType = 'GridData 96';;
    elseif strncmpi(varargin{j},'name',4)
        j = j+1;
        forcename = varargin{j};
    elseif strncmpi(varargin{j},'highpass',6)
        convertarg = {convertarg{:} varargin{j} varargin{j+1}};
        j = j+1;1;
        DATA.keepsmooth = 1;
    elseif strncmp(varargin{j},'noninteract',7)
        DATA.interactive = -1;
        DATA.checkclusters = 0;
    elseif strncmpi(varargin{j},'toplevel',6)
        j = j+1;
        DATA.toplevel = varargin{j};
    elseif strncmp(varargin{j},'tchan',5)
        j = j+1;
        DATA.probe = varargin{j};
        ispk = varargin{j};
    elseif strncmpi(varargin{j},'toptag',6)
        j = j+1;
        DATA.tag.top = varargin{j};
    elseif strncmp(varargin{j},'usealltrials',7)
        DATA.usealltrials = 1;
    elseif strncmp(varargin{j},'verbose',4)
        verbose = 2;
        DATA.verbose = 2;
    end
    j = j+1;
end

if ischar(V)
    if strcmp(V,'nowarn');
        warning('off','stats:gmdistribution:FailedToConverge');
        warning('off','stats:gmdistribution:MaxIterations');
        return;
    end
    if strcmp(V,'quickop'); %timesaver for working on a new routine.
        it = findobj('type','figure','tag',DATA.tag.top);
        DATA = get(it,'UserData');
        PCCluster(DATA,[],'PCmultiple');
        return;
    end
    if strncmp(V,'tmplspace',8) %set space for manual ND auto cut
        it = findobj('type','figure','tag',DATA.tag.top);
        DATA = get(it,'UserData');
        DATA.tmplspace(2,1:length(varargin{1})) = varargin{1};
        if length(varargin) > 1 && isnumeric(varargin{2})
            DATA.ncelltotry = varargin{2}; %number os GM clusters to trya, and which ond to use tof plotting points
        end
        fprintf('Dimensions are:\n');
        for k = 1:length(DATA.TemplateLabels)
            fprintf('%d: %s\n',k,DATA.TemplateLabels{k});
        end
        if strncmp(V,'tmplspacefit',11)
           PCCluster(DATA,[], 'NCellTemplate');
        else
            set(DATA.toplevel,'UserData',DATA);
        end
        return;
    end
    if strncmp(V,'plottrial',8);
        id = varargin{1};
        it = findobj('type','figure','tag',DATA.tag.top);
        DATA = get(it,'UserData');
        if strncmp(V,'plottrialtime',13)
            DATA.onespiketime = id;
            tid = find([DATA.Expt.Trials.Start] < id);
            if isempty(tid)
                tid = find([DATA.Expt.Trials.End] > id);
                id = tid(1);
            else
                id = tid(end);
            end
        end
        [DATA.currenttrial, DATA.spklst] = PlotTrialSpikes(DATA,id,'showall');
        tic;
        set(DATA.toplevel,'UserData',DATA);
        return;
    end
    if exist(V,'file')
        fullVname = V;
        toplevel = findobj('type','figure','tag',DATA.tag.top);
        if isempty(toplevel)
            X.fullvname = '';
            b = 'NOFILE';
        else
            X = get(toplevel,'UserData');
            b = regexprep(X.fullvname,'\.p[0-9]*FullV','FullV');
        end
        a = regexprep(V,'\.p[0-9]*FullV','FullV');
        if strcmp(V,X.fullvname) && isappdata(toplevel,'Vall')
            DATA = GetGuiState(DATA, toplevel);
            Vall = getappdata(toplevel, 'Vall');
            DATA.Expt = X.Expt;
            DATA.name = X.name;
            DATA.interval = X.interval;
            DATA.fullvname = X.fullvname;
        elseif strcmp(a,b) && isappdata(toplevel,'Vall')
            DATA = GetGuiState(DATA, toplevel);
            if verbose > 1
                fprintf('%s Loading %s\n',IDStr(DATA,V));
            end
            FullV = LoadFullV(V, convertarg{:});
            if verbose > 1
                fprintf('Loaded %s\n',IDStr(DATA,V));
            end
            FullV = SetFullVNames(DATA, FullV);
            DATA.Expt = X.Expt;
            DATA.name = X.name;
            DATA.interval = X.interval;
            DATA.fullvname = fullVname;
            V = FullV;
            clear FullV;
        elseif 1 || ~autocutall %need this for autocutall too
            DATA.fullvname = V;
            tt(ttn).time = now;
            tt(ttn).str = sprintf('Loading %s (No GUI)',fullVname);
            ttn = ttn+1;

            %fi figure is already up, clear existing Full first to reduce memory use
            if length(toplevel) == 1 
                if isappdata(toplevel,'Vall')
                rmappdata(toplevel,'Vall');
                end
                if isappdata(toplevel,'ClusterDetails')
                rmappdata(toplevel,'ClusterDetails');
                end
                DATA.loadedClusterDetails = 0;
            end
            if verbose
                tic;
                fprintf('Loading %s %s',fullVname,datestr(now,'HHMM:ss'));
            end
            FullV = LoadFullV(V, convertarg{:});
            FullV = SetFullVNames(DATA, FullV);
    
            if verbose
                fprintf(' took %.2f (%d bytes/sec)\n',toc,FullV.readrate);
            end
            maxl = size(FullV.V,2)-32;
            chspk = 1:size(FullV.V,1);
            V = FullV;
            clear FullV;
        else
            Vall = V;
        end
        clear X;
    elseif strfind(V,'Spikes') %load AllV from SpikeFiles
        loadfromspikes = 1;
    else 
        it = findobj('tag',V,'type','figure');
        if length(it) == 1
        F = SetFigure(V);
        DATA = get(F,'UserData');
        Vall = getappdata(F,'Vall');
        vt = DATA.t;
        else
            fprintf('Can''t Find %s\n',V);
            return;
        end
    end
end

DATA.loadfromspikes = loadfromspikes;
if loadfromspikes
    DATA.name = V;
    DATA.fullvname = '';
    clear V;
    [DATA, V] = ReadSpikeFiles(DATA, DATA.name);
end

if isstruct(V)
    newdata = 1;
    Vall = V;

    if isfield(Vall,'loadname')
        [a,fname] = fileparts(Vall.name);
        DATA.name = [fileparts(Vall.loadname)];
        if size(Vall.V,1) ==1
            p = GetProbeFromName(Vall.loadname);
            if p > 0 && ispk(1) > 0 && p ~= ispk(1);
                newname = sprintf('%s/Expt%d.p%dFullV.mat',DATA.name,Vall.exptno,ispk(1));
                fprintf('Loading %s\n',newname);
                Vall = LoadFullV(newname,convertarg{:});
            end
        end
    else
        DATA.name = Vall.name;
    end
    
    if isfield(V,'Spikes')
    elseif ~isfield(Vall,'V') || isempty(Vall.V)
        return;
    elseif isinteger(Vall.V) && strcmp(convertarg{1},'convert');
        fprintf('Vall Needs converting -> double\n');
        Vall.V = double(Vall.V) .* Vall.intscale(1)/Vall.intscale(2);
    end
%new data but figure is already up, 
    if DATA.interactive >= 0 || ~isfield(DATA,'Expt')
    DATA.Expt = [];
    end
    it = findobj('type','figure','tag',DATA.tag.top);
    if length(it) == 1  & useguidata
        DATA = get(it,'UserData');
        DATA = ResetDataForNewProbe(DATA);
    elseif length(it) == 1 
        DATA = GetGuiState(DATA, it);
        DATA.toplevel = it;
    end
    if length(forcename)
        DATA.name = forcename;
    elseif length(forcedrive)
        DATA.name = regexprep(DATA.name,'[A-Z]:/Spike2/data',forcedrive);
        DATA.name = regexprep(DATA.name,'[A-Z]:/smr',forcedrive);
    end
    if isfield(Vall,'exptno')
        DATA.exptno = Vall.exptno;
        if isfield(DATA,'toplevel') && isfigure(DATA.toplevel) && isappdata(DATA.toplevel,'Expts')
            if verbose > 1
                fprintf('%s Using Expt appdata\n',IDStr(DATA));
            end
            Expts = getappdata(DATA.toplevel,'Expts');
        else
            Expts = {};
        end
        if DATA.interactive < 0 && isfield(DATA.Expt,'Trials')
            DATA.Expt.Header.trialdur = sum([DATA.Expt.Trials.dur]);
            fprintf('%s Expt already loaded\n',IDStr(DATA));
        elseif isfield(DATA.Expt,'exptno') && DATA.Expt.exptno == Vall.exptno && ~isempty(Expts)
            Expts = Expts;
        elseif length(Expts) >= Vall.exptno
            DATA.Expt = Expts{Vall.exptno};
        elseif isfield(Vall,'matfile')
            
            onefile = strrep(Vall.matfile,'.mat',['Expt' num2str(Vall.exptno) '.mat']);
            if exist(onefile,'file')
                DATA.Expt = LoadExptA(DATA,onefile,0);
            else
                DATA.Expt = LoadExptA(DATA,Vall.matfile,Vall.exptno);
            end
            if isfield(DATA.Expt.Stimvals,'po') && DATA.Expt.Stimvals.po > DATA.postperiod
                DATA.postperiod = DATA.Expt.Stimvals.po;
                DATA.Expt.Header.postperiod = DATA.postperiod .*10000;
            end
            if isfield(DATA.Expt.Stimvals,'pr') && DATA.Expt.Stimvals.pr > DATA.preperiod
                DATA.preperiod = DATA.Expt.Stimvals.pr;
                DATA.Expt.Header.preperiod = DATA.preperiod .*10000;
            end
            DATA.matfile = Vall.matfile;
        else
            [DATA.Expt, DATA.matfile] = LoadExpt(DATA,Vall.exptno);
        end
    if verbose > 1
        fprintf('%s Expt loaded\n',IDStr(DATA));
    end
        
        if isempty(DATA.Expt)
            DATA.cluster.exptreadmethod = -1;
%            DATA.plotspk.bytrial = 0;
        elseif isfield(DATA.Expt.Header,'ReadMethod')
            DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;
        else
            DATA.cluster.exptreadmethod = 0;
        end
        if isfield(DATA.Expt,'Header') && isfield(DATA.Expt.Header,'DataType')
            DATA.DataType = DATA.Expt.Header.DataType;
        else
            DATA.DataType = 'Default';
        end
        DATA.Expt.exptno = Vall.exptno;
        SetTrialList(DATA);
    end
    if isfield(Vall,'firstblk')
        DATA.Expt.blk = Vall.firstblk;
        DATA.Expt.exptno = DATA.Expt.exptno+0.1;
    else
        DATA.Expt.blk = 0;
    end
    if isfield(Vall,'start')
    DATA.starttime = Vall.start;
    end
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
    Vall = getappdata(V,'Vall');
    vt = DATA.t;

elseif ~isfield(DATA,'name')
    Vall.V = V;
    newdata = 1;
end

if ~isdir(DATA.name)
    DATA.datadir = fileparts(DATA.name);
else
    DATA.datadir = DATA.name;
end

Array = GetArrayConfig(DATA.datadir);
if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
    DataClusters = getappdata(DATA.toplevel,'Clusters');
else
    DataClusters = {};
end


if exist('Vall','var') && ~ischar(Vall) && ~isempty(Vall);
    if isfield(Vall,'Spikes')
        DATA.nprobes = 1;
        DATA.allnprobes = length(Array.X);
        chspk = 1:DATA.allnprobes;
    else
        maxl = size(Vall.V,2)-32;
        chspk = 1:size(Vall.V,1);
        DATA.nprobes = size(Vall.V,1);
        DATA.exptno = Vall.exptno;
        DATA.samplerate = 1./Vall.samper;
    end
    if isfield(Vall,'NSlabel')
        DATA.probelabel = deblank(Vall.NSlabel);
    end
    DATA.allnprobes = DATA.nprobes;
    if strncmp(DATA.DataType,'GridData',6)
        DATA.allnprobes = sscanf(DATA.DataType,'GridData %d');
    end
end

if verbose > 1
    fprintf('%s Start args\n',IDStr(DATA));
end


passonargs = {};

j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'adcplot',3)
        DATA.plottype = 2;  
    elseif strncmp(varargin{j},'addmean',7)
        DATA.addmean = 1;
    elseif strncmp(varargin{j},'allspikes',5)
        DATA.usetrials = 0;
    elseif strncmp(varargin{j},'autocutall',10)
        autocutall = 1;
    elseif strncmp(varargin{j},'quickautocutall',13)
        autocutall = 2;
        DATA.autocutmode = 'quick';
    elseif strncmp(varargin{j},'quickautocut',9)
        autocutone = 1;
        DATA.autocutmode = 'quick';
    elseif strncmp(varargin{j},'autocut',7)
        autocutone = 1;
    elseif strncmp(varargin{j},'checklast',8)
        checklast = 2;
        passonargs = {passonargs{:} varargin{j}};
    elseif strncmp(varargin{j},'cutmode',7)
        j = j+1;
        DATA.autocutmode = varargin{j};
    elseif strncmp(varargin{j},'distancematrix',7)
        calcdistancematrix = 6;
        savedistancematrix = 1;
        passonargs = {passonargs{:} varargin{j}};
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
        passonargs = {passonargs{:} varargin{j-1} varargin{j}};

    elseif strncmp(varargin{j},'fontsize',6)
        j = j+1;
        DATA.gui.fontsize = varargin{j};
    elseif strncmp(varargin{j},'forcetrigger',5)
        j = j+1;
        forcetrigger = varargin{j};
    elseif strncmp(varargin{j},'fixerrs',5)
        fixerrs = 1;
    elseif strncmpi(varargin{j},'Gridonline',8)
        DATA.DataType = 'GridData 96';
        DATA.expttype = 'gridonline';
    elseif strncmpi(varargin{j},'GridData',5)
        DATA.DataType = 'GridData 96';;
    elseif strncmp(varargin{j},'ignoretrials',10)
        DATA.plotspk.bytrial = 0;
    elseif strncmp(varargin{j},'config',5)
        DATA.readlayout = 2;
        if length(varargin) > j && ischar(varargin{j+1})
            j = j+1;
            if regexp(varargin{j},'^[A-Z]:') | strfind(varargin{j},'/') %real path
                DATA.configfile = varargin{j};
            else
                DATA.configfile = [prefsdir '/' varargin{j} '.config'];
            end
        end
        DATA = ReadConfig(DATA, DATA.configfile);
    elseif strncmp(varargin{j},'layout',5)
        DATA.readlayout = 2;
         if strncmp(varargin{j},'layoutapply',8)
             DATA.readlayout = 2;
         end
        if length(varargin) > j && ischar(varargin{j+1})
            j = j+1;
            if regexp(varargin{j},'^[A-Z]:') | strfind(varargin{j},'/') %real path
                DATA.layoutfile = varargin{j};
            else
                DATA.layoutfile = [prefsdir '/' varargin{j} '.layout.mat'];
            end
        end
    elseif strncmp(varargin{j},'newexpt',7)
        j = j+1;
        eid = varargin{j};
        if length(varargin) > j
            args = varargin(j+1:end);
        else
            args = {};
        end
        set(DATA.toplevel,'name',['Switching to Expt', num2str(eid)]); 
        drawnow;
        p = ProbeNumber(DATA);
        if DATA.loadfromspikes
            name = [DATA.datadir '/Expt' num2str(eid) 'Spikes'];
        elseif strncmp(DATA.DataType,'GridData',8)
            name = [DATA.name '/Expt' num2str(eid) '.p' num2str(p) 'FullV.mat'];
        else
            name = [DATA.name '/Expt' num2str(eid) 'FullV.mat'];
        end
        if isappdata(DATA.toplevel,'Vall')
            rmappdata(DATA.toplevel,'Vall');
        end
        if isappdata(DATA.toplevel,'ClusterDetails')
            rmappdata(DATA.toplevel,'ClusterDetails');
        end
        clear Vall;
        DATA.trigtimes = {};
        if DATA.loadfromspikes
            AllVPcs(name, 'tchan', p, args{:});
            SetFigureName(DATA.toplevel,DATA.tag.top);
            return;
        end
        if ~exist(name) 
            fprintf('%s Does Not Exist\n',name);
            return;
        end
        if DATA.auto.uselastcluster
            args = {args{:} DATA.LastClusters{p}};
        end
        FullV = LoadFullV(name, convertarg{:});
        FullV = SetFullVNames(DATA, FullV);
        fprintf('Loading %s took  %.2f (%.2f disk) %.1fMb\n',name,FullV.loadtime,FullV.initialloadtime,FullV.size);
        AllVPcs(FullV, 'tchan', p, args{:});
        clear FullV;
        SetFigureName(DATA.toplevel,DATA.tag.top); 
        return;
    elseif strncmp(varargin{j},'nomean',5)
        DATA.addmean = 0;
    elseif strncmp(varargin{j},'noninteract',7)
        DATA.interactive = -1;
        %figure('visible','off');
        DATA.watchplots = 0;
        passonargs = {passonargs{:} varargin{j}};
    elseif strncmp(varargin{j},'nocheck',5)
        checklast = 0;
        DATA.checkclusters = 0;
    elseif strncmp(varargin{j},'nolog',5)
        DATA.nolog = 1;
    elseif strncmp(varargin{j},'ncomponents',6)
        j = j+1;
        DATA.ncomponents = varargin{j};
    elseif strncmp(varargin{j},'ptsz',4)
        j = j+1;
        DATA.ptsz = varargin{j};
        if length(DATA.ptsz) ==1
            DATA.ptsz(2) = DATA.ptsz(1);
        end
    elseif strncmp(varargin{j},'refineall',9) %applies cluster exactly, so should use same data
        if length(varargin) > j && iscell(varargin{j+1})
            j = j+1;
            UseClusters = varargin{j};
        else
            UseClusters = {};
        end
       DATA.autorefine = 3;
       refineall = 1;
    elseif strncmp(varargin{j},'refclusters',9) 
%        j = j+1;
%        rfile = varargin{j};
%default is to use refcluster only on files without previous cut, or only
%an autocut.  Use 'refclusterforce' to force use of the refcluster.
        rfile = [DATA.name '/RefClusters.mat'];
        if exist(rfile,'file')
            load(rfile);
            DATA.RefClusters = CondenseClusters(Clusters,0);
            passonargs = {passonargs{:} varargin{j}};
            if strncmp(varargin{j},'refclusterforce',11) 
                userefcluster = 1;
            end
        end
        
    elseif strncmp(varargin{j},'refinecrit',8)
        j = j+1;
        DATA.autorefine = varargin{j};
    elseif strncmp(varargin{j},'refinemode',8)
        passonargs = {passonargs{:} varargin{j} varargin{j+1}};
        j = j+1;
        DATA.refinemode = varargin{j};
    elseif strncmp(varargin{j},'refine',6)
        refinecluster = 1;
        DATA.autorefine = 3; %only refine in origninal scluster dp > 3
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            DATA.autorefine = varargin{j};
        elseif forcecluster ==1
            DATA.autorefine = 0.1; %if only named 1 cluster, force refine
        end
    elseif strncmp(varargin{j},'submean',6)
        j = j+1;
        DATA.subtractmeanV = 1;
        meanV = varargin{j};
    elseif strncmp(varargin{j},'subdir',6)
        j = j+1;
        DATA.clustersubdir = varargin{j};
    elseif strncmp(varargin{j},'saveautocut',10)
        saveautocut = 1;
    elseif strncmp(varargin{j},'savespikes',6)
        saveclusters = 1;
        if strncmp(varargin{j},'savespikesifsafe',16)
            saveclusters = 2;
        end
        if strncmp(varargin{j},'savespikestoman',14)
            saveautocut = 2;
        end
        DATA.savespikes = 1;
        passonargs = {passonargs{:} varargin{j}};
        if j ==1 && saveautocut ~= 2 %savespikes before anything else = save current GUI spikes
            it = findobj('type','figure','tag',DATA.tag.top);
            if length(it) == 1
                DATA = get(it,'UserData');
            end
            if DATA.cluster.auto == 1 && DATA.recluster == 0
                outname = ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);
            else
                outname = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
            end
            DATA =  SaveClusters(DATA, outname);
            SaveSpikes(DATA, DATA.savespkid);
           return;
        end
    elseif strncmp(varargin{j},'savetrigger',6)
        DATA.savetrigger = 1;
    elseif strncmp(varargin{j},'switchmode',6)
        j = j+1;
        DATA.probeswitchmode = varargin{j};
    elseif strncmp(varargin{j},'subtractadj',9)
        subtractadj = 1;
    elseif strncmp(varargin{j},'summary',6)
        gettrigtimes = 2;     
        plotsummary = 1;
        DATA.probeswitchmode = 'usecluster';
    elseif strncmp(varargin{j},'density',5)
        DATA.clplot = 1;
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
        tryall = 0;
    elseif strncmpi(varargin{j},'accel',4)
        DATA.dvdt = 3;
        tryall = 0;
    elseif strncmpi(varargin{j},'clear',4)
        muscale = 0.8;
    elseif strncmpi(varargin{j},'csda',4)
        DATA.csd = 2;
        tryall = 0;

    elseif strncmpi(varargin{j},'csd',3)
        DATA.csd = 1;
        tryall = 0;

    elseif strncmpi(varargin{j},'logfid',6)
        j = j+1;
        DATA.logfid = varargin{j};
    elseif strncmpi(varargin{j},'iteratefit',8)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            DATA.iteratefit = varargin{j};
        else
            DATA.iteratefit = 1;
        end
    elseif strncmpi(varargin{j},'mine',4)
        j = j+1;
        minenergy = varargin{j};
    elseif strncmpi(varargin{j},'minvar',6)
        j = j+1;
        minvar = varargin{j};
    elseif strncmpi(varargin{j},'newdata',7)
        newdata = 1;
    elseif strncmpi(varargin{j},'ndive',5)
        j = j+1;
        DATA.autofit.maxthriter = varargin{j};
    elseif strncmpi(varargin{j},'showwarn',6)
        DATA.showerrs = 1;
    elseif strncmpi(varargin{j},'plotrv',6)
        DATA.plotrv = 1;
    elseif strncmpi(varargin{j},'postperiod',6)
        j = j+1;
        DATA.postperiod = varargin{j};
    elseif strncmpi(varargin{j},'preperiod',6)
        j = j+1;
        DATA.preperiod = varargin{j};
    elseif strncmp(varargin{j},'tagsuff',5)
        j = j+1;
        DATA.tag.spikes = [DATA.tag.spikes  varargin{j}];
    elseif strncmp(varargin{j},'testclassify',6)
        tic;
        E = BoundaryFromCluster([],DATA.cluster,1);
        E.boundarytype = 0;
        profile on;
        PCCluster(DATA,E,1);
        profile  viewer;
%        ClassifySpikes(DATA,DATA.cluster,'quick');
        toc;
        return;
    elseif strncmp(varargin{j},'tryall',5)
        tryall = ones(size(csds));
    elseif strncmp(varargin{j},'tchan',5)
        if strncmp(varargin{j},'tchannew',7)
            newdata = 2;
            DATA = ResetDataForNewProbe(DATA);
        else
            newdata = 1;
        end
        j = j+1;
        if strcmp(varargin{j},'all')
            ispk = 1:DATA.nprobes;
        else
            ispk = varargin{j};
        end
        if length(ispk) > 1
            addch = 1;
            if strncmp(varargin{j-1},'tchanuse',8)
                addch = 2;
            end
        end
        if DATA.loadfromspikes && newdata == 2
            DATA.probe = ispk;
            [DATA,Vall] = ReadSpikeFiles(DATA, DATA.name);
        end
    elseif strncmp(varargin{j},'ichan',5)
        j = j+1;
        ispk = varargin{j};
        DATA.TemplateLabels = TemplateLabels(DATA,0);
    elseif strncmp(varargin{j},'matchcounts',8)
        matcheventcounts = 1;
        passonargs = {passonargs{:} varargin{j}};
    elseif strncmp(varargin{j},'maxrate',6)
        j = j+1;
        maxrate = varargin{j};
    elseif strncmp(varargin{j},'maxspikes',6)
        j = j+1;
        maxspikeset = varargin{j};
    elseif strncmp(varargin{j},'muscale',6)
        j = j+1;
        muscale = varargin{j};
    elseif strncmp(varargin{j},'nprobepc',6)
        j = j+1;
        nprobepc = varargin{j};
    elseif strncmp(varargin{j},'name',4)
        j = j+1;
        DATA.name = varargin{j};
        passonargs = {passonargs{:} varargin{j-1} varargin{j}};
    elseif strncmp(varargin{j},'noplot',4)
        plottype = 0;
    elseif strncmp(varargin{j},'nspk',4)
        j = j+1;
        setnspk = varargin{j};
    elseif strncmp(varargin{j},'oldcl',4)
        oldcluster = 1;
        if  ~isempty(DataClusters{DATA.probe(1)})
            C = DataClusters{DATA.probe(1)};
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
        nprobepc = -1;
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
            PlotFullV(DATA, [0 100]);
            return;
        end
    elseif strncmp(varargin{j},'recut',5) %recut operates on whats in DATA
        DATA.cluster = DataClusters{DATA.probe(1)};
        DATA.csd = DATA.cluster.csd;
        DATA.dvdt = DATA.cluster.dvdt;
        DATA = ReClassify(DATA,'newbound');
        set(DATA.toplevel,'UserData',DATA);
        return;
    elseif strncmp(varargin{j},'quantifyall',11)
        reclassifyall = 2;
    elseif strncmp(varargin{j},'reclassifyall',13)
        reclassifyall = 1;
    elseif strncmp(varargin{j},'reclassify',6) %reclassify applies cluster space to new data (thr, crit, etc differen)
        recluster= 2;
        if length(varargin) > j && isfield(varargin{j+1},'dropi')
            j = j+1;
            DATA.cluster = varargin{j};
            forcecluster = 1;
        end
    elseif strncmp(varargin{j},'readfromlog',8) %reclassify applies cluster space to new data (thr, crit, etc differen)
        readclusterfromlog = 1;
    elseif strncmp(varargin{j},'uselst',6) %reclassify applies cluster space to new data (thr, crit, etc differen)
        useoldlst = 1;
    elseif strncmp(varargin{j},'usecluster',6) %Uses cluster ids, does not reapply boudnary = quicker
        recluster = 4;
    elseif strncmp(varargin{j},'keepcluster',6) %applies cluster exactly, using current DATA.cluster
        recluster= 1;
        forcecluster = 2;
    elseif strncmp(varargin{j},'reapply',6) %applies cluster exactly, to whatever data is here
        recluster= 1;
        if length(varargin) > j && iscell(varargin{j+1})&& isfield(varargin{j+1}{ispk(1)},'dropi')
            j = j+1;
            DATA.cluster = rmfields(varargin{j}{ispk(1)},'quick');
            forcecluster = 1;
            forceclusterexpt = DATA.cluster.exptno;
        elseif length(varargin) > j && isfield(varargin{j+1},'dropi')
            j = j+1;
            DATA.cluster = varargin{j};
            forcecluster = 1;
            forceclusterexpt = DATA.cluster.exptno;
            if DATA.cluster.probe(1) ~= ispk(1)
                pshift = ispk(1)-DATA.cluster.probe;
                DATA.cluster.templateshift = pshift;
                DATA.cluster.MeanSpike.ms =  circshift(DATA.cluster.MeanSpike.ms,pshift);
                DATA.cluster.TemplateUsed =  circshift(DATA.cluster.TemplateUsed,pshift);
                DATA.cluster.probe = ispk(1);
                p = size(DATA.cluster.MeanSpike.vdprime,1);
                if pshift > 0 && p > 2
                    DATA.cluster.MeanSpike.vdprime(p+pshift-2:p+pshift,:) = DATA.cluster.MeanSpike.vdprime(p-2:p,:);
                end
                for k = 1:length(DATA.cluster.next)
                    if isfield(DATA.cluster.next{k},'MeanSpiike')
                    DATA.cluster.next{k}.MeanSpike.ms =  circshift(DATA.cluster.next{k}.MeanSpike.ms,pshift);
                    DATA.cluster.next{k}.TemplateUsed =  circshift(DATA.cluster.next{k}.TemplateUsed,pshift);
                    DATA.cluster.next{k}.probe = ispk(1);
                    p = size(DATA.cluster.next{k}.MeanSpike.vdprime,1);
                    if pshift > 0 && p > 2
                        DATA.cluster.next{k}.MeanSpike.vdprime(p+pshift-2:p+pshift,:) = DATA.cluster.next{k}.MeanSpike.vdprime(p-2:p,:);
                    end
                    end
                end
            end
            DATA.cluster = rmfields(DATA.cluster,'quick'); %if previous cut was quick, want to start over
        end
    elseif strncmp(varargin{j},'rejectbydiff',8)
        j = j+1
        rejectbydiff = varargin{j};
    elseif strncmp(varargin{j},'spool',4)
        spoolspikes = 1;
    elseif strncmp(varargin{j},'spts',4)
        j = j+1;
        spts = varargin{j};
    elseif strncmp(varargin{j},'vsmps',4)
        j = j+1;
        DATA.vsmps = varargin{j};
    elseif strncmpi(varargin{j},'trigtimes',8)
        gettrigtimes = 1;
    elseif strncmp(varargin{j},'vsmooth',4)
        j = j+1;
        smoothv = varargin{j};
    elseif strncmp(varargin{j},'smooth',4)
        j = j+1;
        smoothsd = varargin{j};
    elseif strncmp(varargin{j},'spkrate',4)
        j = j+1;
        spkrate = varargin{j};
        th = NaN;
    elseif strncmp(varargin{j},'th+',3)
        thsign = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            th = varargin{j};
        end
    elseif strncmp(varargin{j},'th-',3)
        thsign = 0;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            th = varargin{j};
        end
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
    elseif strncmp(varargin{j},'plotold',7)
        DATA.clst = DATA.oldclst;
        ReplotPCs(DATA,[]);
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmp(varargin{j},'plottemplate',8)
        TemplatePlot(DATA);
        DATA = get(DATA.toplevel,'UserData');
    elseif strncmp(varargin{j},'dtthreshold',5)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            th = varargin{j};
        end
         if strncmp(varargin{j},'dtthr2',6)
             DATA.trigdt = 2;
         elseif strncmp(varargin{j},'dtthr3',6)
             DATA.trigdt = 3;
         else
            DATA.trigdt = 1;
         end
    elseif strncmp(varargin{j},'triggertemplate',3)
        if strncmp(varargin{j},'triggertemplatedt',17)
            DATA.trigdt = 5;
        else
            DATA.trigdt = 4;
        end
        j = j+1;
        if length(varargin{j}) == 1
            if isfield(varargin{j},'MeanSpike')
                DATA.TriggerTemplate = varargin{j}.MeanSpike.ms;
            else
                DATA.TriggerTemplate = DATA.StdTemplate(1,:);
            end
        else
            DATA.TriggerTemplate = varargin{j};
        end
    elseif strncmp(varargin{j},'threshold',3)
        if isnumeric(varargin{j+1})
            j = j+1;
            th = varargin{j};
        else
%            fprintf('ERRROR - Input argument ''thr'' needs a value\n');
        end
        DATA.trigdt = 0;
    elseif strncmp(varargin{j},'verbose',4)
        verbose = 2;
        passonargs = {passonargs{:} varargin{j}};
    elseif strncmp(varargin{j},'vpts',3)
        DATA.plottype =2;
    elseif strncmp(varargin{j},'watch',3)
        DATA.watcharg{1} = 'front';
        DATA.watchplots = 2;
        passonargs = {passonargs{:} varargin{j}};
    elseif strncmp(varargin{j},'winfront',4)
        FiguresToFront(DATA.tag);
    elseif strncmp(varargin{j},'nowatch',5)
        DATA.watcharg = {};
        DATA.watchplots = 0;
    elseif strncmp(varargin{j},'xcorr',3)
        DATA.plot.xcorr = 1;
    end
    j = j+1;
end


if (ispk == 0 & reclassifyall) | strcmp(DATA.probe,'all')
    DATA.probe = 1:DATA.nprobes;
    ispk = DATA.probe;
end
if DATA.profiling
    profile on;
end
tt(ttn).time = now;
tt(ttn).strs = 'After varargs';
ttn= ttn+1;


if isfield(Vall,'V')
    if isfield(Vall,'t')  && abs(length(Vall.t) - size(Vall.V,2)) > 1
        PrintMsg(DATA.logfid,sprintf('Error! Vall length (%d) does not match t (%d)',size(Vall.V,2),length(Vall.t)));
        res.err ='time length mismatch';
        return;
    end
end

if verbose > 1
    fprintf('%s Args processed at %s\n',IDStr(DATA),datestr(now,'HH:MM:SS'));
end
if forcecluster == 1
    DATA.cluster.src.probe = DATA.cluster.probe;
    DATA.cluster.src.exptno = DATA.cluster.exptno;
    DATA.cluster.src.ctime = DATA.cluster.ctime;
elseif forcecluster == 2
    X = get(DATA.toplevel,'UserData');
    DATA.cluster = X.cluster;
    DATA.cluster.auto = 0;
    clear X;
end
if isempty(DATA.fullvname) && isfield(Vall,'loadname')
    DATA.fullvname = Vall.loadname;
end
%want this after processing varargin so that can do autocutall
if autocutall
    if length(forcedrive) > 2
        Vall.name = regexprep(Vall.name,'[A-Z]:/Spike2/data',forcedrive);
    end
F = SetFigure(DATA.tag.top, DATA);
DATA.toplevel = F;
set(DATA.toplevel,'UserData',DATA)
usearg = [];
for j = 1:length(varargin)
    if sum(strcmp(varargin{j},{'autocutall' 'quickautocutall'})) == 0 
        usearg(j) = true;
    end
    j = j+1;
end
res = AutoCutAll(ispk,  F, Vall, DATA, varargin(find(usearg)));
res.toplevel = F;
res.logfid = DATA.logfid;
return;
end
if reclassifyall == 2
        DataClusters = LoadCluster(DATA.name,DATA.exptno);
        for j = 1:length(ispk)
            if length(DataClusters) >= ispk(j)
            C = DataClusters{ispk(j)};
            needc = 0;
            if DataClusters{ispk(j)}.quick || DataClusters{ispk(j)}.trigdt == 4 || ...
                    DataClusters{ispk(j)}.manual == 2 ||  DataClusters{ispk(j)}.dropi(3) == 0
                needc = 1;
            end
            if DataClusters{ispk(j)}.auto == 1 && length(DATA.RefClusters) >= ispk(j);
                needc = 1;
            end
            for k = 1:length(DataClusters{ispk(j)}.next)
                if ~isfield(DataClusters{ispk(j)}.next{k},'space')
                    ; %do nothing
                elseif (isfield(DataClusters{ispk(j)}.next{k},'quick') && DataClusters{ispk(j)}.next{k}.quick) ...
                        || (isfield(DataClusters{ispk(j)}.next{k},'manual') && DataClusters{ispk(j)}.next{k}.manual == 2)...
                    || DataClusters{ispk(j)}.next{k}.dropi(3) == 0
                    needc = 1;
                end
            end
            if needc == 0
                DATA.probe = ispk(j);
                d = dir(SpkFilename(DATA));
                if isempty(d) || d.datenum < C.savetime(1)-0.01;
                    fprintf('Spike File %s is older than cluster\n',d.name);
                    needc = 1;
                end
            end
            if DataClusters{ispk(j)}.trigdt == 4
                DataClusters{ispk(j)}.trigdt = 0;
            end
            else
                needc = 1;
            end
            if needc
                fprintf('Reclassifying E%d P%d\n',DATA.exptno,ispk(j));
                try
                res{j} = AllVPcs(Vall, 'tchan', ispk(j), passonargs{:},'reclassify');
                catch ME
                    res{j}.err = lasterr;
                    t = getCurrentTask();
                    res{j}.workerid = t.ID;
                    res{j}.errstate = ME;
                end
            end
            res{j}.needed = needc;
            res{j}.exptno = DATA.exptno;
            res{j}.probes = ispk(j);
        end
        return;
end
if reclassifyall || refineall
    if userefcluster && isfield(DATA,'RefClusters')
        UseClusters = DATA.RefClusters;
    end
    if refineall && isempty(UseClusters)
        UseClusters = LoadCluster(DATA.name,DATA.exptno);
    end
    for j = 1:length(ispk)
        try
        if refineall
%in this loop, want to force p to == ispk(j) regardless of Vall
%if its differente, want to force loading of new V when AllVPcs is called
            if size(Vall.V,1) ==1 && ispk(j) ~= Vall.chspk
                p = ispk(j);
            else
                p = ispk(j);
            end
            if isfield(UseClusters{p},'space')
            res{j} = AllVPcs(Vall, 'tchan', ispk(j), passonargs{:},'reapply',UseClusters{p},'refine', DATA.autorefine);
            res{j}.refinemode = DATA.refinemode;
            end
        elseif reclassifyall == 2
            if Clusters{ispk(j)}.quick
                res{j} = AllVPcs(Vall, 'tchan', ispk(j), passonargs{:},'reclassify');
            end
                
        else
            res{j} = AllVPcs(Vall, 'tchan', ispk(j), passonargs{:},'reclassify');
        end
        if fixerrs && res{j}.err == 1 && res{j}.auto == 1
            AllVPcs(Vall,'savespikes');
        end
        catch ME
                    res{j}.err = ME.message;
                    res{j}.errstate = ME;
        end
        res{j}.exptno = Vall.exptno;
        res{j}.probes = ispk(j);
%Gets existing handle from Gui now - no need to close        
        if 0 && isfield(res{j},'logfid') && res{j}.logfid > 2
            try
                fclose(res{j}.logfid);
            catch
                fprintf('Couldnt close handle %d\n',res{j}.logfid);
            end
        end
        drawnow;
    end
    if savedistancematrix
        clname = ClusterFile(DATA.name, DATA.Expt, 'DistanceMatrix');
        for j = 1:length(res)
            if isfield(res{j},'DistanceMatrix')
                DistanceMatrix{j} = rmfields(res{j}.DistanceMatrix,'cid');
            else
                DistanceMatrix{j} = {};
                fprintf('E%dP%d Missing Distance Matrix\n',DATA.exptno,j);
            end
        end
        fprintf('Saving Distance Matrices to %s\n',clname);
        try
        save(clname,'DistanceMatrix');
        catch ME
            fprintf('Error writing %s\n',clname);
            res{1}.errstate = ME;
        end
    end
    return;
end

if ~isempty(vt) && ispk(1) == 0 %Called with current figure, just to set a variable
    set(DATA.toplevel,'UserData',DATA);
    return;
end

if DATA.showerrs
        warning('on','stats:gmdistribution:FailedToConverge');
        warning('on','stats:gmdistribution:MaxIterations');
else
        warning('off','stats:gmdistribution:FailedToConverge');
        warning('off','stats:gmdistribution:MaxIterations');
end

if isfield(Vall,'Spikes')
    if isempty(DataClusters)
        [DataClusters, DATA.FullVData] = LoadDataClusters(DATA);
    end
    DATA.duration = DataClusters{1}.duration;
    DATA.probelist = 1:DATA.allnprobes;
    if strncmp(DATA.DataType,'Grid',4)
        DATA.probelist = DATA.probe(1);
        DATA.probe = 1;
        ispk = 1;
    end
else
    if isempty(Vall) && isempty(DATA.fullvname)
        DATA.fullvname = regexprep(DATA.name,'(Expt[0-9]*).*','$1FullV.mat');
        Vall = LoadFullV(DATA.fullvname,convertarg{:});
    end
    [DATA, Vall, ispk, newdata] = SetupVall(DATA, Vall, ispk, newdata);
    vt = Vall.t;
end

if length(ispk) > 4 && ~ autocutall
    SetFigure(DATA.tag.top,DATA);
    args = {};
    clname = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    if exist(clname,'file')
        load(clname);
        if calcclscores
            args = {args{:},'Clusterscores',Clusters};
        end
    end
    [nr,nc] = Nsubplots(length(ispk));
    for j = 1:length(ispk)
        res{j} = AllVPcs(Vall, varargin{:},'tchan',ispk(j),'noplot',args{:});
        if thsign == 2
            res{j}.pres = AllVPcs(Vall, varargin{:},'tchan',ispk(j),'th+','noplot',args{:});
        end
        SetFigure(DATA.tag.top,DATA);
        subplot(nr,nc,j);
        id = [];
        PlotPCs(res{j}.pcs,1,2,DATA.clplot,res{j}.id,DATA.colors);
        if DATA.showdipvals
        title(sprintf('%.1f dp%.1f,%.1f',max(res{j}.dipvals)*100,res{j}.dp,res{j}.edp));
        end
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
    if isfield(res,'dp')
    SetFigure('Dprimes');
    hold off;
    plot([x.dp]);
    hold on;
    plot([x.edp],'r');
    if isfield(x,'pres')
        p = cat(1,x.pres);
        plot([x.dp],'o-');
        plot([x.edp],'ro-');
    end
    end
    return;
end


[F, isnew] = SetFigure(DATA.tag.top,DATA);
DATA.toplevel = F;
if isnew
if isempty(DATA.configfile)
    DATA = ReadConfig(DATA,DATA.defaultconfig,'print');
end
if isempty(DATA.layoutfile)
    DATA.layoutfile = DATA.defaultlayout;
    ApplyLayout(DATA,'print');
end
end
res.toplevel = F;

if newdata %things to reset if its a new file, OR just a new probe
        DATA.probe = ispk;
        SetTrialList(DATA);
        if isfield(DATA,'toplevel')
        for f = {'AllCSD' 'AllV'} %can clog up memory
        if isappdata(DATA.toplevel,f{1})
            rmappdata(DATA.toplevel,f{1});
        end
        end
        end
        DATA = LoadCellFile(DATA);
        AddCellMenu(DATA);
        DATA = LoadComments(DATA);
        if DATA.auto.uselastcluster && length(DATA.LastClusters) >= DATA.probe(1)
            DATA.cluster = DATA.LastClusters{DATA.probe(1)};
            forcecluster = 1;
        end
end

DATA.idstr = IDStr(DATA);

cname = ClusterFile(DATA.name,'clusterlog','subdir',DATA.clustersubdir);
cname = strrep(cname,'ClusterLog',sprintf('ClusterLogExpt%d',DATA.exptno));
if isfield(DATA,'logname') && ~strcmp(DATA.logname,cname) && DATA.logfid > 0
  try
      fclose(DATA.logfid);
  catch
      fprintf('ERROR closing logfid\n');
  end
    DATA.logfid = -1;
end
if DATA.logfid < 0 && DATA.nolog == 0
    DATA.logfid = fopen(cname,'a');
    DATA.logname = cname;
    if verbose
    disp(sprintf('Opened Log %s id=%d\n',cname,DATA.logfid));
    end
elseif verbose
    disp(sprintf('No Log set (%s). logfid is %d',IDStr(DATA),DATA.logfid));
end
res.logfid = DATA.logfid;

if DATA.logfid > 0
    try
        fprintf(DATA.logfid,'Ex%d (fid%d) Args %s\r\n',Vall.exptno,DATA.logfid,callstring);
    catch
        fprintf('Invalid handle %d. Opening %s again\n',DATA.logfid,cname);
        DATA.logfid = fopen(cname,'a');
    end
end

if newdata == 1 && isfield(DATA,'name');
    DATA.msg = {};
    DATA.plotspk.subtrigger = 0;
    DATA.plotspk.submean = 0;
    DATA.plotspk.submax = 0;
    DATA.plotspk.submin = 0;
    DATA.plotspk.oneprobe = 0;
    DATA.plotspk.includeprepost = 0;

    DATA.currenttrial = 1;
    DATA.plotspk.muscale = muscale;
    DataClusters = {};
    DATA.clst = [];
    DataClusters = LoadDataClusters(DATA);
elseif newdata == 2 %%Keep DataClusters if newdata == 2 - new probe
else
    DataClusters = {};
end

if newdata == 1 || autocutone;
    if spkrate && setnspk == 0
       setnspk = round(DATA.duration .* spkrate);
    else
        spkrate = setnspk./DATA.duration;
    end
end


if recluster && length(DataClusters) < DATA.probelist(DATA.probe(1)) && forcecluster == 0
    if length(DATA.RefClusters) >= DATA.probe(1)
        userefcluster = 1;
        fprintf('DataClusters is empty using refcluster\n', DATA.probelist(DATA.probe(1)));
    else
        
    DATA = AddErr(DATA,'Cant Apply  Cluster for %d - DataClusters only %d elements\n', DATA.probelist(DATA.probe(1)),length(DataClusters));
    recluster = 0; %just do simple cut if no cluster defined
    savespikes = 0;
    end
end

if forcecluster && isempty(DataClusters)
end
setappdata(DATA.toplevel,'Clusters',DataClusters);
if recluster
%    DATA.currentcluster = 1;
    p = ProbeNumber(DATA);
    DATA.recluster = recluster;
    if isfield(DATA,'cluster') && isfield(DATA.cluster,'auto') && DATA.cluster.auto ==1
        wasauto = 1;
    else
        wasauto = 0;
    end
    if recluster == 4 %load event times, dont calc trigger
        cfile = ClusterFile(DATA.name,DATA.Expt,'details','subdir',DATA.clustersubdir);
    end
    if forcecluster == 0 %%used saved cluster
        
        if (userefcluster || wasauto) && length(DATA.RefClusters) >= DATA.probelist(DATA.probe(1))
            DATA.cluster = DATA.RefClusters{DATA.probelist(DATA.probe(1))};
            DATA.cluster.manual = 4;
            fprintf('Using RefClusters for probe %d\n',p)
        else
            DATA.cluster = DataClusters{DATA.probelist(DATA.probe(1))};
        end
        if ~isfield(DATA.cluster,'minvar')
            DATA.cluster.minvar = 0;
        end
        if (recluster == 2 ||newdata == 2) && isfield(DATA.cluster,'Trigger');
            th = DATA.cluster.Trigger;
        end
    end
    DATA.cluster = FixCluster(DATA.cluster);
    if ~isfield(DATA.cluster,'space')
        fprintf('Cant recluster  - old cluster is empty\n');
        return;
    end

%if applying an old boundary to new data, use the old trigger only if none
%was set
   if recluster ==1 && th == 0
       th = DATA.cluster.Trigger;
   end
   if isnan(th)
       th = 0;
   end
    tryall = 0;
    if isempty(spts)
        spts = DATA.cluster.spts;
    end
%only set CSd/DVDY for PC calculation if the cluster is in that space
%otherwise this just wastes time
%check - could go wrong if cluster 2 is cut in PC space
    if isfield(DATA.cluster,'csd') && ismember(DATA.cluster.space(1), [1 6])
        DATA.csd = DATA.cluster.csd;
        DATA.dvdt = DATA.cluster.dvdt;
    end
%if there is smoothing on the command line, it supercedes the save cluster
    if isfield(DATA.cluster,'tsmooth') && smoothsd == 0
        smoothsd = DATA.cluster.tsmooth;
    end
    if isfield(DATA.cluster,'trigdt') && DATA.trigdt == 0
        DATA.trigdt = DATA.cluster.trigdt;
    end
    if isfield(DATA.cluster,'triggerchan') && length(DATA.cluster.triggerchan) > 1
        ispk = DATA.cluster.triggerchan;
        if isfield(DATA.cluster,'triggertype')
            id = strmatch(DATA.cluster.triggertype,{'sum' 'or' 'sumandreplace'});
            if ismember(id,[1 3])
                addch = 1;
            end
        else
            addch = 1;
        end
    end

    if ismember(recluster, [1 2]) 
        if DATA.cluster.minenergy > 0
        th = DATA.cluster.Trigger;
        minenergy = DATA.cluster.minenergy;
        minvar = DATA.cluster.minvar;
        end
        if isfield(DATA.cluster,'vsmps')
            DATA.vsmps = DATA.cluster.vsmps;
        end
    end
    if recluster == 2
        DATA.clst = [];
    end
    if ~isfield(DATA.cluster,'excludetrialids')
        DATA.cluster.excludetrialids = [];
    end
else
    DATA.recluster = 0;
end

if ~isfield(DATA,'cluster') || ~ClusterIsSet(DATA.cluster, DATA.currentcluster)
    DATA.currentcluster = 1;
end


if ~isfield(DATA.Expt,'Header')
    DATA.cluster.exptreadmethod = -1;
elseif isfield(DATA.Expt.Header,'ReadMethod')
    DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;
else
    DATA.cluster.exptreadmethod = 0;
end



tt(ttn).time = now;
tt(ttn).strs = 'Start Trig';
ttn= ttn+1;
DATA.triggertype = 'one';
p = ProbeNumber(DATA);
if forcetrigger
    ispk = forcetrigger;
end

if ~isfield(Vall,'Spikes')
if verbose >1
    fprintf('Making Trig reference for E%dP%d %s\n',DATA.exptno,ProbeNumber(DATA),datestr(now,'HHMM:ss'));
end
if length(ispk)  > 1 && addch
    rV = mean(double(Vall.V(ispk,:)));
    if isinteger(Vall.V)
        rV = rV .* Vall.intscale(1)./Vall.intscale(2);
    end
    DATA.triggertype = 'sum';
    if addch == 2
        DATA.triggertype = 'sumandreplace';
        Vall.V(ispk(1),:) = rV;
    end
else
    if length(ispk) > 1
        DATA.triggertype = 'or';
    end
    if isinteger(Vall.V)
        rV = double(Vall.V(ispk,:))  .* double(Vall.intscale(1))./Vall.intscale(2);
    else
        rV = Vall.V(ispk,:);
    end
end

if length(smoothsd) > 1 %% explose effect of smoothing on the
if verbose >1
    fprintf('Smoothing  %s\n',datestr(now));
end
    colors = mycolors;
    SetFigure(DATA.tag.vhist, DATA);
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
elseif smoothsd > 0.1
        smv = smooth(rV,smoothsd,'gauss');
else
    smoothsd = 0;
    smv = rV;
end
clear rV; 
if isinteger(smv)
  smv = double(smv) .* Vall.intscale(1)./Vall.intscale(2);
end
if subtractadj && ispk > 1 && ispk <= DATA.nprobes
    cp = [ispk-1 ispk+1];
    subv = mean(Vall.V(cp,:)) .* Vall.intscale(1)./Vall.intscale(2);
    smv = smv - subv;
    Vall.V = double(Vall.V) .* Vall.intscale(1)./Vall.intscale(2);;
    Vall.V(ispk,:) = Vall.V(ispk,:) - subv;
end
if DATA.trigdt
    if DATA.trigdt == 3 %10pt moving average of dvdt = energy
        smv = smooth(diff(smv).^2,10);
    elseif ismember(DATA.trigdt,[4 5]) %convolve with template
        if recluster && isfield(DATA.cluster,'TriggerTemplate')
            DATA.TriggerTemplate = DATA.cluster.TriggerTemplate;
        end
        if size(DATA.TriggerTemplate,1) > length(ispk)
            for j = 1:length(ispk)
            x(j,:) = conv(double(Vall.V(ispk(1),:)),fliplr(DATA.TriggerTemplate(ispk(j),:)));
            end
            x = mean(x,1);
            DATA.TriggerTemplate = DATA.TriggerTemplate(ispk,:);
        elseif length(ispk)  > 1
            for j = 1:length(ispk)
            x(j,:) = conv(Vall.V(ispk(j),:),fliplr(DATA.TriggerTemplate(j,:)));
            end
            x = mean(x,1);
        else
            if DATA.trigdt == 4
                x = conv(smv,fliplr(DATA.TriggerTemplate));
            elseif DATA.trigdt == 5
                x = conv(diff(smv),diff(fliplr(DATA.TriggerTemplate)));
            end
        end
        np = round(length(DATA.TriggerTemplate)/2);
        if mod(length(DATA.TriggerTemplate),2) == 1
            smv = x(np-1:end-np);
        else
            smv = x(np:end-np);
        end
    elseif DATA.trigdt == 2 %% acceleration
        smv = diff(smv,2);
    else
        smv = diff(smv);
    end
end
DATA.triggersmooth = smoothsd;
DATA.triggerchan = ispk;


if smoothv && newdata
    for j = 1:size(Vall.V,1)
        Vall.V(j,:) = smooth(double(Vall.V(j,:)),smoothv,'gauss');
    end
end


tt(ttn).time = now;
tt(ttn).strs = 'Start Trig';
ttn= ttn+1;
if verbose > 1
    fprintf('%s %s %s\n',IDStr(DATA),tt(ttn-1).strs,datestr(now,'HHMM:ss'));
end
DATA.nprobes = size(Vall.V,1);
if ~isfield(DATA,'cluster')
    DATA.cluster = [];
end
F = SetFigure(DATA.tag.top,DATA);
DATA.toplevel = F;
res.toplevel = F;

if DATA.interactive >= 0
    set(F,'UserData',DATA);
    setappdata(F,'Vall',Vall);
    if newdata
        if exist('AutoClusters','var')
            setappdata(DATA.toplevel,'AutoClusters',AutoClusters);
        end
    end
else
    DATA.Vall = Vall;
end
%load trig times before thresholding current channel. Then set trig times
%for current channel to match new values after threshold
if gettrigtimes 
    if gettrigtimes == 2
        [DATA, DataClusters] = LoadTrigTimes(DATA,1:DATA.nprobes,'savexy');
    else
        [DATA, DataClusters] = LoadTrigTimes(DATA,1:DATA.nprobes);
    end
    if recluster == 2
        DATA.lastxy = DATA.xy{1};
        DATA.oldtrigtimes = Vall.t(DATA.trigtimes{DATA.probe(1)});
        DATA.oldclst = DATA.clst;
        PlotXY(DATA.xy{1},DATA.clst);
        th = DATA.cluster.Trigger;
    end
elseif recluster == 2
    if userefcluster == 0
    [DATA, DataClusters, oldxyloaded] = LoadTrigTimes(DATA,[]);
    if forcecluster == 0 && isfield(DataClusters{p},'Evec')
        DATA.cluster.Evec = DataClusters{p}.Evec;
    end
    setappdata(DATA.toplevel,'Custers',DataClusters);
    end
    if oldxyloaded && length(DATA.trigtimes) >= DATA.probe(1)
        SetFigure(DATA.tag.oldxy, DATA);
        DATA.lastxy = DATA.xy{1};
        if DATA.trigtimes{DATA.probe(1)}(end) > Vall.t(end)
            DATA.oldtrigtimes = Vall.t(DATA.trigtimes{DATA.probe(1)});
        else
            DATA.oldtrigtimes = DATA.trigtimes{DATA.probe(1)};
        end
        DATA.oldclst = DATA.clst;
        if DATA.interactive > 0
            PlotXY(DATA.xy{1},DATA.clst);
        end
    elseif DATA.interactive > 0 
%faster to just use recluster 1. Good if the was called from GUI
% But need 2 if doing batch reclassify and want quantification
          recluster = 1;
    end
        th = DATA.cluster.Trigger;
end

if isempty(spts)
    spts = [-12:27];
end

if recluster == 4  %use trig times from saved cluster
    DATA.xy = {[]}; %in case its not set in LoadTrigTimes
    [DATA, DataClusters] = LoadTrigTimes(DATA,DATA.probe(1));
    DATA.trigcheck(DATA.probe(1)) = 1;
    DATA.Trigger = DATA.cluster.Trigger;
    DATA.autoth = 1;
    DATA.spts = DATA.cluster.spts;
    DATA.setnspk = DATA.cluster.nspks;
    DATA.gmcid = [];
    id = DATA.trigtimes{DATA.probe(1)};
else
%in ver 1.1 and before, some triggers (trigdt) were save to clsueter without
%converting to double
    if (th > max(smv) || th < min(smv)) && DATA.cluster.version <= 1.1
        a = th;
        th = th .* Vall.intscale(1)./Vall.intscale(2);
        DATA = AddErr(DATA,'Converting Trigger to double %.0f -> %.3f',a, th);
    end
    DATA.Trigger = th;
    DATA.autoth = 1;
    DATA.spts = spts;
    DATA.thsign = thsign;
    if setnspk > 0
        DATA.setnspk = setnspk;
    end
    if recluster == 1 && matcheventcounts
        DATA.Trigger = 0;
        DATA.thsign = (1+sign(th))/2;
        if isfield(DATA.cluster,'eventrate')
            DATA.setnspk = round(DATA.cluster.eventrate .* length(smv));
        else
            DATA.setnspk = DATA.cluster.nspks;
        end
    end
    [id,  DATA.Trigger, nev] = TriggerV(DATA, smv);
    if length(id) > 1e6
        if DATA.interactive  < 0 
            if ~isfield(DATA,'chspk')
                DATA.chspk = 0;
            end
            fprintf('E%dP%d> 1 million events. (Trigger %f) Exiting AllVPcs\n',DATA.exptno,DATA.chspk(1),DATA.Trigger);
            return;
        end
        yn =  questdlg('>1Million events. Proceed?','trigger','Yes','No','No');
        if strcmpi(yn,'no')
            return;
        end
    end
    DATA.cluster.eventrate = nev./length(smv);
end
if spkrate > 0
    DATA.spkrate = spkrate;
else
    DATA.spkrate = DATA.cluster.eventrate./Vall.samper;
end

clear sgn;
%remove any spikes at very beginning or end where there isn't enough
%data to include the whole spike
if size(id,1) > 1
    id = id';
end

if verbose > 0 fprintf('E%dP%d %dEvents\n',DATA.exptno,ProbeNumber(DATA),length(id)); end
ignoreid = [];
missedtrials = [];
if isfield(DATA, 'Expt') && isfield(DATA.Expt,'Trials') && DATA.usetrials
    if isfield(DATA.Expt.Header,'expname')
        res.expname = DATA.Expt.Header.expname;
    end
    iid = [];
    piid = [];
     if isfield(DATA.cluster,'excludetrialids')
         xct = DATA.cluster.excludetrialids;
     else
         xct = [];
     end
     missedtrials = [];
     blkend = (Vall.blkstart + Vall.blklen.*Vall.samper).*10000;
     blkstart = Vall.blkstart .* 10000;

    for j = 1:length(DATA.Expt.Trials)
        tid = find(blkstart < DATA.Expt.Trials(j).Start(1) & blkend > DATA.Expt.Trials(j).End(end));
        if isempty(tid)
                missedtrials = [missedtrials DATA.Expt.Trials(j).id];
        end
        oid = find(vt(id) > DATA.Expt.Trials(j).Start(1)./10000 - DATA.preperiod & ...
            vt(id) < DATA.Expt.Trials(j).End(end)/10000 + DATA.postperiod);
        pid = find(vt(id(oid)) > DATA.Expt.Trials(j).End(end)/10000); %in postperiod
        iid = [iid oid];
        piid = [piid oid(pid)];
    end
    ntrials = j;
    uid = unique(iid);
    if length(DATA.clst) == length(id)
        DATA.clst = DATA.clst(uid);
        if recluster == 4
            DATA.xy{1} = DATA.xy{1}(uid,:);
        end
    end
    ignoreid = setdiff(id,id(uid));
    id = id(uid);
    [a, pid] = ismember(piid,uid);
DATA.postevents = pid; %index to event index, not FullV times
if isfield(DATA.Expt.Header,'trialdur')
DATA.duration = (DATA.Expt.Header.trialdur+ntrials*(DATA.preperiod+DATA.postperiod))./10000;
else
    disp(sprintf('Header Missing trialdur E%dP%d\n',DATA.exptno,DATA.probe(1)));
    DATA.duration = ntrials * (2+DATA.preperiod+DATA.postperiod);
end
setnspk = DATA.duration * spkrate;
else
DATA.postevents = [];
end
DATA.missedtrials = missedtrials;
if length(missedtrials)
    DATA = AddErr(DATA,'Missing Trials%s\n',sprintf(' %d',missedtrials));
end

if maxrate > 0
    maxspikeset = DATA.duration .* maxrate;
end


if isempty(id)
    DATA = AddErr(DATA,'No Spikes in Expt\n');
    res.t = [];
    res.Trigger = DATA.Trigger;
    return;
end

DATA.rV = smv(:,id);
DATA.trigtimes{DATA.probe(1)} = id;
DATA.trigcheck(DATA.probe(1)) = 1;
res.t = vt(id);
if DATA.plotrv
    allid = repmat(id,length(spts),1) + repmat(spts',1,length(id));
    AllrV = smv(allid);
end


if checklast == 2
    setappdata(DATA.toplevel,'TriggerV',smv);
end

clear smv;
end

DATA.nvpts = length(spts);
res.toplevel = DATA.toplevel;


res.th = th;
allid = [];


res.maxspksused = 0;
if isfield(Vall,'Spikes')
    S = Vall.Spikes;
    AllV(DATA.probe(1),:,:) = double(S.values') .* S.maxv./S.maxint;
    if isfield(S,'xchans')
        AllV(S.xchans,:,:)= double(S.xvalues) .* S.maxv./S.maxint;
        chspk = sort([p S.xchans]);
    else
        chspk = p;
    end
    res.t = S.times./10000;
    spts = DataClusters{p}.spts;
    tpt = find(spts ==0);
    DATA.spts = spts;
    DATA.rV = mean(AllV(DataClusters{p}.triggerchan,tpt,:),1);
    DATA = CopyFields(DATA, DataClusters{p},{'Trigger'});
%    DATA.Trigger = DataClusters{p}.Trigger;
    DATA.triggersmooth = DataClusters{p}.tsmooth;
    DATA.triggerchan = DataClusters{p}.triggerchan;
    DATA.missedtrials = DataClusters{p}.missingtrials;
    DATA.chspk = 1:DATA.nprobes;
    if recluster == 2 %can;'t do reclassify with just spikes
        recluster = 1;
    end
elseif minenergy
    if verbose
        fprintf('Applying min energy at %s\n',datestr(now,'HHMM:ss'));
    end
    allid = repmat(id,length(spts),1) + repmat(spts',1,length(id));
    AllV = reshape(Vall.V(ispk(1),allid),[size(allid)]);
    clear allid;
    energy = sqrt(squeeze(sum(diff(AllV,1,1).^2)));
    spkvar = std(AllV);
    clear AllV;
    ie = [];
    adjustenergy = 0;
   %if the user has set minenergy, don't mess with it. With autocutting seems like sometimes
   %want to adjust minenergy, but not sure when. When it comes up, set
   %adjustenergy
    if length(energy) <= setnspk
        ie = [1:length(energy)];
    elseif adjustenergy
        while length(ie) < setnspk
            ie = find(energy > minenergy & spkvar > minvar);
            minenergy = minenergy * 0.8;
            minvar = minvar * 0.8;
        end
    else
            ie = find(energy > minenergy & spkvar > minvar);
    end
   
    if  length(ie)  > maxspikeset
        if th(1) < 0
            prc = 100 * maxspikeset./length(ie);
        else
            prc = 100-(100 * maxspikeset./length(ie));
        end
        th = prctile(DATA.rV,prc);
    end
    while length(ie) > maxspksallowed 
        minenergy = minenergy * 1.1;
        minvar = minvar * 1.1;
        ie = find(energy > minenergy & spkvar > minvar);
        res.maxspksused = length(ie);
    end
        
    DATA.rV = DATA.rV(ie);
    clear energy;
    clear spkvar;
    clear AllV;
    allid = repmat(id(ie),length(spts),1) + repmat(spts',1,length(ie));
    AllV = reshape(Vall.V(:,allid),[size(Vall.V,1) size(allid)]);
    res.t = vt(id(ie));
    [a,b] = ismember(DATA.postevents,ie);
    DATA.postevents = b;
    id = id(ie);
    DATA.trigtimes{DATA.probe(1)} = id;
else
allid = repmat(id,length(spts),1) + repmat(spts',1,length(id));
AllV = reshape(Vall.V(:,allid),[size(Vall.V,1) size(allid)]);
end
DATA.nevents = size(AllV,3);

if DATA.Trigger(1) <0
    res.xsd = std(DATA.rV(1,:));
else
end

if isinteger(AllV)
    AllV = double(AllV) .* Vall.intscale(1)/Vall.intscale(2);
end
res.spkrate = DATA.nevents./DATA.duration;
if verbose > 0 fprintf('E%dP%d Built AllV %d events\n',DATA.exptno,ProbeNumber(DATA),DATA.nevents); end

if isfield(Vall,'meanV')
    DATA.meanV = Vall.meanV(allid);
    if isfield(Vall,'meanV') && DATA.addmean
%  slightly faster, but might cause memory swapping
%        mV = shiftdim(repmat(Vall.meanV(allid),[1 1 24]),2);
%        AllV = AllV+mV;
        for j = 1:size(AllV,1)
            AllV(j,:,:) = squeeze(AllV(j,:,:)) + Vall.meanV(allid);
        end
    end
else
    DATA.meanV = squeeze(mean(AllV,1));
    DATA.addmean = 1;
end
a = whos('Vall');
DATA.fullvsize = memsize(Vall);
y = max(max(abs(AllV(DATA.probe(1),:,:))));
DATA.voffset = [1:DATA.nprobes].*y;

clear Vall;
if rejectbydiff
    a = max(abs(diff(DATA.meanV)));
    bid = find(a > rejectbydiff);
elseif DATA.addmean && DATA.loadfromspikes == 0
    allbid = [];
    gid = 1:size(AllV,3);
    avar = sum(DATA.meanV(:,gid).^2);
    bid = find(avar > prctile(avar,99) * 2);
    while length(bid)
        allbid = [allbid bid];
        gid = setdiff(1:size(AllV,3),allbid);
        avar = std(DATA.meanV(:,gid));
        avar = sum(DATA.meanV(:,gid).^2);
        bid = find(avar > prctile(avar,99) * 2);
        bid = gid(bid);
    end


    avar = sum(abs(diff(DATA.meanV(:,gid))));
    bvar = smooth(avar,10);
    bid = find(avar > prctile(avar,99) * 1.5);
    while length(bid)
        allbid = [allbid bid];
        gid = setdiff(1:size(AllV,3),allbid);
        avar = sum(abs(diff(DATA.meanV(:,gid))));
        bid = find(avar > prctile(avar,99) * 1.5);
        bid = gid(bid);
    end
    bid = allbid;
elseif size(AllV,1) == 24
avar = squeeze(std(mean(AllV(1:16,:,:))));
bvar = squeeze(std(mean(AllV(17:24,:,:))));
%[da, aid] = sort(avar);
%[db, bid] = sort(avar);
%avar and bvar are really completely correlated because the channels sum to 0
bid = find(avar > prctile(avar,99) * 2 & bvar > prctile(bvar,99) * 2);
else
    bid = [];
end
clear avar;
clear bvar;
if length(bid)
    DATA.artifacttimes = res.t(bid);
    fprintf('%s Made AllV  %d spikes. Removed %d suspicous events at %s,\n',IDStr(DATA),size(AllV,3),length(bid),datestr(now,'HHMM:ss'));
    if DATA.interactive >= 0
    SetFigure(DATA.tag.vare,DATA, DATA.watcharg{:});
    if DATA.nprobes == 24
    hold off;
    plot(squeeze(mean(AllV(17:24,:,bid))));
    hold on;
    plot(squeeze(mean(AllV(1:16,:,bid))));
    else
        plot(DATA.meanV(:,bid));
    end
    drawnow;
    end
    gid = setdiff(1:size(AllV,3),bid);
    [a,b] = ismember(DATA.postevents,gid);
    DATA.postevents = b;
    if length(DATA.clst) == size(AllV,3)
        DATA.clst = DATA.clst(gid);
    end
        
    AllV = AllV(:,:,gid);
    DATA.meanV = DATA.meanV(:,gid);
   
    allid = allid(gid);
    res.t = res.t(gid);
    if DATA.plotrv
    AllrV = AllrV(:,gid);
    end
    DATA.rV = DATA.rV(gid);
    nerr = nerr+1;
    errs{nerr} = sprintf('Removed %d suspicous events, ',length(bid));
    DATA.trigtimes{DATA.probe(1)} = DATA.trigtimes{DATA.probe(1)}(gid);
    ignoreid = [ignoreid id(bid)];
elseif verbose
    fprintf('E%dP%d Made AllV  %d spikes at %s\n',DATA.exptno,ProbeNumber(DATA),size(AllV,3),datestr(now,'HHMM:ss'));
end
res.trigartifacts = length(bid);
if DATA.plotrv && DATA.interactive >= 0
    setappdata(DATA.toplevel,'AllrV',AllrV);
    clear AllrV ;
end
clear allid;
fullcov = 1;
DATA.nevents = size(AllV,3);
res.spkrate = DATA.nevents./DATA.duration;
if DATA.interactive >= 0
    setappdata(DATA.toplevel,'MeanV', DATA.meanV);
    DATA = rmfields(DATA,'meanV');
end

if recluster == 2  % 'reclassify' = compare with last also
    if checklast && isfield(DATA,'oldtrigtimes')
    [t, ida] = setdiff(DATA.oldtrigtimes,res.t);
    if checklast == 2
        spkfile = SpkFilename(DATA);
        try
          a = load(spkfile);
          if isfield(a.Spikes,'maxv')
              a.Spikes.values = double(a.Spikes.values) * a.Spikes.maxv./a.Spikes.maxint;
          end
          DATA.oldSpikes = a.Spikes;
        catch
            DATA = AddErr(DATA,'Error loading %s\n',spkfile);
        end

        goodid = [];
        for j = 1:length(ida)
            d = min(abs(res.t-DATA.oldtrigtimes(ida(j))));
            if d < 0.0005
                goodid = [ida(j) goodid];
                baddiff(j) = 0;
            else
                baddiff(j) = 1;
            end
            tdiff(j) = d;
        end
        ida = setdiff(ida, goodid);
        t = DATA.oldtrigtimes(ida);
    end
            
    n = length(ida); %number of old times not here
    [a, ida] = ismember(t,vt);
    ida = ida(ida > 0);
    bid = find(~ismember(ida,ignoreid)); %missing and not actively removed
    if isempty(ignoreid)
        ignoreid = 0;
    end
    for j = 1:length(bid)
        xid(1,j) = min(abs(id-ida(bid(j))));
        xid(2,j) = min(abs(ignoreid-ida(bid(j))));
    end
    if length(bid)
        nid = find(xid(1,:) >5 & xid(2,:) > 5);
    else
        nid = [];
    end
% reclustern(1) = #of old trig times not in new list
% reclustern(2) = #of those that are not on the ignore list (ignoreid)
% reclustern(4) = # of events in old trig times
    DATA.reclustern = [n length(bid) length(nid) length(DATA.oldtrigtimes)];
    fprintf('Triggers %d/%d old ones gone. (%d not ignored, %d >0.5ms from closest)\n',n,length(DATA.oldtrigtimes),length(bid),length(nid));
    else
        DATA.reclustern = [0 0 0 0];
    end
    try
            res.VCdatediff = DATA.cluster.savetime(end)-FullVData.builddate;
    catch
    end
end

pcplots = [1 2; ...
           1 3; ...
           1 4; ...
           1 5; ...
           2 3; ...
           2 4; ...
           2 5; ...
           3 5];

%need to reset vpts for new probes, because probe # is added   
vsmps = DATA.vsmps;
DATA.dvpts = [0 9 0 13; ...
           0 9 0 16; ...
           0 9 -1 9; ...
           0 13 1 13; ...
           0 20 -1 14; ...
           0 10 -1 5; ...
           0 10 0 30; ...
           0 9 0 12];

 if newdata == 1
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
    DATA.tmplots(14,:) = [2 16];
    DATA.tmplots(15,:) = [2 17];
    DATA.tmplots(16,:) = [2 18];
    DATA.tmplots(17,:) = [1 2];
    DATA.tmplots(18,:) = [3 4];
    DATA.tmplots(19,:) = [1 3];
    DATA.tmplots(20,:) = [2 4];
    DATA.tmplots(21,:) = [1 4];
    DATA.tmplots(22,:) = [2 3];
    DATA.tmplots(23,:) = [1 2];
    DATA.tmplots(24,:) = [3 4];

    if DATA.nprobes == 1
        DATA.tmplots(4,:) = [1 10];
        DATA.tmplots(6,:) = [1 11];
        DATA.tmplots(7,:) = [1 12];
    end
    
    DATA.gmtypes = [1 0 2 3];
    DATA.gmtypelabels = {'PCs', 'Var-E', 'ADC', 'Template' 'Template2' 'Reserved' 'StdTemplate'};
    DATA.pcspace = [1:4];  %N-D fits use first four PCs
    DATA.tmplspace  = [1 2 8 10 5 12];
    DATA.tmplspace(2,:)  = [1 2 3 4 0 0]; %for Template group 2
    DATA.tmplspace(3,:)  = [1 2 3 4 0 0]; %for fixed std templates
    DATA.tmplspace(4,:)  = [1 2 3 4 0 0]; %for fixed std templates
    DATA.vspace = [6 11 15 20];
    DATA.restricttimerange = [];
    DATA.excludetrialids = [];
    DATA.selectprobe = zeros(1,DATA.nprobes);

 end
DATA.minenergy = minenergy;
DATA.minvar = minvar;
DATA.energy = [];
DATA.energy(1,:) = sqrt(squeeze(sum(diff(AllV(ispk(1),:,:),1,2).^2,2)));
DATA.uid = 1:size(AllV,3);


DATA.probe = ispk;
np = size(AllV,1);
nv = size(AllV,2);
DATA.nvpts = nv;
DATA.vpts = SetVsamples(vsmps,ispk(1),np,nv);

DATA.dvpts(:,1) = DATA.dvpts(:,1) + DATA.probe(1);
DATA.dvpts(:,3) = DATA.dvpts(:,3) + DATA.probe(1);
id = find(DATA.dvpts(:,1) < 1);
DATA.dvpts(id,1) = 1;
id = find(DATA.dvpts(:,1) > np);
DATA.dvpts(id,1) = np;
id = find(DATA.dvpts(:,3) < 1);
DATA.dvpts(id,3) = 1;
id = find(DATA.dvpts(:,3) > np);
DATA.dvpts(id,3) = np;
id = find(DATA.dvpts(:,2) > nv);
DATA.dvpts(id,2) = nv;
id = find(DATA.dvpts(:,4) > nv);
DATA.dvpts(id,4) = nv;

 DATA.pcplots = pcplots;
DATA.t = res.t;

DATA.spksperview = 100;
DATA.spklst = 1:100;
DATA.pcprobes = chspk;

id = find(DATA.vspace > length(DATA.spts));
DATA.vspace(id) = length(DATA.spts);
DATA.TemplateLabels = TemplateLabels(DATA,0);
if recluster == 0
    DATA.cluster.spts = DATA.spts;
end

if nprobepc < 0
    DATA.chspk = chspk;
elseif  length(ispk) > 1
    DATA.chspk = min(ispk)-nprobepc:max(ispk)+nprobepc;
elseif nprobepc <= 1
    DATA.chspk = [ispk(1)-1:ispk(1)+1];
    if DATA.chspk(end) > DATA.nprobes
        DATA.chspk(end) = DATA.chspk(1)-1;
    elseif DATA.chspk(1) < 1
        DATA.chspk(1) = DATA.chspk(end)+1;
    end
else
    DATA.chspk = [ispk(1)-nprobepc:ispk(1)+nprobepc];
end
DATA.chspk = DATA.chspk(DATA.chspk > 0 & DATA.chspk <= DATA.nprobes);
if isempty(DATA.plotspk.probes)
DATA.plotspk.probes = DATA.chspk;
end

xspk = DATA.chspk(DATA.chspk ~= ispk(1));
for j = 1:length(xspk)
    DATA.energy(j+1,:) = sqrt(squeeze(sum(diff(AllV(xspk(j),:,:),1,2).^2,2)));
end
F = SetFigure(DATA.tag.top,DATA);
if DATA.interactive >= 0
    setappdata(F,'plotcsd',0);
end

DATA.toplevel = F;
set(F,'UserData',DATA);
%xsd is the SD of all local minima. Closer that 1sd to zero doesnt'
%get explored automatically - roo much risk of running out ouf memory
%clear AllV;

F = SetFigure(DATA.tag.top,DATA);
SetGUI(DATA);
if verbose >1
    fprintf('Cacluating PCs %s %s\n',IDStr(DATA),datestr(now,'HHMM:ss'));
end
DATA.errs = {DATA.errs{:} errs{:}};

if length(tryall) > 1 % when reapplying, just use CSD/Dvdt as before

    np = size(DATA.pcplots,1);
    for j = 1:length(tryall)
        if tryall(j)
        DATA.dvdt = dvdts(j); 
        DATA.csd = csds(j);
        [Cs{j}, Evec{j}, pcs{j}, dips{j}, chsppks{j}, errs, fits{j}] = CalcPCs(DATA,AllV,nprobepc);
        DATA.errs = {DATA.errs{:} errs{:}};
        else
            dips{j} = 0;
        end
    end
    res.alldips = cell2mat(dips');
    if DATA.usebmi
    [a,b] = max(max(res.alldips(:,1:np)'));
    [c,d] = max(res.alldips(:,9)); %mahal distance for first 4 
    else
        [c,d] = max(res.alldips);
        b = d;
    end
    if d ~= b  %different answer
        cb = max(res.alldips(d,1:8));
        am = res.alldips(b,9);
        DATA.msg = {DATA.msg{:} sprintf('Best PC type %d (bm%.2f,mahal%.1f) or %d(%.2f,%.1f)',b,a,am,d,cb,c)};
        if verbose > 1
        fprintf('%s\n',DATA.msg{end});
        end
        if b > 2 && a < 0.25
            b = d;
        end
    end
    res.pcs = pcs{b};
    res.dipvals = dips{b};
    DATA.dvdt = dvdts(b);
    DATA.csd = csds(b);
    DATA.Evec = Evec{b};
    DATA.pcfit = fits{b};
    res.Evec = Evec{b};
%    DATA.chspk = chspk;
    C = Cs{b};
    DATA.alldips = res.alldips;
else
    if recluster %no need to calcuate GMfits etc
        [C, res.Evec, res.pcs, dip, chspk, errs] = CalcPCs(DATA,AllV,nprobepc);
        DATA.errs = {DATA.errs{:} errs{:}};
        if isfield(DATA.cluster,'Evec') && size(res.Evec.Evec,1) == size(DATA.cluster.Evec.Evec,1)
            sgn = sign(sum(res.Evec.Evec .* DATA.cluster.Evec.Evec,1));
            nid = find(sgn < 0);
            res.Evec.Evec(:,nid) = -res.Evec.Evec(:,nid);
            res.pcs(:,nid) = -res.pcs(:,nid);
        end
        DATA.pcfit.took = 0;
        res.dipvals = 0;
    else
        [C, res.Evec, res.pcs, dip, chspk, errs, DATA.pcfit] = CalcPCs(DATA,AllV,nprobepc);
        res.dipvals = dip;
    end
DATA.pcs = res.pcs;
DATA.Evec = res.Evec;
DATA.errs = {DATA.errs{:} errs{:}};
end
if isfield(res,'alldips')
DATA.dipvals = res.alldips;
else
    DATA.dipvals = res.dipvals;
end
DATA.pcs = res.pcs;
DATA.spkvar = [];
for j = DATA.chspk
DATA.spkvar(j,:) = squeeze(std(AllV(j,:,:),[],2));
end
if DATA.profligate
res.edip = HartigansDipTest(sort(DATA.energy(1,:)));
end
DATA.nprobepc = nprobepc;

DATA.cluster.trigparams.preperiod = DATA.preperiod;
DATA.cluster.trigparams.postperiod = DATA.postperiod;
DATA.clid = [];
DATA.nid = [];
y = max(max(abs(AllV(DATA.probe(1),:,:))));
DATA.voffset = [1:size(AllV,1)].*y;
DATA.plotcsd = 0;
xres = rmfield(res,'pcs');
res.chspk = DATA.chspk;
res.memsz = [0 DATA.fullvsize];
str = 'manual';
if saveautocut || autocutone
    str = 'auto';
end
if verbose
    fprintf('Making Cut %s %s (%s)\n',IDStr(DATA),datestr(now,'HHMM:ss'),ClusterFile(DATA.name,DATA.Expt,str,'subdir',DATA.clustersubdir));
end

F = SetFigure(DATA.tag.top,DATA);
DATA.toplevel = F;
if DATA.interactive >= 0
    setappdata(F,'AllV',AllV);
else
    DATA.AllV = AllV;
end

if saveautocut == 1
    [E, res.cluster] = CutAndSave(DATA,'refine');
    DATA = get(DATA.toplevel,'UserData');
    DATA.MeanSpike = res.cluster.MeanSpike;
    res.chspk = DATA.chspk;
    res.cluster.minspke = prctile(DATA.energy(1,DATA.clid),1) .* 0.95;
    res.cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),DATA.clid),1) .* 0.95;
    set(DATA.toplevel,'UserData',DATA);
    if DATA.savespikes 
        SaveSpikes(DATA, DATA.savespkid);
    end
    tt = TimeMark(tt,'Finish');
    res.times = tt;
    res.cluster = SmallCluster(res.cluster);
    res.memsz = [memsize(DATA) DATA.fullvsize];
    return;
end



if autocutone && recluster == 2
    autocutone = 2;  %call autocut after reclassify
end
if autocutone ==1
    DATA.cluster = FixCluster(DATA.cluster);
    if strcmp(DATA.autocutmode,'quick')
        DATA.cluster.probe = ProbeNumber(DATA);
        DATA = QuickAutoCut(DATA);
        res.cluster = DATA.cluster;
        E = BoundaryFromCluster([],DATA.cluster,1);
    else
        [E, res.cluster] = CutAndSave(DATA,'nosave','refine');
    end
    DATA = get(DATA.toplevel,'UserData');
elseif recluster
    plottype = 0;  %don't do other things after classify
    DATA.cluster = FixCluster(DATA.cluster);
    if isfield(DATA.cluster,'excludetrialids') && ~isempty(DATA.cluster.excludetrialids)
        DATA.excludetrialids = DATA.cluster.excludetrialids;
        SetTrialList(DATA);
        if length(DATA.clst) < size(DATA.pcs,1);
            DATA.clst = ones(size(DATA.pcs,1),1);
        end
        DATA = ExcludeTrials(DATA,DATA.cluster.excludetrialids,0);
        DATA = SetPCs(DATA, 0);
    elseif isfield(DATA.cluster,'restricttimerange') && ~isempty(DATA.cluster.restricttimerange)
        DATA = RestrictTimeRange(DATA,DATA.cluster.restricttimerange);
        DATA = SetPCs(DATA,0);
    end
    if recluster == 3  %just use template
        if length(DATA.forceclusters) > 1
            for j = 1:length(DATA.forceclusters)
                DATA.clid = [];
                DATA.cluster = DATA.forceclusters{j};
                DATA = ReClassify(DATA,'template');
                res.cluster{j} = rmfield(DATA.cluster,{'r' 'bestcl'});
            end
        else
            DATA = ReClassify(DATA,'template');
            res.cluster = DATA.cluster;
        end
    elseif recluster == 4 %use previous ids
        DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(BoundaryFromCluster([],DATA.cluster, DATA.currentcluster));
        needtemplate = NeedTemplateForCluster(DATA.cluster,1);
        DATA.plottype = WhichPlotType(DATA.cluster, DATA.currentcluster);
        if needtemplate & (~isfield(DATA,'TemplateScores') ||newdata)
            DATA = CalcTemplateScores(DATA);
        elseif isfield(DATA,'TemplateScores')
            DATA = rmfields(DATA, {'TemplateScores' 'TemplateUsed'});
        end
        if DATA.cluster.space(1) == 6
            DATA.cluster.shape = 2;
            [DATA.ndxy, a, b] = ProjectND(DATA,DATA.cluster.space(2),DATA.cluster.gmfit);
            if b.err & ~isfield(DATA,'TemplateScores')
                DATA = CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
                DATA.plottype = 3;
                res.err = 1;
            end
        end
        res.cluster = SmallCluster(DATA.cluster);
        if 0 %% old way. Why need to recalc Mean?? Quicker switching without
        DATA.cluster.MeanSpike = PlotMeanSpike(DATA,'recalc','cluster',1);
        id = unique(DATA.clst);
        for j = 3:length(id)
            DATA.cluster.next{j-2}.MeanSpike = PlotMeanSpike(DATA,'recalc','cluster',j-1);
        end
        else
        id = unique(DATA.clst);
        for j = 2:length(id)
            PlotMeanSpike(DATA,'cluster',j-1);
        end
        end


        DATA.clid = find(DATA.clst == 2);
        DATA.nid = find(DATA.clst == 1);
        ReplotPCs(DATA,[]);
        PlotTriggerHist(DATA,DATA.cluster,'quick');
        if DATA.plotspk.showmean
            PlotSpikes(DATA,1);
        end
    elseif recluster == 2  || recluster == 1 %space, but new data
        res.err = 0;
        DATA.cluster.exptno = DATA.exptno;
        DATA.cluster.clst = ones(DATA.nevents,1);
        [iscll, cellid] =  isacell(DATA, DATA.exptno, ProbeNumber(DATA));
        if iscll
            a = find(cellid > 0);
            res.cell = cellid(a(1));
            res.cellid = cellid;
        else
            res.cell = 0;
        end
        if ~isfield(DATA.cluster,'clusterprog')
            DATA.cluster.progversion = 0;
            DATA.cluster.clusterprog = 'AllVPvs';
        elseif ~isfield(DATA.cluster,'progversion')
            id = strfind(DATA.cluster.clusterprog,' ');
            if isempty(id)
                DATA.cluster.progversion = 0;
            else
                DATA.cluster.progversion = sscanf(DATA.cluster.clusterprog(id(end):end),'%f');
            end
        end
        if readclusterfromlog
            DATA = ReadFromLog(DATA);
        end
        if recluster == 2 && oldxyloaded
            oldid = find(DATA.oldclst == 2);
            xcl = [];
            goodxcl = ones(size(DATA.excludetrialids));
            for j = 1:length(DATA.cluster.excludetrialids)
                t = find([DATA.Expt.Trials.id] == DATA.cluster.excludetrialids(j));
                if length(t) == 1
                    id = find(DATA.oldtrigtimes > DATA.Expt.Trials(t).Start(1)./10000 - DATA.preperiod & ...
                        DATA.oldtrigtimes < DATA.postperiod+DATA.Expt.Trials(t).End(end)./10000);
                    xcl = [xcl id];
                else
                    fprintf('Id %d no longer in Expt!\n',DATA.cluster.excludetrialids(j));
                    goodxcl(j) = 0;
                end
            end
            DATA.cluster.excludetrialids = DATA.cluster.excludetrialids(find(goodxcl));
            DATA.excludetrialids = DATA.cluster.excludetrialids;
            oldid = setdiff(oldid,xcl);
        end
        if DATA.cluster.space(1) ~= 6 && DATA.cluster.shape == 2
            DATA = AddErr(DATA,'Error Shape 2 but space %s\n',sprintf('%d ',DATA.cluster.space));
            DATA.cluster.shape = 1;
        end
        
        if DATA.cluster.auto == 2 %was set in plotcluster
            DATA.cluster.auto = 0;
        end
        E = BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);
        if DATA.profiling
            fprintf('Return from BoundaryCluster %.4f at %.4f\n',mytoc(E.endtime),mytoc(tt(1).time));
        end
        DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(E);
        [needtemplate, plottypes] = NeedTemplateForCluster(DATA.cluster,1);
        DATA.plottype = WhichPlotType(DATA.cluster, DATA.currentcluster);
        if needtemplate & (~isfield(DATA,'TemplateScores') ||newdata)
            if isfield(DATA.cluster,'TemplateUsed')
                if size(DATA.cluster.TemplateUsed,1) == 2 && DATA.plottype ~= 7 && 0
                    DATA.plottype = 7;
                    DATA.cluster.space = [6 7];
                    DATA.cluster.shape = 2;
                end
                DATA = SetTemplateData(DATA,DATA.currentcluster);
            elseif sum(isnan(DATA.cluster.MeanSpike.ms(DATA.probe(1),:)))
                DATA.clst = ones(length(DATA.t(DATA.uid)),1);
                [a, id] = find(ismember(DATA.t,DATA.oldtrigtimes(oldid)));
                if length(id)
                    DATA.clst(id) = 2;
                    TemplatePlot(DATA,'nodip');
                    DATA = get(DATA.toplevel,'UserData');
                else
                    DATA = AddErr(DATA,'NAN Meanspike\n');
                    DATA.cluster.space = [1 1 2];
                end
            else
                DATA = CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
%used to call this lager and remove filesd form dta here.  WHy?   can cause
%ProjectND to fail
%            elseif isfield(DATA,'TemplateScores')
%                DATA = rmfields(DATA, {'TemplateScores' 'TemplateUsed'});
            end
            if ~isfield(DATA.cluster,'TemplateUsed')
                DATA.cluster.TemplateUsed = DATA.TemplateUsed;
            end
        end
        
%if this clsuter previously used predefined eigenvectors, use those again        
        if isfield(DATA.cluster,'forceevec')
            forceevec = DATA.cluster.forceevec;
        else
            forceevec = 0;
        end
        if DATA.cluster.space(1) == 6
            clusterspace = DATA.cluster.space(2);
        else
            clusterspace = DATA.cluster.space(1);
        end
        
        if (ismember(recluster,[1 2]) || forceevec) && clusterspace == 1
            
            if isfield(DATA.cluster,'Evec') && recluster == 2 %recalssify used old PCs. 
            [a,DATA.Evec, DATA.pcs] = CalcPCs(DATA, AllV, DATA.nprobepc, DATA.cluster.Evec);
            DATA.cluster.forceevec = 1;
            else
            [a,DATA.Evec, DATA.pcs] = CalcPCs(DATA, AllV, DATA.nprobepc);
            DATA.cluster.forceevec = 0;
            end
        end
        
        if DATA.cluster.space(1) == 6
            if strncmp(DATA.cluster.clusterprog,'AllVpcs',7) && DATA.cluster.progversion <= 1.2
                DATA.cluster.shape = 2;
            end
            [DATA.ndxy, a, b] = ProjectND(DATA,DATA.cluster.space(2),DATA.cluster.gmfit);
            DATA.xy{1} = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),DATA.cluster.angle);
            res.err = b.err;
            if b.err & ~isfield(DATA,'TemplateScores')
                DATA = CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
                DATA.plottype = 3;
            end
            if b.err ==1 && DATA.cluster.auto == 1 %redo autocut
                PrintMsg(DATA.logfid,'Redoing AutoCut to fix dimension error');
                autocutone = 1;
            end
        elseif DATA.cluster.space(1) == 2 && DATA.cluster.shape ==2 && DATA.cluster.space(2) == 4 %earlier bug
            DATA.ndxy = ProjectND(DATA,DATA.cluster.space(2),DATA.cluster.gmfit);
            DATA.xy{1} = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),DATA.cluster.angle);
        elseif DATA.cluster.space(1) == 1
            DATA.xy{1} = xyrotate(DATA.pcs(:,DATA.cluster.space(2)),DATA.pcs(:,DATA.cluster.space(3)),DATA.cluster.angle);
        elseif DATA.cluster.space(1) == 3
            DATA.xy{1} = xyrotate(DATA.TemplateScores(:,DATA.cluster.space(2)),DATA.TemplateScores(:,DATA.cluster.space(3)),DATA.cluster.angle);
        elseif DATA.cluster.space(1) == 2 && length(DATA.cluster.space) > 4 %adc values
            p = DATA.cluster.space(2:end);
            xy(DATA.uid,1) = AllV(p(1),p(2),DATA.uid);
            xy(DATA.uid,2) = AllV(p(3),p(4),DATA.uid);
            %            x =  squeeze(AllV(a(1),a(2),:));
            DATA.xy{1} = xyrotate(squeeze(AllV(p(1),p(2),:)),squeeze(AllV(p(3),p(4),:)),DATA.cluster.angle);
        end
    flipstr = [];
    flip = 1;
    DATA.cluster = CheckScoreScaling(DATA, DATA.cluster);
    if DATA.cluster.shape > 0 %not an ellipse - check sign
        if DATA.cluster.crit < 0
            ncut = sum(DATA.xy{1}(DATA.uid,1) < DATA.cluster.crit);
            nalt =  sum(DATA.xy{1}(DATA.uid,1) > -DATA.cluster.crit);
        else
            ncut = sum(DATA.xy{1}(DATA.uid,1) > DATA.cluster.crit);
            nalt =  sum(DATA.xy{1}(DATA.uid,1) < -DATA.cluster.crit);
        end
        if recluster == 2 && exist('oldid','var') 
            nspk = length(oldid);
            if ncut < nspk/10 && nalt > nspk/5
                DATA.cluster.crit = DATA.cluster.crit .* -1;
                DATA.cluster.sign = DATA.cluster.sign .* -1;
                flipstr = 'flip';
                flip = -1;
            end
        end
    end

    if autocutone
        [E, scores, dips, xy, details] = AutoCut(DATA);
        DATA = get(DATA.toplevel,'UserData');
        DATA.ndxy = xy;
        DATA.cluster = ClusterFromBoundary(E, DATA.cluster);
        CheckClusters(DATA.cluster,'CheckFitSpace')
        %            DATA.ndxy = ProjectND(DATA, E.space(2), E.gmfit);
    end
    if length(DATA.uid) == 0 && isfield(DATA.cluster,'restricttimerange')
        DATA = UseAllEvents(DATA);
    end


    ctimes(1) = DATA.cluster.ctime;
    for j = 1:length(DATA.cluster.next)
        if isfield(DATA.cluster.next{j},'ctime')
            ctimes(j+1) = DATA.cluster.next{j}.ctime;
        else
            if isempty(DATA.cluster.next{j})
                ctimes(j+1) = 0;
            else
                errordlg(sprintf('Missing ctime Cluster %d',j+1),'Cluster Cut Time','modal');
                ctimes(j+1) = j+1;
            end
        end
    end
    ctimes = ctimes(ctimes > 0);
    [a, tlist] = sort(ctimes);


    %check to see if there are missing "next" clusters
    wasempty = 0;
    if recluster == 2 && DATA.reclustern(1) == 0 && oldxyloaded %trigger times match - can use old clst
        DATA.clst = DATA.oldclst;
        cc = DATA.currentcluster;
        for j = 1:length(DATA.cluster.next)
            if isempty(DATA.cluster.next{j})
                wasempty = 1;
                if sum(DATA.clst ==  j+2) > 1
                    DATA.currentcluster = j+1;
                    DATA.cluster.next{j}.MeanSpike = PlotMeanSpike(DATA,'recalc');
                end
            elseif wasempty
                wasempty = 2;
            end
        end
        DATA.currentcluster = cc;
    end
    if useoldlst == 0
        DATA.clst = [];
    end
    if recluster == 1 && DATA.watchplots < 2
        DATA = ClassifyAll(DATA,1,'quick','recluster');
    else
        DATA = ClassifyAll(DATA,1,'recluster');
    end
    if DATA.profiling
        fprintf('Return from ClassifyAll %.4f at %.4f\n',mytoc(DATA.endtime),mytoc(tt(1).time));
    end

    %if an error messed up the classification, don't copy meanspike
    res.cluster = SmallCluster(DATA.cluster);
    ReplotPCs(DATA,[]);
    if recluster ~= 1 %will have been called in ClassifySpikes
        PlotTriggerHist(DATA,DATA.cluster);
    end
    if DATA.profiling
        fprintf('Replotted at %.4f\n',mytoc(tt(1).time));
    end
    if recluster == 2 && isfield(DataClusters{DATA.probe(1)},'space') && exist('oldid','var')
        oldC = DataClusters{DATA.probe(1)};
        if length(oldC.space) == 1
            oldC.space(2) = 0;
        end
        if max(oldid) > length(DATA.oldtrigtimes)
            oldid = oldid(oldid <= length(DATA.oldtrigtimes));
        end
        %idxa is matching times for classified spikes
        [ts, idxa, idxb] = intersect(round(DATA.t(DATA.clid).*1000), round(DATA.oldtrigtimes(oldid).*1000));
        unsafetosave = 0;
        [ida, iida] = setdiff(round(DATA.t(DATA.clid).*1000), round(DATA.oldtrigtimes(oldid).*1000));
        [idb, iidb] = setdiff( round(DATA.oldtrigtimes(oldid).*1000),round(DATA.t(DATA.clid).*1000));
        if checklast > 1
            tb = DATA.oldtrigtimes(oldid(iidb));
            ta = DATA.t(DATA.clid(iida));
            for j = 1:length(ida)
                [c,d] = min(abs(ta(j)-tb));
                if c > 0.0005
                    goodid(j) = 0;
                else
                    goodid(j) = 1;
                end
            end
            ida = ida(goodid == 0);
            iida = iida(goodid == 0);
            goodid = [];
            for j = 1:length(idb)
                [c,d] = min(abs(tb(j)-ta));
                if c > 0.0005
                    goodid(j) = 0;
                else
                    goodid(j) = 1;
                end
            end
            idb = idb(goodid == 0);
        end
                
        if DATA.interactive >= 0 && ~isempty(DATA.clid) && ~isempty(oldid)
            SetFigure(DATA.tag.oldxy, DATA);
            subplot(2,1,2);
            PlotXY(flip .*DATA.lastxy,DATA.oldclst);
            title(sprintf('Mahal 1D %.2f, 2D %.2f Space%s',oldC.mahal(4),oldC.mahal(1),...
                sprintf(' %d',oldC.space)));
            subplot(2,1,1);
            hold off;
            h(1) = plot(DATA.oldtrigtimes,flip.*DATA.lastxy(:,1),'k.','buttondownfcn',{@ShowFullV, 1}); %old events
            labels{1} = 'old triggers';
            hold on;
            h(2) = plot(DATA.t(DATA.uid),DATA.xy{1}(:,1),'c.','buttondownfcn',{@ShowFullV, 2});
            labels{2} = 'New Triggers';
            h(3) = plot(DATA.oldtrigtimes(oldid),flip.*DATA.lastxy(oldid,1),'.','buttondownfcn',{@ShowFullV, 1}); %old cell
            labels{3} = 'old cell';
            h(4) = plot(DATA.t(DATA.clid),DATA.xy{1}(DATA.clid,1),'r.','buttondownfcn',{@ShowFullV, 2});
            labels{4} = 'New Cell';
            yl = get(gca,'ylim');
            yl(2) = yl(1) + diff(yl)/8;
            if isfield(DATA.Expt,'Trials')
                for j = 1:length(DATA.Expt.Trials)
                    t = DATA.Expt.Trials(j).Start(1)./10000;
                    plot([t t],yl,'k-');
                end
            end
            if ~isempty(idxa)
            h(5) = plot(DATA.t(DATA.clid(idxa)),DATA.xy{1}(DATA.clid(idxa),1),'g.');
            end
            labels{5} = '?';
            if ishandle(h)
                legend(h,labels);
            end
            
            
            SetFigure(DATA.tag.hist, DATA);
            subplot(2,1,2);
            h = ishold;
            hold on;
            if ~isempty(iida)
                plot(DATA.xy{1}(DATA.clid(iida),1),DATA.xy{1}(DATA.clid(iida),2),'x');
            end
            if h == 0
                hold off;
            end
        end

        str = sprintf('P%d last cut at %s: %d/%d or %d/%d spikes are different (%s)',...
            DATA.probe(1),datestr(oldC.ctime),length(ida),...
            length(DATA.clid),length(idb),length(oldid),flipstr);
        res.matchcounts = [length(ida) length(DATA.clid) length(idb) length(oldid) NaN];
        if DATA.reclustern(1) == 0 %trigger times match - can compare lsts
            if length(DATA.uid) == length(DATA.oldclst)
                res.matchcounts(5) = sum(DATA.clst(DATA.uid) ~= DATA.oldclst);
            elseif length(DATA.clst) == length(DATA.oldclst)
                res.matchcounts(5) = sum(DATA.clst ~= DATA.oldclst);
            end
        end

        if abs(length(DATA.clid)-length(oldid))/length(oldid)  > 0.2
            unsafetosave = 1;
            unsafelabels{1} = 'classified id length mismatch';
        end
        PrintMsg(DATA.logfid,str);
        ida = setdiff( DATA.oldtrigtimes,DATA.t);
        fprintf('Triggers %d/%d old ones gone\n',length(ida),length(DATA.oldtrigtimes))
        
        id = find(DATA.cluster.clst ==2);
        missid = setdiff(oldC.times,DATA.t(id));
        res.matchcounts(6) = length(missid);
        nx = min([length(DATA.cluster.next) length(oldC.next)]);
        for k = 1:nx
            if isfield(DATA.cluster.next{k},'times') && isfield(oldC.next{k},'times')
                id = find(DATA.cluster.clst ==2+k);
                missid = setdiff(oldC.next{k}.times,DATA.t(id));
                res.matchcounts(6+k) = length(missid);
            end
        end
        if DATA.interactive >= 0
        SetFigure(DATA.tag.oldxy, DATA);
        title(str);
        SetFigure(DATA.tag.comparexy, DATA);
        hold off;
        plot(DATA.xy{1}(DATA.clid(idxa),1),DATA.lastxy(oldid(idxb),1),'.');
        hold on;
        plot(DATA.xy{1}(DATA.clid(idxa),2),DATA.lastxy(oldid(idxb),2),'r.');
        refline(1);
        end
        if length(idxa) > 1
            xc = corrcoef(DATA.xy{1}(DATA.clid(idxa),1),DATA.lastxy(oldid(idxb),1));
            yc = corrcoef(DATA.xy{1}(DATA.clid(idxa),2),DATA.lastxy(oldid(idxb),2));
            if abs(xc(1,2)) < 0.5
                unsafetosave = unsafetosave+2;
            end
        else
            xc = [];
            if length(oldid) > 5 % if there are no spikes in the older file, there is no problem
                unsafetosave = unsafetosave+2;
            end
        end
        
        if autocutone %'autocut' on command line
            unsafetosave = 0;
        end
        PrintMsg(DATA.logfid,'%s Trigger match %d %d %d %d',IDStr(DATA),DATA.reclustern(1),DATA.reclustern(2),DATA.reclustern(3),DATA.reclustern(4));
 %reclustern(3) is missing triggers that dont have a match at +-0.5ms
        if DATA.reclustern(3)./DATA.reclustern(4) > 0.1
            unsafetosave = unsafetosave+4;
        end
        %manual == 4 means using reflcuster, don't worry about matching old
        %space
        if DATA.cluster.manual ~= 4
            DATA.cluster.auto = DataClusters{DATA.probe(1)}.auto;
        end
        if DATA.cluster.manual ==2 %done in plotcluster but not ye quantified. = Can't check
            unsafetosave = -1000;
        elseif unsafetosave & DATA.cluster.auto == 1  && ~isacell(DATA, DATA.exptno, DATA.probe(1));
            if DATA.savespikes
                fprintf('Mismatched, so redoing autocut and saving\n');
                [E, res.cluster] = CutAndSave(DATA);
            else
                fprintf('Mismatched, so redoing autocut\n');
                [E, res.cluster] = CutAndSave(DATA,'nosave');
            end
            DATA = get(DATA.toplevel,'UserData');
            CheckClusters(DATA.cluster,'CheckFitSpace')
            unsafetosave  = -unsafetosave;
        end
        if wasempty == 2  %empty clusters were saved
            unsafetosave = unsafetosave + 128;
        end
 % sometimes triggerchan has 2 probes but trigger clearly was really just
 % one.  Mark this so can check if necessary
        if unsafetosave > 0 && length(DATA.cluster.triggerchan) > 1
            unsafetosave = unsafetosave + 64;
        end

%if saved cluster is quick, then can't use saved xy,clst
%so just look at times of classified events
        if oldC.quick 
            if sum(res.matchcounts(6:end)) == 0
                fprintf('List of clustered spike times matches\n');
                unsafetosave = 0;
            end
            res.wasquick = 1;
        else
            res.wasquick = 0;
        end
        res.overlapn = length(ts);
        res.cutfraction(1) = length(oldid)./length(DATA.oldclst);
        res.cutfraction(2) = length(DATA.clid)./length(DATA.clst);
        res.trigmatch = DATA.reclustern;
        res.stds  = cat(1,std(DATA.xy{1}),std(DATA.lastxy));

        if length(xc) > 1
            res.xcorr = [xc(1,2) yc(1,2)];
        else
            res.xcorr = [0 0 ];
        end
    elseif recluster == 1 && refinecluster
% Would like to fit old ellipse, then use this starting point to fit GM
% model in this space, and adjust ellipse to catpure those points.  But GM
% model may go off and fit something else entirelys, esp whern there are
% two clusters.
%lemM211 Expt 35 P6 is  a good example
          if length(DataClusters)< DATA.probe(1) || ~isfield(DataClusters{DATA.probe(1)},'mahal') 
                 unsafetosave = 0;
          elseif DataClusters{DATA.probe(1)}.mahal(1) > DATA.cluster.mahal(1)
             if forcecluster == 0
                 unsafetosave = 1;
             end
         elseif DataClusters{DATA.probe(1)}.auto == 0 && DATA.cluster.exptno ~= forceclusterexpt
             %for now only replace auto cuts with refined manual cuts
             %if auto == 0 but manual ==3, means this is a previous
             %"refinement", so overwrite. 
             if isfield(DataClusters{DATA.probe(1)},'manual') && ...
                     DataClusters{DATA.probe(1)}.manual == 1
             unsafetosave = 2;
             end
         end
    end
        res.auto = DATA.cluster.auto;
        res.unsafetosave = unsafetosave;
        if DATA.interactive > 0
            drawnow;
        end
    elseif recluster == 10  %old 2, diabled for now
        DATA = ReClassify(DATA);
        res.cluster = SmallCluster(DATA.cluster);
    else
        DATA = ReClassify(DATA);
        res.cluster = SmallCluster(DATA.cluster);
    end

    E = DATA.clusterboundary{DATA.currentcluster};
    if DATA.profiling
        fprintf('Return from BoundaryCluster %.4f at %.4f\n',mytoc(E.endtime),mytoc(tt(1).time));
    end
    if isfield(E,'space')
    DATA.plottype = WhichPlotType(E,1);
    end
    PlotHistogram(DATA,E, DATA.quickcutmode);
    if DATA.profiling
        fprintf('Return from PlotHistoGram at %.4f\n',mytoc(tt(1).time));
    end
    if recluster == 2 && unsafetosave >= 0 && DATA.interactive > 0 %<0 = redone autocut, so may not match
        subplot(2,1,2);
        hold on;
        plot(DATA.xy{1}(DATA.clid(iida),1),DATA.xy{1}(DATA.clid(iida),2),'gx');
        hold off;
    end
    if DATA.autorefine
        res.refinemode = DATA.refinemode;
    end
    if calcdistancematrix
        res.DistanceMatrix = CalcDistanceMatrices(DATA, calcdistancematrix);
        DATA.errs = {DATA.errs{:} res.DistanceMatrix.errs{:}};
        DATA.errstates = {DATA.errstates{:} res.DistanceMatrix.errstates{:}};
        if DATA.interactive >= 0
            b = res.DistanceMatrix;
            DATA.gmcid = res.DistanceMatrix.cid;
            DATA.usegmcid = 1;
            GetFigure('DistanceMatrix');
            [E,V] = eig(squeeze(b.D(:,:,1)));
            [c,d] = sort(E(:,1));
            D = MatrixPermute(b.D(:,:,1),d);
            imagesc(D);
            caxis([0 5]);
            ReplotPCs(DATA,[]);
%            set(DATA.toplevel,'UserData',DATA);
        end
    end
elseif DATA.interactive >= 0  %this is inteactive for sure
    DATA.watchplots = 1;
    DATA.watcharg = {'front'};
    DATA.interactive = 1;

    plotdprimemax = 1;
%    [d, details] = FindDip(DATA.pcs(:,1),DATA.energy(1,:));
    [d, details] = GMDip(DATA.pcs(:,1:2),DATA.energy(1,:),'label',DATA.idstr);
    dcrit = d(1);
    E = BoundaryFromCluster([],DATA.cluster,DATA.currentcluster);
    E.pos(1) = dcrit;
    E.pos(3) = dcrit;
    E.pos(2)= min(DATA.pcs(:,2));
    E.pos(4)= max(DATA.pcs(:,2));
    p = E.pos;
    E.xyr = [mean(p([1 3])) mean(p([2 4])) abs(diff(p([1 3]))/2) abs(diff(p([2 4]))/2)];
    E.shape = 1;
    E.space = [1 1 2];
    E.h = [];
    E.pcplot = [1 2];
    if dcrit < 0
        DATA.clid = find(DATA.pcs(:,1) < dcrit);
        DATA.nid = find(DATA.pcs(:,1) >= dcrit);
        cid = find(DATA.pcs(:,1) < d(2));
        nid = find(DATA.pcs(:,1) >= dcrit);
    else
        DATA.clid = find(DATA.pcs(:,1) > dcrit);
        DATA.nid = find(DATA.pcs(:,1) <= dcrit);
    end
    res.id = DATA.clid;
%    res.dprime = abs(details.dprime);



    res.dvdt = DATA.dvdt;
    res.csd = DATA.csd;

end


res.errs = DATA.errs;

if DATA.profiling
    fprintf('recluster done at  %.4f\n',mytoc(tt(1).time));
end
if calcclscores
    k = 0;
    for j = 1:length(Clusters)
        if ~isempty(Clusters{j});
            k = k+1;
          res.Clusterscores(k,:) = Clusters{j}.MeanSpike.ms * squeeze(AllV(ispk(1),:,:)); 
        end
    end
end

    if plotsummary
        SetFigure(DATA.tag.allspikes, DATA);
        PlotQuickSpikes(DATA,500);
        PlotAllProbes(DATA,'xy')
    end
if saveclusters == 1 || (saveclusters == 2 && unsafetosave <= 0)
%Calls from teh command line mean that this is automatic clustering. 
%Unless its a recluster
    if saveautocut == 2
        DATA.cluster.auto=1;
        DATA.savespikes = 2;
    end
    if (recluster == 0 && DATA.cluster.auto ~= 1) || recluster >  0 || saveautocut == 2
    outname = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    else
    outname = ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);
    end
    DATA =  SaveClusters(DATA, outname);
    set(DATA.toplevel,'UserData',DATA);
    Cluster = DATA.cluster;
% Can't do this any more - DATA.MeanSpike might be for cl 2
%    cluster.MeanSpike = DATA.MeanSpike;
    if DATA.savespikes
        SaveSpikes(DATA, DATA.savespkid);
    end
end
res.memsz = [memsize(DATA) DATA.fullvsize];


if DATA.csd && DATA.interactive >= 0
    AllCSD = diff(AllV,DATA.csd,1);
    setappdata(DATA.toplevel,'AllCSD',AllCSD);
end


if DATA.interactive >= 0
SetFigure(DATA.tag.covar, DATA);
imagesc(C);
end
DATA.nid = 1:size(AllV,3);

if DATA.readlayout == 2
    ApplyLayout(DATA);
end



if recluster
    if DATA.profiling
        fprintf('QuickSpks at %.4f\n',mytoc(tt(1).time));
    end
    if DATA.quickcutmode.plotspikes && DATA.interactive >=0
        QuickSpks(DATA, 1000);
    end
    fprintf('recluster %s took %.2f (at %s)\n',IDStr(DATA),mytoc(tt(1).time),datestr(now));
    if DATA.profiling
        for j = 1:length(tt)
            fprintf('  %s took %.2f \n',tt(j).strs,mytoc(tt(j).time));
        end
    end
end

if DATA.profiling
    profile viewer
end

if DATA.plot.exptfit > 0 || DATA.plot.expt
    DATA.Expt = PlotExptCounts(DATA);
end
if DATA.elmousept.shape < 0 && DATA.interactive > 0
    DATA = SetEllipseDrawing(DATA, 1, 'cluster', DATA.currentcluster);
end
if plottype == 0
%may not get logfid from GUI if not interactive 
%which means each call will open a new handle. So close here
    DATA = CloseLog(DATA);
    res.errs = {res.errs{:} DATA.errs{:}};
    res.logfid = DATA.logfid;
    res.errstates = DATA.errstates;
    if DATA.interactive >= 0
        drawnow; % force draing of any other figurs to be sure this sets active figure
        set(0,'currentfigure',DATA.toplevel);
        SetGUI(DATA);
        set(DATA.toplevel,'UserData',DATA);
    end
    return;
end


DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
SetFigure(DATA.tag.vare, DATA);
hold off;
c = DATA.probe(1);
plot(DATA.energy(1,:),DATA.spkvar(c,:)./DATA.energy(1,:),'.');


           
set(DATA.toplevel,'UserData',DATA);
ShowTaggedProbes(DATA);



SetFigure(DATA.tag.spikes, DATA);
DATA.plotdvdt = 0;
PlotSpikes(DATA,1:100);

if oldscores & ~isempty(DataClusters{DATA.probe(1)})
    DATA = CalcTemplatesFromMean(DATA, DataClusters{DATA.probe(1)}.MeanSpike);
    DATA.plottype = WhichPlotType(DataClusters{DATA.probe(1)},1);
    E = BoundaryFromCluster(E,DataClusters{DATA.probe(1)}, DATA.currentcluster);
    E.pcplot = DataClusters{DATA.probe(1)}.space(2:end);
    DATA.cluster = DataClusters{DATA.probe(1)};
end

if ~isfield(DATA,'TemplateUsed')
    DATA.TemplateUsed = [];
end

if recluster == 0 && DATA.interactive >= 0
    %if recluster i s0, templates may not have been calcualted.  Don't use
    %existing "TemplateScores" field - might belong to earlier probe
    if ismember(DATA.plottype,[3 4])
        DATA.plottype = 1;
    end
    DATA.cluster.space = [1 1 2];
PlotHistogram(DATA,E,'quick');
if autocutone == 0
    [cl, Cluster] = ClassifySpikes(DATA,E,'quick');
    cl.MeanSpike = PlotMeanSpike(DATA);
    if oldscores == 0
        DATA.cluster = rmfield(Cluster,'r');
    else
    end
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    Cluster.MeanSpike = cl.MeanSpike;
    Cluster.chspk = DATA.chspk;
    Cluster.minspke = min(DATA.energy(1,DATA.clid));
end
DATA = ReplotPCs(DATA,E);
tt(ttn).time = now;
tt(ttn).str = 'Finish';
ttn = ttn+1;
res.times = tt;
end
if spoolspikes
    SpoolSpikes(DATA,DATA.watcharg{:});
end
if DATA.checkclusters
CheckClusters(DataClusters,'Start');
CheckClusters(DataClusters,'CheckNexts','Start');
CheckClusters(DataClusters,'CheckFitSpace');
end
res.toplevel = DATA.toplevel;
if ~isfield(DATA.cluster,'next')
    DATA.cluster.next = {};
end
DATA = CloseLog(DATA);
set(F,'UserData',DATA);

SetGUI(DATA);


function DATA = CloseLog(DATA)
    if DATA.interactive < 0 
        if DATA.logfid > 0
            try 
                fclose(DATA.logfid); 
                fprintf('Closed log (fid%d)\n',DATA.logfid);
            catch
                fprintf('Error Closing log (fid%d)\n',DATA.logfid);
            end
            DATA.logfid = -1;
        end
    end

function C = CheckClusterFields(C)
    if ~isfield(C,'auto')
        C.auto = 0;
    end
    if ~isfield(C,'manual')
        C.manual = 0;
    end
    
function [DATA, V] = ReadSpikeFiles(DATA, name)
    V.name = name;
    id = regexp(V.name,'Expt[0-9]*');
    DATA.exptno = sscanf(V.name(id+4:end),'%d');
    DATA.Expt.exptno = DATA.exptno;
    s = SpkFilename(DATA);
    V.Spikes = ReadSpikeFile(s,'allprobes');
    V.exptno = DATA.exptno;
    [monk, monkey, mdir] = GetMonkeyName(DATA.name);
    if isfield(V.Spikes,'matfile')
        V.matfile = V.Spikes.matfile;
    else
        [a,b] = fileparts(DATA.name);
        V.matfile = [a '/' monk mdir '.mat'];
    end

function str = IDStr(DATA)
str = sprintf('E%dP%d',DATA.exptno,ProbeNumber(DATA));

function C = CheckScoreScaling(DATA, C)
%Template scores are all normalized so that spaces are roughly isotropic
%This can cause a porblem with reapplying when the ratio of
%spikes/nonspikes changes because the SD of the scores changes. So 
%rescale the cluster to compensate for any change in the scaling

return;

if DATA.currentcluster == 1
if isfield(C,'rescaled') && C.rescaled == 1
    return;
end
    if C.space(1) == 3 && isfield(C,'spacescale')
        newscale = DATA.TemplateScaling(C.space(2:end));
        dscale = sqrt(sum(newscale .^2)./sum(C.spacescale.^2));
        C.xyr = C.xyr ./ dscale;
        C.pos = C.pos ./ dscale;
        C.rescaled = 1;
    end
else
    k = DATA.currentcluster -1;
if isfield(C.next{k},'rescaled') && C.next{k}.rescaled == 1
    return;
end
    if C.next{k}.space(1) == 3 && isfield(C.next{k},'spacescale')
        newscale = DATA.TemplateScaling(C.next{k}.space(2:end));
        dscale = sqrt(sum(newscale .^2)./sum(C.next{k}.spacescale.^2));
        C.next{k}.xyr = C.next{k}.xyr ./ dscale;
        C.next{k}.pos = C.next{k}.pos ./ dscale;
        C.next{k}.rescaled = 1;
    end
end
    
    
function ClusterFromPoints(DATA)
    if DATA.cluster.shape == 0  %%
        G = GMfit(DATA.xy{1},2,1,'idlist',DATA.clst);
        id = cluster(G, DATA.xy{1});
        DATA.clst = id;
        if isfield(DATA.cluster,'aspectratio')
            params = FindEnclosingEllipse(DATA.xy{1},DATA.clst,2,DATA.cluster.aspectratio);
        else
            params = FindEnclosingEllipse(DATA.xy{1},DATA.clst,2,1);
        end
        DATA.cluster.xyr = params(1:4);
        DATA.cluster.angle = params(5);
        DATA = ClassifyAll(DATA,1);
    end

function type = WhichPlotType(C,clnum)
    type = 1;
    [a, C] = ClusterIsSet(C,clnum);
    if isfield(C,'space')
      if C.space(1) == 6;
          type = C.space(2);
          if type == 4   %%template
              type = 3;
          end
      elseif C.space(1) == 2 && C.shape == 2 && C.space(2) == 4 %bug in old
          type = C.space(2);
      else
          type = C.space(1);
          if type ==3 && C.space(3) > 12
              type = 4;
          end
      end
    else
        type = 0;
    end
if type == 4
    type = 4; %used to force this to 3. But need 4 if want to plot this automatcially
end


function [DataClusters, FullVData] = LoadDataClusters(DATA)
    DataClusters = {};
    FullVData = [];

    afile = ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);
    if exist(afile,'file')
        load(afile);
       AutoClusters = CondenseClusters(Clusters,0);
    else
        AutoClusters = {};
    end

    cfile = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    if exist(cfile,'file')
        load(cfile);
        DataClusters = Clusters;
        if isfield(DATA,'TemplateScores') %out of date
            DATA = rmfields(DATA,'TemplateScores','TemplateUsed');
        end
        if DATA.plottype == 3 & isfield(Clusters{DATA.probe(1)},'mean') %need to recalculate
            oldscores = 1;
        end
    end


    useauto = zeros(1,DATA.nprobes);
    for j = 1:length(AutoClusters)
        if j > length(DataClusters) ||  ~isfield(DataClusters{j},'mahal')
            DataClusters{j}.auto = 1;
            useauto(j) = 1;
        end
        if isfield(DataClusters{j},'savetime')
            if DataClusters{j}.auto == 1 && AutoClusters{j}.savetime(1) > DataClusters{j}.savetime(1)
                DataClusters{j} = AutoClusters{j};
                DataClusters{j}.auto = 1;
                useauto(j) = 1;
            end
        elseif isfield(DataClusters{j},'mahal')
            DATA = AddErr(DATA,'Probe %d no savetime\n',j);
        else
            useauto(j) = 1;
            if ~isfield(DataClusters{j},'mahal')
                DataClusters{j} = AutoClusters{j};
                DataClusters{j}.auto = 1;
                useauto(j) = 1;
            else
            DataClusters{j}.savetime = [0 0];
            end
        end
    end
    id = find(useauto);
    for j = 1:length(id)
        fprintf('Using AutoCluster for probe %d\n',id(j));
    end
    for j = 1:length(DataClusters)
        if isfield(DataClusters{j},'next') && ~iscell(DataClusters{j}.next) %old style
            last = DataClusters{j}.next;
            DataClusters{j} = rmfield(DataClusters{j},'next');
            DataClusters{j}.next{1} = rmfields(last,'next'); %get rid of next.next
        elseif ~isfield(DataClusters{j},'next') && isfield(DataClusters{j},'mahal')
            DataClusters{j}.next = {};
        end
        if ~isfield(DataClusters{j},'excludetrialids')
            DataClusters{j}.excludetrialids = [];
        end
        if ~isfield(DataClusters{j},'clusterprog')
            DataClusters{j}.clusterprog = '';
        end
        if isfield(DataClusters{j},'next')
        for k = 1:length(DataClusters{j}.next)
            if isfield(DataClusters{j}.next{k},'mahal') && ~isfield(DataClusters{j}.next{k},'chspk') && isfield(DataClusters{j},'chspk')
                DataClusters{j}.next{k}.chspk = DataClusters{j}.chspk;
%                DataClusters{j}.next{k}.quick = 0; %should always be zero when saved
            end
        end
        end
%
% clear any fields that shouldn't be there, but might have been saved by
% earlier versions
        DataClusters{j} = rmfields(DataClusters{j},'clst','r','xy');
%        DataClusters{j}.quick = 0; %should always be zero when saved
%Mar 2012 allowed clsuters toe be save "quick" to allow subsequent
%quantification to save time.
        if ~isfield(DataClusters{j},'marked')
            DataClusters{j}.marked = 0;
        end
    end
    if DATA.checkclusters
    CheckClusters(DataClusters,'Loaded');
    CheckClusters(DataClusters, 'CheckNexts','Loaded');
    CheckClusters(DataClusters,'CheckFitSpace');
    end


function [DATA, Vall, ispk, newdata] = SetupVall(DATA,Vall, ispk, newdata)

if isfield(Vall, 'Spikes')
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
    if length(vt) > size(Vall.V,2)
        PrintMsg(DATA.logfid,'sum blklen > length(V)\n');
        vt = vt(1:size(Vall.V,2));
    end
    Vall.t = vt;
else
    vt= [1:size(Vall.V,2)] .* DATA.interval;
end

forceispk = 1; 
if isfield(DATA,'DataType')
    if strncmp(DATA.DataType,'GridData',8)
        DATA.probelist = 1:96;
    end
end
if size(Vall.V,1) == 1 && ispk(1) > 0
    newdata = 1;
    if ispk > 1
        DATA.probelist = ispk;
        if ispk == Vall.chspk
            ispk = 1;
        elseif forceispk
            Vall.chspk = ispk;
            ispk = 1;
        end
    else
        DATA.probelist = Vall.chspk;
    end
else
    DATA.probelist = 1:DATA.nprobes;
end

if DATA.subtractmeanV
    scale = (meanV*double(Vall.V'))./(meanV*meanV');
    Vall.V = Vall.V - int16(meanV.*scale);
end
pres = {};
if size(Vall.V,1) == 1 && ispk == 0 
    ispk = 1;
    if Vall.chspk < 96
        DATA.probelist = Vall.chspk;
    else
        DATA.probelist = GetProbeFromName(Vall.loadname);
    end
    newdata = 1;
end
if newdata == 1
    DATA.duration = size(Vall.V,2) .* DATA.interval;
end



 function  [true, each] =  NeedTemplateForCluster(C, all)

 needothermeans = 0;
 if ~isfield(C,'space') %%No cluster
     true = 0;
     each = [];
     return;
 end
 each(1) = WhichPlotType(C, 1);
 if sum(ismember(C.space,[17 18]))
     needothermeans = 1;
 end
 if all && isfield(C,'next')
 for j = 1:length(C.next)
         each(j+1) = WhichPlotType(C,j+1);
         if isfield(C.next{j},'space') && sum(ismember(C.next{j}.space,[17 18]))
             needothermeans = 1;
         end
 end
 end
 true = sum(ismember(each,[3 4])) > 0;
 if needothermeans 
     true = 2;
 end
 
 %Find flags/states set in the GUI andcopy them to a new DATA struct
%
function DATA = GetGuiState(DATA, F)
    X = get(F,'UserData');
    if ~isfield(X,'toplevel') || X.toplevel ~= F
        return;
    end
    if X.currentcluster > 1 %&& length(DATA.cluster.next) >= X.currentcluster-1 ... 
            %&& ~isempty(DATA.cluster.next{X.currentcluster-1})
        DATA.currentcluster = X.currentcluster;
    end
    DATA.plotspk.bytrial = X.plotspk.bytrial;
    DATA.plotspk.showmean = X.plotspk.showmean;
    DATA.plotspk.showfullv = X.plotspk.showfullv;
    DATA.plotspk.includeprepost = X.plotspk.includeprepost;
    DATA.elmousept = X.elmousept;
    DATA.logfid = X.logfid;
    DATA.maintitle = X.maintitle;
    DATA.clustericon = X.clustericon;
    DATA.quickcutmode = X.quickcutmode;
    DATA.comparecell = X.comparecell;

    for f = {'handles' 'loadedClusterDetails' 'DataType' 'Expt' 'probeswitchmode'...
            'plotexptfit' 'logname' 'auto' 'checkclusters' 'profiling' 'ptsz'...
            'setspkrate' 'Comments' 'plot' 'LastClusters' 'gui' 'clplot'}
        if isfield(X,f{1})
            DATA.(f{1}) = X.(f{1});
        end
    end
    
    
function DATA = ResetDataForNewProbe(DATA)
    DATA.xy = {[]};
    DATA.clst = [];
    DATA.plotspk.probes = [];  %so that its reset
    DATA.restricttimerange = [];
    DATA.excludetrialids = [];
    DATA.usestdtemplates = 0;
    DATA.trigdt = 0;

function DATA = ReadFromLog(DATA)
    
    logfile = sprintf('%s/ClusterLogExpt%d.mat',DATA.name,DATA.exptno);
    if exist(logfile,'file')
        load(logfile);
        probeids = CellToMat(ClusterLog,'probe');
        exptids = CellToMat(ClusterLog,'exptno');
        reclusters = squeeze(CellToMat(ClusterLog,'recluster'))';
        id = find(probeids == DATA.probe & reclusters(1,:) == 0 & exptids == DATA.exptno )
        cluster = DATA.cluster;
        for j = length(id):-1:1
            C = ClusterLog{id(j)};
            f = {'crit' 'shape' 'space' 'xyr' 'angle',};
            DATA.currentcluster = 1;
            for j = 1:length(f)
                if isfield(C,f{j})
                    cluster.(f{j}) = C.(f{j});
                end
            end
            [a,b] = ClassifySpikes(DATA, cluster);
            [c, d] = Counts(DATA.oldclst(a.id));
            [e,f] = max(c);
            fprintf('%s seems like C%d\n',datestr(C.savetime),d(f));
        end
    end
        
    
function sz = memsize(X)
    x = whos('X');
    sz = x.bytes ./(1024 * 1024 * 1024);
     

function S = SmallCluster(C)
%remove fields from C that use memory
S = rmfields(C,'r','xy');
for j = 1:length(S.next)
    S.next{j} = rmfields(S.next{j},'r','clst');
end


function [true, cluster] = ClusterIsSet(C, cl)
%check that cluster cl is defined for a give clusterstruct
    true = 0;
    cluster = C;
    if cl == 1
        cluster = C;
        if isfield(C,'space')
            true = 1;
        end
    elseif ~isfield(C,'next') && isfield(C,'cluster') && C.cluster == cl
        cluster = C;
        true = 1;
    elseif isfield(C,'next') && length(C.next) >= cl-1 && isfield(C.next{cl-1},'space')
        true = 1;
        cluster = C.next{cl-1};
    end
       

function [DATA, DataClusters, success] = LoadTrigTimes(DATA, checktimes, varargin)
    tic;
    forcecheck = 0;
    clst = {};
    p = DATA.probe(1);
    savexy = 0;
    success = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'savexy',6)
            savexy = 1;
        end
        j = j+1;
    end
    
    fprintf('Loading Trigger times.....');
    if ~isfield(DATA,'trigcheck') || forcecheck
        DATA.trigcheck = zeros(1,DATA.nprobes);
    end
    np = length(DATA.trigcheck);
    if np < DATA.nprobes
        DATA.trigcheck(np+1:DATA.nprobes) = 0;
    end
    if sum(DATA.trigcheck(checktimes) ==0) == 0  && DATA.loadedClusterDetails && isappdata(DATA.toplevel,'ClusterDetails');
        ClusterDetails = getappdata(DATA.toplevel,'ClusterDetails');
        DATA.xy{1} = ClusterDetails{DATA.probe(1)}.xy;
        DATA.clst = ClusterDetails{DATA.probe(1)}.clst;
        fprintf('preloaded\n');
        success = 1;
        return;
    end
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    
    afile = ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);
    %first load details from AutoClusersDetails
    dfile = strrep(afile,'.mat','Details.mat');
    xys = {};
    p = ProbeNumber(DATA);
    if exist(dfile,'file')
        load(dfile);
        for j = 1:length(ClusterDetails)
            if isfield(ClusterDetails{j},'t')
                trigtimes{j} = ClusterDetails{j}.t;
                if isfield(ClusterDetails{j},'clst')
                    clst{j} = ClusterDetails{j}.clst;
                end
                xys{j} = ClusterDetails{j}.xy;
            end
        end
        AutoClusterDetails = ClusterDetails;
    else
        AutoClusterDetails = {};
    end
    cfile = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    dfile = strrep(cfile,'.mat','Details.mat');
    if exist(dfile,'file')
        load(dfile);
        for j = 1:length(AutoClusterDetails)
            if j > length(ClusterDetails) || isempty(ClusterDetails{j}) || ...
                    (DataClusters{j}.auto && AutoClusterDetails{j}.ctime > ClusterDetails{j}.ctime(1))
                ClusterDetails{j} = AutoClusterDetails{j};
            end
        end
        for j = 1:length(ClusterDetails)
            if isfield(ClusterDetails{j},'t')
               xys{j} = ClusterDetails{j}.xy;
               if isfield(ClusterDetails{j},'clst')
                   if diff(size(ClusterDetails{j}.clst)) > 0
                       ClusterDetails{j}.clst = ClusterDetails{j}.clst';
                   end
                   clst{j} = ClusterDetails{j}.clst;
               else
                   clst{j} = ones(length(ClusterDetails{j}.t),1);
               end
               trigtimes{j} = ClusterDetails{j}.t;
               
            end
        end
    end
    if length(clst) >= DATA.probe(1);
        DATA.clst = clst{DATA.probe(1)};
        if isfield(DataClusters{p},'xtimes') & ~isempty(DataClusters{p}.next) & ~isempty(DataClusters{p}.next{1})
            id = find(ismember(trigtimes{p},DataClusters{p}.xtimes{1}));
            DATA.clst(id) = 3;
            DATA.usedxtimes = 1;
        end
    end
    if length(xys) < p %Couldnt load thie probe
        return;
    end
if strncmp(DataClusters{p}.clusterprog,'AllVPcs',7) && isfield(DataClusters{p},'progversion')
    progversion = DataClusters{p}.progversion;
else
    progversion = 100;
end
%Usually if Clustershas an Evec structure, its an error. Only use it if its newer
%But if versions before 1.11 if if forceevec was set, the Evec used for
%calculating was in Clusters, while CluserDetails got new Evec
%lemM239 Expt20 P 3 wants Clusters.Evec even though no field forceevec.
%?why
    if isfield(ClusterDetails{p},'Evec') %% can get Evec in Clusters that is different. ? how? 
        if isfield(DataClusters{p},'Evec')
            if isfield(DataClusters{p}.Evec,'calctime') && DataClusters{p}.calctime > ClusterDetails{p}.calcltime
                DataClusters{p}.Evec = ClusterDetails{p}.Evec;
            elseif (~isfield(DataClusters{p},'forceevec') || DataClusters{p}.forceevec ==0) && progversion < 1.11
                DataClusters{p}.Evec = ClusterDetails{p}.Evec;                
            else
            end
        else
            DataClusters{p}.Evec = ClusterDetails{p}.Evec;
        end
    end
    DATA.xy{1} = xys{DATA.probe(1)};
    if savexy
        setappdata(DATA.toplevel,'ClusterDetails',ClusterDetails);
        DATA.loadedClusterDetails = 1;
    end
    
    Vall = getappdata(DATA.toplevel,'Vall');
    if isempty(checktimes)
        DATA.trigtimes{DATA.probe(1)} = trigtimes{DATA.probe(1)};
    end
%DATA.trigtimes only has times converted to indexes in the FullV array
    for j = checktimes;
        redo = 0;
        if j > length(DATA.trigtimes) || length(DATA.trigtimes{j}) == 0
            redo = 1;
        elseif length(DATA.trigtimes{j}) ~= length(trigtimes{j})
            redo = 1;
            DATA = AddErr(DATA,'Probe %d Trigtime mismatch\n',j);
        end
        if DATA.trigcheck(j) == 0 || forcecheck || redo
            DATA.trigtimes{j} = find(ismember(Vall.t,trigtimes{j}));
            DATA.trigcheck(j) = 1;
            if length(DATA.trigtimes{j}) < length(trigtimes{j})
                id = find(ismember(trigtimes{j},Vall.t(DATA.trigtimes{j})));
                if DATA.probe(1) == j
                DATA.clst = DATA.clst(id);
                DATA.xy{1} = DATA.xy{1}(id,:);
                end
            end
        end
    end
    if ~isempty(DATA.xy{1})
    success = 1;
    end

    fprintf('....Done\n');
toc
    


function [id,  th, nevents] = TriggerV(DATA, rV)
id = [];
th = DATA.Trigger;
for j = 1:size(rV,1)
    sgn(j,:) = diff(sign(diff(rV(j,:),1,2)),1,2);
end
if th(1) < 0
%    id = find(sgn(1,:) > 0)+1;

    id = find(sgn(1,:) > 0 & rV(2:end-1) < th(1))+1;
elseif th(1) > 0
    id = find(sgn(1,:) < 0 & rV(2:end-1) > th(1))+1;
else
    autoth = 1;
end
if isempty(id)
    if DATA.trigdt == 4
        GetFigure(DATA.tag.vhist);
        hold off;
        id = find(abs(sgn) > 1);
        hist(rV(id),500);
    end
    if DATA.thsign == 1 || DATA.trigdt == 3
        id = find(sgn(1,:) < 0)+1;
        prc = DATA.setnspk .* 100./length(id); % get 1000 spikes
        th(1) = prctile(rV(id),100-prc);
        id = id(rV(1,id) > th(1));
    else
        id = find(sgn(1,:) > 0)+1;
        prc = DATA.setnspk .* 100./length(id); % get 1000 spikes
        if prc > 100 %can happen if nspk > # minima
            th(1) = max(rV(1,id));
        else
            th(1) = prctile(rV(1,id),prc);
        end
        id = id(rV(1,id) < th(1));
        if size(rV,1) > 1
            xid = find(sgn(2,:) > 0)+1;
            prc = DATA.setnspk .* 100./length(xid); % get 1000 spikes
            th(1) = prctile(rV(2,xid),prc);
            xid = xid(rV(2,xid) < th(1));
            id = union(xid,id);
        end
    end
end

%if the ISI is very  short, and the trigger channel does not go back to
%near zero (th/3) between two trigger points, then throw away the one
%witht the smaller triggger
nevents = length(id);
sid = find(diff(id) < DATA.isicheck(1));
if length(sid) < length(rV) .* 0.01 %must be a low trigger
okid = [];
for j = 1:length(sid)
    if min(rV(id(sid(j)):id(sid(j)+1)) .*sign(th(1))) < abs(th(1))/DATA.isicheck(2)
        okid = [okid sid(j)];
    end
end
sid = setdiff(sid,okid);
v = cat(1,rV(id(sid)),rV(id(sid+1)));
if th(1) > 0
[a,b] = min(v);
else
[a,b] = max(v);
end    
xid = id(sid+b-1);
fprintf('%s Removing %d double Triggers\n',IDStr(DATA),length(xid));
yid = id(sid+2-b);
id = setdiff(id,xid);
else
    DATA = AddErr(DATA,'Too many double Triggers (%d/%d)\n',length(sid),length(rV));
end
id = id(id > -DATA.spts(1) & id < size(rV,2)-DATA.spts(end));

if DATA.trigdt == 3  %spk energy trigger
    FullV = getappdata(DATA.toplevel,'Vall');
    for j =1:length(id)
        if DATA.thsign == 1
        [a,b] = max(FullV.V(DATA.probe(1),id(j)-5:id(j)+5));
        else
        [a,b] = min(FullV.V(DATA.probe(1),id(j)-5:id(j)+5));
        end
        id(j) = id(j)+b-7;
    end
end

function res = AutoCutAll(ispk, toplevel, Vall,DATA, args)
    quickfollow = 0;
    maxthriter = 1;
    nolog = 0;
    j = 1;
    while j <= length(args)
        if strncmpi(args{j},'quickfollow',8)
            quickfollow = 1;
        elseif strncmpi(args{j},'nolog',5)
            nolog = 1;
        end
        j = j+1;
    end
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    autofit = DATA.autofit;
    if length(ispk) < 3 && 0 %what was this for
        ispk = 1:size(Vall.V,1);
    end
    tstart = now;
    allcuts = {};
    if nolog == 0
    logfile = ClusterFile(Vall.name,'log','subdir',DATA.clustersubdir,DATA.Expt);
    logfid = fopen(logfile,'a');
    else
        logfid = -1;
    end
    if isfield(Vall,'exptno')
        name = sprintf('%s Expt %d',Vall.name,Vall.exptno);
    else
        name = Vall.name;
    end
    PrintMsg(logfid,'Start on %s at %s\r\n',name,datestr(now));
    if size(Vall.V,1) == 1 && strncmp(DATA.DataType,'Grid',4)
        exname = Vall.loadname;
        for j = 1:length(ispk)
            istart = now;
            name = regexprep(exname, '\.p[0-9]*',sprintf('.p%d',ispk(j)));
            if strcmp(DATA.autocutmode,'quick')
                res.cls{j} =  AllVPcs(name ,'nowatch',args{:},'tchan',ispk(j),'logfid',logfid,'quickauto','toplevel',toplevel);
            else
                res.cls{j} =  AllVPcs(name ,'nowatch',args{:},'tchan',ispk(j),'logfid',logfid,'saveautocut','toplevel',toplevel);
            end
            res.cls{j} = rmfields(res.cls{j},'pcs','t');
        end
        DATA = get(toplevel,'UserData');
        DATA.logfid = -1;
        res.Clusters = DataClusters; 
        return;
    end
    if size(Vall.V,1) < max(ispk)
        PrintMsg(logfid,'Only %d/%d probes',size(Vall.V,1),max(ispk));
    end
    ispk = ispk(find(ispk <= size(Vall.V,1)));
    for j = ispk
        istart = now;
        a =  AllVPcs(Vall ,'nowatch',args{:},'tchan',j,'logfid',logfid,'saveautocut');
        t = a.cluster.Trigger;
        spksd = std(a.cluster.MeanSpike.ms,0,2);
        musd = std(a.cluster.MeanSpike.mu,0,2);
        [sv, svid] = max(spksd(a.chspk));
           a.cluster.good = GoodCluster(a.cluster);
           a.maxmean = a.chspk(svid);
           allcuts = {allcuts{:} a};
        otherps = setdiff(a.chspk, j);
%
%  Explore lowering the trigger level if necessary. This can be for two
%  reasons:
%            1)  Dropping spikes. Produces a clipped distribution of
%            trigger point voltages, measured with dropi
%            2) cluster not dtopping too many spikes, but there are very
%            few non-cluster events, so the real cell is divided into tw
%            similar groups. NeedMore() Checks for this
%only explore lower triggers if the spike is biggest on this channel. Otherwise
%it will go to very low triggers to get all spikes that are really on another
%channel, and run out of memory
           evi = a.Evec.Eval(1)./sum(a.Evec.Eval);
           nloop = 0;
           while NeedMore(a.cluster, a.Evec) && a.chspk(svid) == j && a.cluster.nspks < 500000 && nloop < autofit.maxthriter
               nloop = nloop+1;;
               a.cluster.needmore = NeedMore(a.cluster,a.Evec);
               PrintMsg(logfid, 'P%d: NeedMore%d %.2f(%.2f), bmi %.3f mahal%.2f, increasing events to %d\n',...
                   j,a.cluster.needmore,a.cluster.dropi(3),a.cluster.trigsd,a.cluster.bmc,a.cluster.mahal(1),a.cluster.nspks*2);
                 a = AllVPcs(Vall, 'nowatch',args{:},'tchan',j, 'spkrate',a.spkrate * 2,'logfid', logfid, 'previous',PrevCluster(a.cluster),'saveautocut');
                 evi = a.Evec.Eval(1)./sum(a.Evec.Eval);
           end
%if big spike is on another probe, can get spike/MU reversed, makeing it look liek the spike is biggest on 
%this probe. So, check that size of spike on this probe is bigger that mu
%on other probes
          if length(a.cluster.chspk) > 1
              peakhere = spksd(j) > max(spksd)/1.1 && spksd(j) > max(musd(otherps))/2;
          else
              peakhere = 1;
          end
           if GoodCluster(a.cluster)  && peakhere
               dropi = a.cluster.dropi(3);
               nloop = 0;
           while ((a.cluster.dropi(3) < 2 && t < -a.xsd && GoodCluster(a.cluster) > 1 ...
                   && a.cluster.dropi(3) >= dropi)) && ...
                    a.maxspksused == 0 && nloop < autofit.maxthriter
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
               PrintMsg(logfid,'P%d: DropSD %.2f(%.2f), bmi %.3f mahal%.2f, lowering threshold to %.3f, mine %.2f, took %.1f size %.1f,%.1f\n',...
                   j,a.cluster.dropi(3),a.cluster.trigsd,a.cluster.bmc,a.cluster.mahal(1),t,a.cluster.minspke,mytoc(tstart),a.memsz(1),a.memsz(2));
               a.cluster.needmore = 4;
%put last cluster in a.cluster.first, so if pass a.cluster to AllVPcs, it is stored               
               a.cluster.first = PrevCluster(a.cluster);
               if isnan(a.cluster.minspke)
                   a.cluster.minspke = 0;
               end
%If following up on previous cut, its quicekst to use the last as a
%template and just use template space
%Its quicker, stops it chasing different clusters each time
               if (a.cluster.space(1) == 6 && a.cluster.space(2) == 4) || quickfollow || autofit.maxthriter == 1
               a = AllVPcs(Vall, 'nowatch',args{:},'tchan',j, 'thr',t,'maxrate',a.spkrate * 2, 'mine',a.cluster.minspke,'minvar',a.cluster.minspkvar,'reapply',a.cluster,'logfid',logfid,'saveclusters');
               else
               a = AllVPcs(Vall, 'nowatch',args{:},'tchan',j, 'thr',t,'mine',a.cluster.minspke,'minvar',a.cluster.minspkvar,'logfid',logfid,'previous',a.cluster.first,'saveautocut');
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
           PrintMsg(logfid,'P%d(%d): %d/%d Spikes, dropi %.2f (%.2f) bmc %.3f, mahal %.1f,%.1f,%.1f G%d took %.2f size %.1f, %.1f at %s',...
               j,a.maxmean,c.ncut,c.nspks,c.dropi(3),c.dropi(4),c.bmc,c.mahal(1),c.mahal(2),c.mahal(4),GoodCluster(a.cluster),...
               mytoc(istart),a.memsz(1),a.memsz(2),datestr(now));
           DataClusters{j}.totaltime = mytoc(istart);
              
    end
        if logfid > 2
            fclose(logfid);
        end
        DATA = get(toplevel,'UserData');
        if DATA.plot.xcorr
        a = PlotSpikeTimes(DataClusters,'xcorr');
        for j = 1:length(DataClusters)
%            DataClusters{j}.synci = a.synci(j,:);
        end
        end
        DATA.logfid = -1;
%        SaveClusters(DATA,ClusterFile(DATA.name,DATA.Expt,'auto'));
        res.Clusters = DataClusters;
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
    first.Trigger = C.Trigger;

function res = QuantifyQuickClusters(DATA, ispk, varargin) 
    
    passonargs = {};
    j = 1;
    while j <= length(varargin)
        passonargs = {passonargs{:} varargin{j}};
        j = j+1;
    end
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    for j = 1:length(ispk)
            if isfield(DataClusters{ispk(j)},'quick') && DataClusters{ispk(j)}.quick
                fprintf('Reclassifying E%d P%d\n',DATA.exptno,ispk(j));
                res{j} = AllVPcs(DATA.toplevel, 'tchannew', ispk(j), passonargs{:},'reclassify');
                res{j}.exptno = DATA.exptno;
                res{j}.probes = ispk(j);
            end
    end

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
                DATA.nid = setdiff(1:DATA.nevents,DATA.clid)';
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
    AllV = GetAllV(DATA);
    
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
    space = DATA.cluster.space;
    oldms = DATA.cluster.MeanSpike.ms;
    if size(DATA.cluster.gmfit.mu,2) == 6
        DATA.cluster.space = [6 4];
    end
    if fittemplate
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
        [a,b, DATA.xy{1}, details] = TemplateSpace(DATA,'template');
        DATA.ndxy = DATA.xy{1};
        fprintf('Template GM %.2f 6D %.2f 2d, %.2f1D\n',a(1),a(2),details.gmdprime)
        DATA.cluster.gmfit = b;
        DATA.cluster.gmdprime = details.gmdprime;
        DATA.cluster.sign = details.cluster.sign;
        DATA.cluster.crit = details.cluster.crit;
        DATA.cluster.shape = 2;
        DATA.cluster.angle = 0;
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.xy = DATA.xy{1}; %best ND;
        DATA.cluster.xy(:,3) = DATA.TemplateScores(:,2); %sum of r
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.bestcl = details.bestcl;
        DATA.cluster.usegmcluster = details.usegmcluster;

    elseif ismember(DATA.cluster.space(1),[3 4]) & ~isfield(DATA,'TemplateScores')
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
        a = DATA.cluster,space
        DATA.xy{1}(:,1) = DATA.TemplateScores(space(2),:);
        DATA.xy{1}(:,2) = DATA.TemplateScores(space(3),:);
    elseif DATA.cluster.space(1) == 6 && DATA.cluster.space(2) == 4
        DATA =  CalcTemplatesFromMean(DATA,DATA.cluster.MeanSpike);
        DATA.plottype = 3;
            [a,b, DATA.xy{1}, details] = BestSpace(DATA,'template');
            DATA.ndxy = DATA.xy{1};
            DATA.clst = details.cid;
            DATA.cluster.gmfit = b;
        DATA = ReplotPCs(DATA,[]);
    elseif DATA.cluster.space(1) == 6
        DATA.ndxy = ProjectND(DATA,DATA.cluster.space(2),DATA.cluster.gmfit);
        if space(2) == 3 %ADC
            DATA.plottype = 2;
        elseif space(2) == 4 %template
            DATA.plottype = 3;
        else
            DATA.plottype = 1;
        end
    elseif space(1) == 1
        DATA.xy{1}(:,1) = DATA.pcs(DATA.uid,space(2));
        DATA.xy{1}(:,2) = DATA.pcs(DATA.uid,space(3));
        DATA.plottype = 1;
    elseif space(1) == 3 || (space(1) ==6 && space(2) == 3)
        DATA.xy{1}(:,1) = AllV(DATA.probe(1),space(2), DATA.uid);
        DATA.xy{1}(:,2) = AllV(DATA.probe(1),space(3), DATA.uid);
        DATA.plottype = 2;
    end
        if newbound
%            [dip, details] = FindDip(DATA.xy{1}(:,1),DATA.energy(1,:),'gmix');
%need to apply rotation to 2-D spaces first
            [dip, details] = GMDip(DATA.xy{1},DATA.energy(1,:),'label',DATA.idstr);
            DATA.cluster.crit = dip(1);
            DATA.cluster.sign = details.sign;
            DATA.cluster.gmdipres = details.dipres;
            DATA.cluster.gmdprime = details.gmdprime;
            DATA.cluster.autodipsize = details.dipsize;
            DATA.cluster.dipsize = details.cdipsize;
        end
    [cl, cluster, DATA.xy{1}] = ClassifySpikes(DATA,DATA.cluster);
    if newbound || fittemplate%if boundary changed, need to record new criterion, threshold/dropi etc.
        DATA.cluster = cluster;
    end
    DATA.Cluster{DATA.probe(1)}.mahal = cluster.mahal;
    DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(BoundaryFromCluster([],cluster, DATA.currentcluster));
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    DATA.MeanSpike = cl.MeanSpike;
    DATA.cluster.MeanSpike = cl.MeanSpike;
    DATA.cluster.minspke = prctile(DATA.energy(1,DATA.clid),1) .* 0.95;
    DATA.cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),DATA.clid),1) .* 0.95;
    if fittemplate
        DATA.cluster.xy = DATA.xy{1};
        DATA.cluster.xy(:,3) = DATA.TemplateScores(:,2);
        xc = corrcoef(oldms(:),DATA.cluster.MeanSpike.ms(:));
        DATA.cluster.templatexc = xc(1,2);
    end

function DATA = CalcTemplateScores(DATA,  varargin)
    if  isfield(DATA.cluster,'TemplateUsed')
        MeanSpike.ms = DATA.cluster.TemplateUsed;
        if size(DATA.cluster.mumeanUsed,1) == length(DATA.cluster.chspk)
        MeanSpike.mu(DATA.cluster.chspk,:) = DATA.cluster.mumeanUsed;
        else
            MeanSpike.mu = DATA.cluster.mumeanUsed;
        end
        MeanSpike.vdprime(DATA.cluster.chspk,:) = DATA.cluster.DprimeUsed;
    else
        MeanSpike = DATA.cluster.MeanSpike;
    end
    DATA = CalcTemplatesFromMean(DATA, MeanSpike, varargin{:});



function DATA = CalcTemplatesFromMean(DATA, MeanSpike, varargin);
        
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'stdtemplate',6)
            DATA.usestdtemplates = 1;
        end
        j = j+1;
    end
    cl = DATA.currentcluster;
    DATA.TemplateScores = zeros(DATA.nevents, 12);
    DATA.TemplateLabels = TemplateLabels(DATA,0);
    [Scores, DATA.TemplateUsed, a] = CalcScores(DATA,MeanSpike);
    if DATA.profiling
        fprintf('Return from Calcscores took %.4f\n',mytoc(a.endtime));
    end
    chspk = DATA.chspk(DATA.chspk <= size(MeanSpike.vdprime,1));
    DATA.DprimeUsed = MeanSpike.vdprime(chspk,:);
    if DATA.currentcluster > 1
        DATA.cluster.next{cl-1}.mumeanUsed = MeanSpike.mu(chspk,:);
        DATA.cluster.next{cl-1}.DprimeUsed = MeanSpike.vdprime(chspk,:);
    else
        DATA.cluster.mumeanUsed = MeanSpike.mu(chspk,:);
        DATA.cluster.DprimeUsed = MeanSpike.vdprime(chspk,:);
    end
    ispk = find(DATA.chspk == DATA.probe(1));

    DATA.tmpdips = zeros(1,8);
    if isfield(DATA,'TemplateScores')
    oldscores = DATA.TemplateScores;
    end
    if DATA.usestdtemplates
        DATA.TemplateScores = squeeze(Scores(1,:,:))';
        return;
    end
    if size(Scores,1) >= ispk
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
    DATA.TemplateScores(:,13)= sum(Scores(:,7,:)); %absdiff
    DATA.TemplateScores(:,14)= sum(Scores(:,6,:)); %mu score
    DATA.TemplateScores(:,15)= sum(Scores(:,8,:)); %C2 score
    DATA.TemplateScores(:,16)= sum(Scores(:,8,:))-sum(Scores(:,2,:)); %P-1 score
    
    if length(DATA.chspk) == 1
        DATA.TemplateScores(:,2) = Scores(1,6,:); %MU r
        DATA.TemplateScores(:,3) = Scores(1,3,:); %std1
        DATA.TemplateScores(:,4) = Scores(1,4,:); %std1
        DATA.TemplateScores(:,11) = Scores(1,7,:); %Trigger
        DATA.TemplateScores(:,10)= Scores(1,5,:); %dprime
        DATA.TemplateScores(:,12)= Scores(1,8,:); %std1 dvdt
    end
    if size(Scores,2) > 9
    DATA.TemplateScores(:,17)= sum(Scores(:,9,:)); %P-1 score
    DATA.TemplateScores(:,18)= sum(Scores(:,10,:)); %P+1 score
    end
    if DATA.tmpnorm
        stds = std(DATA.TemplateScores);
        for j = 1:size(DATA.TemplateScores,2)
            DATA.TemplateScores(:,j) = DATA.TemplateScores(:,j)./stds(j);
            DATA.TemplateScaling(j) = stds(j);
        end
    end
%    DATA.tmpdips = CalculateTemplateDips(DATA);

    
    
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
    f = fields(E);
    for j =  1:length(f)
        C.(f{j}) = E.(f{j});
    end

    
    if isfield(E,'angle') && E.pos(1) == E.pos(3)
        C.angle = E.angle;
        C.crit =E.pos(1);
    else
        C.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),C.angle);
        C.crit = mean(exy(:,1));
    end
    C.len = abs(diff(E.pos([1 3])) + i * diff(E.pos([2 4])));
    C.ctype = 0;
        
function E = BoundaryFromCluster(E, C, n)
    
    if n  > 1 && length(C.next) >= n-1
        C = C.next{n-1};
    end

f = fields(C);
for j = 1:length(f)
        E.(f{j}) = C.(f{j});
end

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
if ~isfield(C,'shape')
    E.shape = 1;
    E.ctype = 1;
    return;
end
if C.shape == 2 %line in ND space
    if ~isfield(C,'angle')
        C.angle = 0;
    end
E.xyr(1) = C.crit(1) .* cos(-C.angle);
E.xyr(2) = C.crit(1) .* sin(-C.angle);
exy = xyrotate([C.crit(1) C.crit(1)],y,-C.angle);
E.pos(1) = exy(1,1);
E.pos(2) = exy(1,2);
E.pos(3) = exy(2,1);
E.pos(4) = exy(2,2);
elseif C.shape == 1 %line in
    if ~isfield(C,'crit')
        exy = xyrotate(C.pos([1 3]),C.pos([2 4]),C.angle);
        C.crit = exy(1,1);
    end
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
    E.xyr = C.xyr;
else
    E.pos = C.pos;
end


angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));

if isfield(C,'sign')
    E.sign = C.sign;
else
    E.sign = 0;
end
if isfield(C,'shape')
    E.shape = C.shape;
    if C.shape == 2
        E.space = C.space;
    elseif C.shape == 1 %% not sure if this needs difffernt treatment
        E.space = C.space;
    else 
        E.space = C.space;
    end
end
if ~isfield(E,'h')
    E.h = [];
end
if ~isfield(E,'shape')
    E.shape = 1;
end
if isfield(C,'next')
    for j = 1:length(C.next)
        if isfield(C.next{j},'shape') %not empyt or just meanspike
            C.next{j} = rmfields(C.next{j},'next'); %in case old/error
        E.next{j} = BoundaryFromCluster([],C,1+j);
        else
        E.next{j} = [];
        end
    end
end
E.ctype = 1;
E.endtime = now;


function DATA = NextPCs(DATA)
   
    DATA.hidecluster = 2;
    DATA.currentcluster = 2;
    AllV = GetAllV(DATA);
    cid = find(DATA.clst == DATA.hidecluster);
    DATA.uid = find(DATA.clst ~= DATA.hidecluster & DATA.clst > 0);
    [C, DATA.Evec, pcs, dip, chspk, errs, DATA.pcfit] = CalcPCs(DATA,AllV,DATA.nprobepc);
    DATA.pcs(DATA.uid,:) = pcs;
    id = setdiff(1:length(DATA.clst),DATA.uid);
    DATA.clst(cid) = -DATA.hidecluster;
    DATA = ReplotPCs(DATA,[]);
    
function DATA = SetPCs(DATA, replot)
    AllV = GetAllV(DATA);
    [C, DATA.Evec, pcs, dip, chspk, errs, DATA.pcfit] = CalcPCs(DATA, AllV,DATA.nprobepc);
    DATA.errs = {DATA.errs{:} errs{:}};
    DATA.pcs(DATA.uid,1:size(pcs,2)) = pcs;
    if replot
    DATA = ReplotPCs(DATA,[]);
    end

function DATA = CalcICA(DATA,ns)    
    
AllV = GetAllV(DATA);
chspk = DATA.chspk;
uid = DATA.uid;
if DATA.csd 
    if DATA.csd == 2
    csd = diff(AllV,2,1);
    else
    csd = diff(AllV,1,1);
    end
    TV = csd(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,csd(chspk(j),:,uid));
    end
elseif DATA.dvdt == 2
    TV = AllV(chspk(1),:,uid);
    TV = cat(2,TV,diff(AllV(chspk(1),:,uid),1,2));
    for j = 2:length(chspk)
        TV = cat(2,TV,AllV(chspk(j),:,uid));
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),1,2));
    end
elseif DATA.dvdt == 3
    TV = diff(AllV(chspk(1),:,uid),2,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),2,2));
    end
elseif DATA.dvdt
    TV = diff(AllV(chspk(1),:,uid),1,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),1,2));
    end
else
    TV = AllV(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,AllV(chspk(j),:,uid));
    end
end
TV = squeeze(TV)';

if diff(size(TV)) > 0
DATA.icas = FastICA(TV,'numofic',ns);
else
DATA.icas = FastICA(TV','numofic',ns)';
end
if ns < 5 %fill with zero for plotting
    DATA.icas(end,5) = 0;
end



 function [C, Evec, pcs, dip, chspk, errs, details] = CalcPCs(DATA, AllV, nprobepc,  varargin)

errs = {};
details = [];
uid = DATA.uid;
Evec = [];
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'Evec')
        Evec = varargin{j};
    end
    j = j+1;
end

if length(nprobepc) > 1 || nprobepc(1) >= 0
    if length(nprobepc) > 1
        chspk = DATA.probe(1) + nprobepc;
    else
        chspk = DATA.probe(1)-nprobepc:1:DATA.probe(1)+nprobepc;
    end
    if DATA.csd == 2
        chspk = chspk-2;
        chspk = chspk(chspk > 0 & chspk < size(AllV,1)-1);
        if isempty(chspk)
            chspk = DATA.probe(1);
        end
    elseif DATA.csd
        chspk = chspk-1;
        chspk = chspk(chspk > 0 & chspk < size(AllV,1));
    else
        chspk = chspk(chspk > 0 & chspk <= size(AllV,1));
    end
elseif nprobepc <0 
    chspk = DATA.chspk;
end

evspace = 0;
if DATA.csd 
    if DATA.csd == 2
        csd = diff(AllV,2,1);
        evspace = 1;
    else
        csd = diff(AllV,1,1);
        evspace = 2;
    end
    TV = csd(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,csd(chspk(j),:,uid));
    end
elseif DATA.dvdt == 2
    TV = AllV(chspk(1),:,uid);
    TV = cat(2,TV,diff(AllV(chspk(1),:,uid),1,2));
    for j = 2:length(chspk)
        TV = cat(2,TV,AllV(chspk(j),:,uid));
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),1,2));
    end
    evspace = 3;

elseif DATA.dvdt == 3
    TV = diff(AllV(chspk(1),:,uid),2,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),2,2));
    end
    evspace = 4;
elseif DATA.dvdt
    TV = diff(AllV(chspk(1),:,uid),1,2);
    for j = 2:length(chspk)
        TV = cat(2,TV,diff(AllV(chspk(j),:,uid),1,2));
    end
    evspace = 5;
else
    evspace = 6;
    TV = AllV(chspk(1),:,uid);
    for j = 2:length(chspk)
        TV = cat(2,TV,AllV(chspk(j),:,uid));
    end
end
TV = squeeze(TV)';

if isempty(Evec) || size(TV,2) ~= size(Evec.Evec,1) 
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
%arbitrary sign convention just so that its consistent when re-applying    
    for j = 1:size(pc,1)
        if sum(pc(:,j)) < 0
            pc(:,j) = -pc(:,j);
        end
    end
    Evec.Evec = pc;
    Evec.calctime = now;
    Evec.chspk = chspk;
    Evec.npts = length(uid);
    Evec.space = evspace;
    Evec.npts = length(uid);
else
    Evec.forced = 1;
    Evec.newspace = evspace;
    C = [];
end
pcs = TV*Evec.Evec;

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
    [P, dip(1), details] = GMfit(pcs(:,1:4),2,1);
    end
end


function OldPlotFullV(V, ispk, DATA)

SetFigure('FullV');
AllV = GetAllV(DATA);
hold off;
y = 0;
for p = 1:length(ispk)
plot(V(ispk(p),:)+y);
hold on;
id = DATA.clid;
spts = DATA.spts;
for j = 1:length(id);
    plot([DATA.t(id(j)):DATA.t(id(j))+length(spts)-1]+spts(1),AllV(ispk(p),:,id(j))+y,'r');
end
nid = DATA.nid;
for j = 1:length(nid);
    plot([DATA.t(nid(j)):DATA.t(nid(j))+length(spts)-1]+spts(1),AllV(ispk(p),:,nid(j))+y,'g');
end
y = y + max(V(ispk(p),:));
end
if isfield(DATA,'onespiketime')
    for j = 1:length(DATA.onespiketime)
        plot(DATA.onspiketime,get(gca,'ylim'),'b:');
    end
end

hold off;


function distance  = gmdistance(G)

    try
        D = mahal(G,G.mu);
        distance = sqrt(2./((1./D(1,2))+(1./D(2,1))));
    catch %in case fit failed
        distance = 0;
    end

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
  
  
function [theta, c, details] = BestAngleGM(xy, G, dipfit, varargin)
%
% Find best angle to cut in a 2D space, using 1-D Gaussian Mixtures

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
      G = GMfit(xy,2,1);
  end
  id = cluster(G, xy);
  if isempty(dipfit)
      [dip, dipfit] = GMDip(xy,[],'idlist',id,'noplot');
  end
  if quickmode
      cid = find(id > 1);
      nid = find(id == 1);
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
          [aa,bb] = GMDip(XY,[],'idlist',id,'noplot');
          details.d(j) = bb.mahal(bb.best);
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
      [aa,bs{2}] = GMDip(XY,[],'idlist',id);
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
%      [aa,bb] = FindDip(xy(:,1),DATA.energy(1,:));
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

 function CheckClusters(Clusters, str, varargin)

     if isempty(Clusters)
         return;
     end
     if strcmp(str,'CheckFitSpace')
         for j = 1:length(Clusters)
             if iscell(Clusters)
                 C = Clusters{j};
             else
                 C = Clusters(j);
             end
             if isfield(C,'gmfit') && isfield(C.gmfit,'mu') %might be empty 
               
             if C.space(1) == 6
                 fitsz = size(C.gmfit.mu,2);
                 spacesize = [4 2 4 6 6 6 4];
                 ssz = spacesize(C.space(2));
                 if fitsz ~= ssz;
                     if C.auto
                         astr = 'auto';
                     else
                         astr = '';
                     end
                     if length(Clusters) > 1
                         cstr = sprintf('%d',j);
                     else
                         cstr = sprnintf('P%d',C.probe);
                     end
                     fprintf('Cluster %s Space %d, %d dims, but fit %d dims %s\n',cstr,C.space(2),ssz,fitsz,astr);
                 end
             end
             end
         end
         return;
     end

     if strcmp(str,'CheckNexts')
         nexterrs = 0;
         for j = 1:length(Clusters)
             if iscell(Clusters)
                 C = Clusters{j};
             else
                 C = Clusters(j);
             end
             
             if isfield(C,'next') %might be empty 
                 wasempty = 0;
                 for k = 1:length(C.next)
                     if isempty(C.next{k})
                         wasempty = 1;
                     elseif wasempty
                         wasempty = 2;
                     end
                 end
                 if wasempty == 1 && k > 1 %trailing empty next
                     fprintf('Expt %d Probe %d Cluster missing next{}\n',C.exptno,j);
                     nexterrs = nexterrs+1;
                 end
             end
         end
         %allowed empty nexts now.
         if nexterrs
%                 errordlg(sprintf('Missing Nexts in %d Clusters at %s',nexterrs,varargin{1}),'Cluster Load error','modal');
         end
         return;
     end

     for j = 1:length(Clusters)
         C = Clusters{j};
         if ~isfield(Clusters{j},'mahal') %empty really
             missing(j) = 1;
             fprintf('Missing Cluster on %d\n',j);
         else
             missing(j) = 0;
         end
         if ~isfield(Clusters{j},'MeanSpike')
             fprintf('Missing Meanspike on %d\n',j);
             missing(j) = 2;
         end
         if isfield(Clusters{j},'restricttimerange') && length(C.restricttimerange) == 2
             if min(C.times) > C.restricttimerange(2) || max(C.times) < C.restricttimerange(1)
                 errordlg(sprintf('Bad Time range for probe %d',j),'restrict time error');
             end
         end
     end
 nerr = sum(missing);
 if nerr
     errordlg(sprintf('Missing %d Clusters at %s',nerr,str),'Cluster Load error','modal');
 end
 nm = sum(missing == 2);
 if nm
     errordlg(sprintf('Missing %d MeanSpikes at %s',nm,str),'Cluster Load error','modal');
 end
            

function SaveMeanSpikeOnly(DATA, outname) 
    
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    if DATA.checkclusters
        CheckClusters(DataClusters,'Save');
    end
%Saving should only change one probe. So read in fiel from disk, and
%only update the current probe, then write out. 
    if exist(outname,'file')
        load(outname);
        for j = 1:length(DataClusters);
            if j > length(Clusters) || isempty(Clusters{j})
                Clusters{j} = DataClusters{j};
            elseif ~isfield(Clusters{j},'mahal') %empty really
                Clusters{j} = DataClusters{j};
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
            if ~isfield(Clusters{j},'auto')
                Clusters{j}.auto = DataClusters{j}.auto;
            elseif strfind(outname,'AutoClusterTimes')
                Clusters{j}.auto = 1;
            end
            if ~isfield(Clusters{j},'mahal') %still empty
                AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
                if length(AutoClusters) >= j
                    Clusters{j} = AutoClusters{j};
                    fprintf('Cluster %d was empty - reverting to Auto\n',j);
                    errordlg(sprintf('Cluster %d was empty Reloaded AutoCluster\n',j),'Cluster Error','modal');                    
                else
                    errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',j),'Cluster Error','modal');                    
                end
                    
            end
        end
    else
        Clusters = DataClusters;
    end

    Clusters{p}.MeanSpike = DATA.cluster.MeanSpike;
    for j = 1:length(DATA.cluster.next)
        if isfield(DATA.cluster.next{j},'MeanSpike')
            Clusters{p}.next{j}.MeanSpike = DATA.cluster.next{j}.MeanSpike;
        end
    end
  
    
function CheckClusterValues(DATA, C)
    if DATA.check.dropi(1) > 0 
        if C.dropi(3) < DATA.check.dropi(1) && -C.fitdprime(1) > DATA.check.dropi(2)
            msg = sprintf('Cluster %d (GM%.2f) May be dropping spikes (dropi %.2f)',C.cluster,C.fitdprime(1),C.dropi(3));
            if DATA.interactive > 0
                errordlg(msg,'Cluster Warning','modal');
            end
            PrintMsg(DATA.logfid,msg);
        end
    end

function [DATA, id] = SaveClusters(DATA, outname,varargin)

    quickmode = 0;
    savexy = 2;

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'quick',5)
            quickmode = 1;
            savexy = 0;
        end
        j = j+1;
    end
    
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    dname = strrep(outname,'.mat','Details.mat');
    logname = [DATA.datadir '/' 'ClusterLogExpt' num2str(DATA.exptno) '.mat'];
    if DATA.toplevel
    set(DATA.toplevel,'name',sprintf('Saving %s',outname));
    drawnow;
    Vall = getappdata(DATA.toplevel,'Vall');
    else
        Vall = DATA.Vall;
    end
    if ~isfield(DATA.cluster,'next') || ~iscell(DATA.cluster.next) %should not need this - check
        DATA = AddErr(DATA,'ERROR!!!!!!!!!!! - missing .next\n');
        DATA.cluster.next  = {};
    end
    if quickmode == 0
          DATA = ClassifyAll(DATA,0); %Check all are up to date
          if DATA.toplevel; set(0,'currentfigure',DATA.toplevel); end
    end
    
    if sum(strcmp(DATA.DataType,{'Spike2'  'Default'}));
        SaveArrayConfig(DATA);
    end
    if DATA.checkclusters
    CheckClusters(DataClusters,'Save');
    CheckClusters(DataClusters,'CheckNexts','Save');
    end
%Saving should only change one probe. So read in file from disk, and
%only update the current probe, then write out. 
    if exist(outname,'file')
        load(outname);
        for j = 1:length(DataClusters);
            if j > length(Clusters) || isempty(Clusters{j})
                Clusters{j} = DataClusters{j};
            elseif ~isfield(Clusters{j},'mahal') %empty really
                Clusters{j} = DataClusters{j};
            elseif isfield(Clusters{j},'mean') %old
                Clusters{j} = rmfield(Clusters{j},'mean');
            elseif isfield(Clusters{j},'r') 
                Clusters{j} = rmfield(Clusters{j},'r');
            elseif isfield(Clusters{j},'clst') 
                Clusters{j} = rmfield(Clusters{j},'clst');
            elseif isfield(Clusters{j},'t') 
                Clusters{j} = rmfield(Clusters{j},'t');
            end
            
            if isfield(Clusters{j},'next') && ~iscell(Clusters{j}.next) %old style
                last = Clusters{j}.next;
                Clusters{j} = rmfield(Clusters{j},'next');
                Clusters{j}.next{1} = rmfields(last,'next'); %get rid of next.next
            elseif ~isfield(Clusters{j},'next') && isfield(Clusters{j},'mahal')
                Clusters{j}.next = {};
            end
            
            
            if ~isfield(Clusters{j},'probe')
                Clusters{j}.probe = j;
            end
            if ~isfield(Clusters{j},'auto')
                if isfield(DataClusters{j},'auto')
                    Clusters{j}.auto = DataClusters{j}.auto;
                else
                    Clusters{j}.auto = 0;
                end
                    
            elseif strfind(outname,'AutoClusterTimes')
                Clusters{j}.auto = 1;
            end
            if ~isfield(Clusters{j},'mahal') && isappdata(DATA.toplevel,'AutoClusters') %still empty
                AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
                if length(AutoClusters) >= j
                    Clusters{j} = AutoClusters{j};
                    fprintf('Cluster %d was empty - reverting to Auto\n',j);
                    errordlg(sprintf('Cluster %d was empty Reloaded AutoCluster\n',j),'Cluster Error','modal');                    
                elseif DATA.checkclusters
                    errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',j),'Cluster Error','modal');                    
                end
                    
            end
        end
    else
        Clusters = DataClusters;
    end
    if length(Clusters) >= DATA.probe(1) && isfield(Clusters{DATA.probe(1)},'auto')
    wasauto = Clusters{DATA.probe(1)}.auto;
    else
        wasauto = 1;
    end
        

    if exist(dname,'file') && quickmode == 0
        load(dname);
    elseif isappdata(DATA.toplevel,'ClusterDetails')
        ClusterDetails = getappdata(DATA.toplevel,'ClusterDetails');
    end
    p = ProbeNumber(DATA);

    if DATA.autorefine > 0
        DATA.cluster.manual = 3;
        res.refinemode = DATA.refinemode;
    elseif DATA.cluster.auto ==1
        DATA.cluster.manual = 0;
    else
        DATA.cluster.manual = 1;
    end

% Check DATA.cluster.next before copying to Clusters{p}
    for j = 1:length(DATA.cluster.next)
        DATA.cluster.next{j} = rmfields(DATA.cluster.next{j},'r','clst','rescaled');
        %Keep empty nexts empty
        if isfield(DATA.cluster.next{j},'space')
            DATA.cluster.next{j}.manual = DATA.cluster.manual;
        end
        CheckClusterValues(DATA, DATA.cluster.next{j});
    end      

%if this is the first lap of an autocluster we don't want to keep the stored clsuter.first -
%what happend last time around
    if isempty(DATA.lastcut) && length(Clusters) >= p && isfield(Clusters{p},'first')
        Clusters{p} = rmfield(Clusters{p},'first');
    end
    
    f = fields(DATA.cluster);
    for j = 1:length(f)
        Clusters{p}.(f{j}) = DATA.cluster.(f{j});
    end
    Clusters{p}.probe = p;
    
    xid = {};
    if isfield(DATA,'clst')
        id = find(DATA.clst(DATA.uid) == 2);
        nc = length(unique(DATA.clst(DATA.uid)));
        for j = 3:nc
        xid{j-2} = find(DATA.clst(DATA.uid) == j);
        end
    else
        id = 1:length(DATA.uid);
    end
    [a,b] = Counts(DATA.clst(DATA.uid));
    if length(a) < 2
        a(2) = 0;
    end
    fprintf('Saving (%s)/%d Spikes (P%d) to %s\n',sprintf('%d ', a(2:end)),length(DATA.uid),p,outname);
% ClusterDetails records all event times, and the classification (clst)
%Clusters just has the times of the classified  events = smallest file ffor
%combine
    Clusters{p}.times = DATA.t(DATA.uid(id));
    Clusters{p}.errs = DATA.errs;
    if quickmode == 0
    Clusters{p}.dpsum = sum(abs(DATA.cluster.MeanSpike.vdprime(DATA.chspk,:)),2);
    end
    if ~isempty(DATA.restricttimerange)
        Clusters{p}.restricttimerange = DATA.restricttimerange;
    elseif isfield(Clusters{p},'restricttimerange')
        Clusters{p} = rmfield(Clusters{p},'restricttimerange');
    end
    if ~isempty(DATA.excludetrialids)
        Clusters{p}.excludetrialids = DATA.excludetrialids;
    elseif isfield(Clusters{p},'excludetrialids')
        Clusters{p} = rmfield(Clusters{p},'excludetrialids');
    end
    Clusters{p}.spkfile = SpkFilename(DATA);
    Clusters{p}.exptno = DATA.Expt.exptno;
    ClusterDetails{p}.Evec = DATA.Evec;
%when calling reclassify, ClusterDetails is loaded, so can get evec from
%there. But its cheap (ish- why not have in in Clusters? 
%    Clusters{p}.Evec = DATA.Evec;
    if isfield(DATA.Expt.Header,'ReadMethod')
        DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;
    else
        DATA.cluster.exptreadmethod = 0;
    end

%remove fields that we don't want saved. N.B. might be inherited from an old Clusters file. 
%xtimes is no longer used
%xy, goes into ClusterDetails
%Evec more htan double sinze, so keep this in ClusterDetails. Only need it
%for Reclassify.
    Clusters{p} = rmfields(Clusters{p},'r','xtimes', 'Evec','rescaled');
    Clusters{p}.trigdt = DATA.trigdt;
    Clusters{p}.tsmooth = DATA.triggersmooth;
    Clusters{p}.triggerchan = DATA.triggerchan;
    Clusters{p}.triggertype = DATA.triggertype;
    Clusters{p}.clst  = DATA.clst(DATA.uid);
    Clusters{p}.version = DATA.version;
    Clusters{p}.duration = DATA.duration;
    Clusters{p}.probe = p;
    if isfield(DATA,'DataType')
        Clusters{p}.DataType = DATA.DataType;
    end
    if ~isfield(DATA.cluster,'auto')
        DATA.cluster.auto = 0;
    end
    if DATA.recluster ~= 2
    Clusters{p}.auto = 0;
    else
    Clusters{p}.auto = DATA.cluster.auto;
    Clusters{p}.recluster = 2;
    end
    Clusters{p}.pcmean = mean(DATA.pcs(:,1:4)); %to check for sign reversal
    Clusters{p}.memsz = [memsize(DATA) DATA.fullvsize];
    if isfield(DATA.cluster,'DprimeUsed')
        Clusters{p}.DprimeUsed = DATA.cluster.DprimeUsed;
    end
    if isfield(DATA.cluster,'mumeanUsed')
        Clusters{p}.mumeanUsed = DATA.cluster.mumeanUsed;
    end
    if isfield(DATA.cluster,'TemplateUsed')
        Clusters{p}.TemplateUsed = DATA.cluster.TemplateUsed;
    elseif isfield(DATA,'TemplateUsed')  && ~isfield(Clusters{p},'TemplateUsed')
        Clusters{p}.TemplateUsed = DATA.TemplateUsed;
        if isfield(DATA,'DprimeUsed') % can happne if templates not calculated for this probe
            Clusters{p}.DprimeUsed = DATA.DprimeUsed;
        end
    end
    Clusters{p}.chspk = DATA.chspk;
    Clusters{p}.vsmps = DATA.vsmps;
    if isfield(DATA.cluster,'bestspace')
    if DATA.cluster.bestspace(2) == 3
        Clusters{p}.vspace = DATA.vspace;
    end
    end
    if DATA.checkclusters
        CheckClusters(Clusters,'CheckFitSpace');
    end
    if isfield(Vall,'builddate')
%        Clusters{p}.FullVbuilddate = Vall.builddate;
        FullVData.builddate = Vall.builddate;
    end
    Clusters{p}.isicheck = DATA.isicheck;
    if length(DATA.artifacttimes)
        Clusters{p}.artifacttimes = DATA.artifacttimes;
    end
    Clusters{p}.clusterprog = sprintf('AllVPcs %.2f',DATA.version);
    Clusters{p}.progversion = DATA.version;
    Clusters{p}.missingtrials = DATA.missedtrials;
    if DATA.autorefine > 0 
        Clusters{p}.manual = 3;
    elseif DATA.cluster.auto ==1
        Clusters{p}.manual = 0;
    else
        Clusters{p}.manual = 1;
    end
    CheckClusterValues(DATA, DATA.cluster);

    
    
    Clusters{p} = rmfields(Clusters{p},{'clst' 'r'});
    Clusters{p}.savetime(1) = now;
    if NeedTemplateForCluster(Clusters{p},1) ==2 && isfield(DATA.Template,'othermeans')
        Clusters{p}.MeanSpike.othermeans = DATA.Template.othermeans(2:end);
    end
    if DATA.trigdt == 4 && isfield(DATA,'TriggerTemplate')
        Clusters{p}.TriggerTemplate = DATA.TriggerTemplate;
    end
    if isfield(DATA,'xy')
            ClusterDetails{p}.xy = DATA.xy{1}(DATA.uid,:);
            ClusterDetails{p}.t = DATA.t(DATA.uid);
            ClusterDetails{p}.clst = DATA.clst(DATA.uid);
        if isappdata(DATA.toplevel,'ClusterDetails')
            setappdata(DATA.toplevel,'ClusterDetails',ClusterDetails);
            PlotAllProbes(DATA, 'allxy', 'probes',p,'linewidth',2);
        end
    end
    ClusterDetails{p}.angle = Clusters{p}.angle;
    ClusterDetails{p}.shape = Clusters{p}.shape;
    if isfield(Clusters{p},'crit')
        ClusterDetails{p}.crit = Clusters{p}.crit;
    else
        ClusterDetails{p}.crit = NaN;
    end
    if savexy && isfield(DATA,'xy')
        if DATA.savetrigger
            ClusterDetails{p}.rV = DATA.rV;
        end
        if savexy ==2
            for j = 1:length(Clusters{p}.next)
                if length(DATA.xy) > j && ~isempty(Clusters{p}.next{j}) ...
                        && ~isempty(DATA.xy{j+1})
                    ClusterDetails{p}.next{j}.xy = DATA.xy{j+1}(DATA.uid,:);
                end
            end
            id = DATA.uid;
        else
            minpts = min([ceil(DATA.duration * 50) DATA.nevents]);
            minpts = min([minpts length(DATA.uid)]);
            if length(DATA.clid) > length(DATA.uid)/2
                npts = length(DATA.xy{1});
            elseif length(DATA.clid) < length(DATA.uid)/10;
                npts = length(DATA.uid)/5;
            else
                npts = length(DATA.clid).*2;
            end
            npts = max([ceil(npts) minpts]);
            if mean(DATA.xy{1}(DATA.clid,1)) > mean(DATA.xy{1}(DATA.nid,1))
                [a,id] = sort(DATA.xy{1}(DATA.uid,1),'descend');
            else
                [a,id] = sort(DATA.xy{1}(DATA.uid,1));
            end
            id = sort(id(1:npts));
            if length(id) > DATA.nevents;
                id = id(id < DATA.nevents);
            end
            id = DATA.uid(id);
        end
        if diff(size(ClusterDetails{p}.clst)) > 0
            fprintf('Clst size 2 is %d\n',size(ClusterDetails{p}.clst,2));
            ClusterDetails{p}.clst = ClusterDetails{p}.clst';
        end
        ClusterDetails{p}.ctime = Clusters{p}.ctime;
    else
        id = [];
    end
    %if making an automatic cut, and a manual one was made already, dont 
    %overwrite the spikes file.  To force overwrite, set savespikes to 2.
    if DATA.savespikes ==1 && wasauto == 0 && DATA.cluster.auto == 1
        fprintf('Will not save spikes(auto) - manual cut was made\n');
        id = [];
    end
    
    if quickmode
        DataClusters = Clusters;
    end
    Clusters{p}.savetime(2) = now;
    
    if DATA.nolog == 0 && quickmode == 0
    if exist(logname,'file') 
        load(logname);
        ncl = length(ClusterLog)+1;
    else
        ncl = 1;
    end
    ClusterLog{ncl}.cluster = DATA.currentcluster;
    if DATA.currentcluster > 1
        C = Clusters{p}.next{DATA.currentcluster-1};
    else
        C = Clusters{p};
    end
    ClusterLog{ncl}.shape = C.shape;
    ClusterLog{ncl}.space = C.space;
    ClusterLog{ncl}.mahal = C.mahal;
    if isfield(C,'crit')
        ClusterLog{ncl}.crit = C.crit;
    else
        ClusterLog{ncl}.crit = NaN;
    end
    ClusterLog{ncl}.angle = C.angle;
    if C.shape == 0
        ClusterLog{ncl}.xyr = C.xyr;
    end
    ClusterLog{ncl}.savetime = now;
    ClusterLog{ncl}.user = DATA.user;
    ClusterLog{ncl}.probe = p;
    ClusterLog{ncl}.exptno = DATA.exptno;
    ClusterLog{ncl}.hostname = gethostname;
    ClusterLog{ncl}.recluster = [DATA.recluster DATA.cluster.auto];
    ClusterLog{ncl}.ncut = [C.ncut DATA.nevents];
    save(logname,'ClusterLog');
    if isfield(DATA.cluster,'gmfit')
        nd = size(DATA.cluster.gmfit.mu,2);
    else
        nd = 0;
    end
    if DATA.logfid > 2 && DATA.logfid < 20
        fprintf(DATA.logfid,'E%dP%d Space%s Fit %d dims Saved (%s)/%d spikes at %s\n',...
            DATA.exptno,p,sprintf(' %d',DATA.cluster.space),nd,datestr(now),...
            sprintf('%d ',Counts(DATA.clst)),DATA.nevents);
    end
    end
    
    nerr = 0;
    badc = 0;
    for j = 1:length(Clusters)
        if ~isfield(Clusters{j},'mahal')
            nerr = nerr+1;
            badc = j;
        end
    end
    if nerr > 0 && DATA.checkclusters
        errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',badc),'Cluster Error','modal');
    end

    if isempty(Vall)
        if isfield(DATA,'FullVData')
            FullVData = DATA.FullVData;
        else
            FullVData = [];
        end
    else
        FullVData.blklen = Vall.blklen;
        FullVData.blkstart = Vall.blkstart;
        FullVData.samper = Vall.samper;
        FullVData.missedtrials = DATA.missedtrials;
        if isfield(Vall,'meanV')
            FullVData.submean = 1;
        end
        if isfield(Vall,'builddate')
            FullVData.builddate = Vall.builddate;
        end
        if isfield(Vall,'highpass')
            FullVData.highpass = Vall.highpass;
            Clusters{p}.highpass = Vall.highpass;
        end
    end
    if quickmode == 0
        save(dname,'ClusterDetails','FullVData');
    end
    Clusters{p}.savetime(3) = now;

    if DATA.auto.uselastcluster
    %Keep track of last saved cluster for each probe
        for j = 1:length(Clusters)
            if j > length(DATA.LastClusters) || isempty(DATA.LastClusters{j})
                DATA.LastClusters{j} = Clusters{j};
            end
        end
        DATA.LastClusters{p} = Clusters{p};
    end

    save(outname,'Clusters','FullVData');
    DATA.savespkid = id;
    if DATA.toplevel
        set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));
    end

    
function C = StripClusters(Clusters)
        
for j = 1:length(Clusters)
    C{j} = Clusters{j};
    if isfield(C{j},'xy')
        C{j} = rmfield(C{j},'xy');
    end
end



function name = SpkFilename(DATA)
    [a,b] = fileparts(DATA.name);
    if ~isempty(regexp(b,'Expt[0-9]*Spikes'))
        [c,b] = fileparts(a);
        b = regexprep(b,'Expt[0-9]*Spikes',c);
    end
    if isempty(b)
        [c,b] = fileparts(a);
    end
    b = regexprep(b,'([GM][0-9]*).*','$1');
    mnk = GetMonkeyName(DATA.name);
    if isempty(strfind(b,mnk))
        b = [mnk b];
    end
    
    if isdir(DATA.name)
        a = DATA.name;
    end
    if length(DATA.clustersubdir)
        spkdir = [a '/' DATA.clustersubdir '/Spikes'];
    else
        spkdir = [a '/Spikes'];
    end
    
    xs='';
    if isfield(DATA,'Expt') && rem(DATA.Expt.exptno,1) > 0.001
        xs = 'a';
    end
    name = sprintf('%s/%s.p%dt%d%s.mat',spkdir,b,ProbeNumber(DATA),floor(DATA.Expt.exptno),xs);
        
function probe = ProbeNumber(DATA)
    if isfield(DATA,'probelist')
        if DATA.probe(1) > length(DATA.probelist)
            probe = DATA.probe(1);
        else
            probe = DATA.probelist(DATA.probe(1));
        end
    elseif isfield(DATA,'probe')
        probe = DATA.probe(1);
    else
        probe = 1;
    end
        
        
function SaveSpikes(DATA, id)
    
    if isempty(id)
        return;
    end
    name = SpkFilename(DATA);
    AllV = GetAllV(DATA);
    a = fileparts(name);
    if ~exist(a,'dir')
        mkdir(a);
    end
    p = ProbeNumber(DATA);
    Spikes.values = squeeze(AllV(DATA.probe(1),:, id));
    if size(Spikes.values,2) > 100
        Spikes.values = Spikes.values';
    end
    Spikes.maxv = max(abs(Spikes.values(:)));
    Spikes.maxint = 32000;
    Spikes.values = int16(Spikes.values .* Spikes.maxint./Spikes.maxv);
    Spikes.times = reshape(DATA.t(id),length(id),1);
    xy = DATA.xy{1}(id,1);
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
    if isfield(DATA,'matfile');
        Spikes.Header.matfile = DATA.matfile;
    end
    fprintf('Saving %d Spikes to %s\n',length(id),name);
    save(name,'Spikes');
    if DATA.saveallspikes
        xname = regexprep(name,'.p([0-9])*t','.p$1xt');
        if ~strcmp(name,xname) && size(AllV,1) > 1
        chspk = setdiff(DATA.chspk, DATA.probe(1));
        Spikes.values = squeeze(AllV(chspk,:, id));
        Spikes.TriggerV = DATA.rV(id);
        Spikes.maxv = max(abs(Spikes.values(:)));
        Spikes.maxint = 32000;
        Spikes.values = int16(Spikes.values .* Spikes.maxint./Spikes.maxv);
        fprintf('Saving %d Spikes to %s\n',length(id),name);
        Spikes.chspk = chspk;
        save(xname,'Spikes');
        end
    end
    
function tcut = IsTemplateCut(E)    
   tcut = 0;
   if isfield(E,'space') && E.space(1) ==6 && ismember(E.space(2), [4 7])
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
    [E, Scores, tmpdips, xy, details] = AutoCut(DATA, args{:});
    DataClusters = getappdata(DATA.toplevel,'Clusters');

    if details.pctook > DATA.pcfit.took * 1.1
        fprintf('PC fits %.2f vs %.2f\n',details.pctook, DATA.pcfit.took);
    end

    if ~isempty(Scores) && (~isfield(DATA,'TemplateScores') || E.plottype == 3)
        DATA.TemplateScores = Scores;
        DATA.tmpdips = tmpdips;
    elseif E.newscores
        DATA = get(DATA.toplevel,'UserData');
    end
    if E.angle ~= 0 %auto rotation caused by 1d/2d mismatch
        DATA.xy{1} = xyrotate(xy(:,1),xy(:,2),E.angle);
    else
    DATA.xy{1} = xy;
    end
    DATA.ndxy = xy;
    DATA.cboundary = E;
    if IsTemplateCut(E)
        DATA.plottype = 3;
    else
        DATA.plottype = E.plottype;
    end
    if DATA.watchplots
%    PlotHistogram(DATA,E); %should be called in Classifyspikes
    end
    [cl, DATA.cluster, DATA.xy{1}] = ClassifySpikes(DATA,E);
    DATA.cluster.auto = 1;
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
    if isfield(DATA.Expt.Header,'ReadMethod')
        DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;
    else
        DATA.cluster.exptreadmethod = 0;
    end
    % don't watndot do this unless saving, in which case its done below
% This wipes out fields like savetimes
%    DATA.Clusters{DATA.probe(1)} = DATA.cluster;
    
    
    if nosave == 0
        outname = ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);
        DATA =  SaveClusters(DATA, outname);
        if length(DATA.savespkid) > length(DATA.xy{1})
            fprintf('Save Spike Mismatch %d vs %d\n',length(DATA.savespkid),length(DATA.xy{1}));
        end
        if DATA.savespikes
            SaveSpikes(DATA,DATA.savespkid);
        end

    else
            id = find(cl.clst(DATA.uid) > 1);
% ClusterDetails records all event times, and the classification (clst)
%Clusters just has the times of the classified  events = smallest file ffor
%combine
        DataClusters{DATA.probe(1)}.times = DATA.t(DATA.uid(id));
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
        spaces(1,:) = DATA.cluster.space;
        spaces(2,:) = [3 1 8];
        spaces(3,:) = [3 2 10];
        spaces(4,:) = [3 1 8];
    elseif DATA.usestdtemplates
        spaces(1,:) = [3 1 2];
    else
        spaces(1,:) = [3 1 8];
        spaces(2,:) = [3 2 10];
        spaces(3,:) = [3 1 8];
    end
    for j = 1:size(spaces,1);
        [objs{j+1}, distance(j+1), a] = GMfit(DATA.TemplateScores(:,spaces(j,2:3)),nc,ntr,'idlist',idlist);
    end
    if DATA.usestdtemplates
        nd = find(DATA.tmplspace(3,:) > 0);
        [objs{1}, distance(1), a] = GMfit(DATA.TemplateScores(:,DATA.tmplspace(3,nd)),nc,ntr,'idlist',idlist);
        btype = 7;
    else
        [objs{1}, distance(1), a] = GMfit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),nc,ntr,'idlist',idlist);
        btype = 4;
    end
%        objs{5} = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),3,'Options',statset('MaxIter',1000));
        details.Converged(1) = objs{1}.Converged;
        if mean(a.fail) == 1
            distance(1) = 0;
            C.sign = 0;
        elseif ~isfield(DATA,'clid') || isempty(DATA.clid) || recalc
            xy = ProjectND(DATA, btype, objs{1});
            %        [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
            [C.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
            details.usegmcluster = 0;
            C.sign = details.sign;
            if details.sign >= 0
                DATA.clid = find(xy(:,1) > C.crit(1));
                DATA.nid = find(xy(:,1) <= C.crit(1));
            else
                DATA.clid = find(xy(:,1) < C.crit(1));
                DATA.nid = find(xy(:,1) >= C.crit(1));
            end
            clst(DATA.clid) = 2;
            clst(DATA.nid) = 1;
            SetFigure(DATA.tag.tmplscore, DATA);
            PlotXY(xy,clst);
            set(gca,'UserData',[NaN 6 btype]);
        else
            C.sign = 0;
        end
        obj = objs{1};
        details.Gs = objs;
        [a,b] = max(distance);
        if distance(b) > distance(1)
            [C.crit, details] = GMDip(DATA.TemplateScores(:,spaces(b-1,2:3)),DATA.energy(1,:),'label',DATA.idstr);
            details.bestcl = cluster(objs{b},DATA.TemplateScores(:,spaces(b-1,2:3)));
            details.space = spaces(b-1,:);
            xy = DATA.TemplateScores(:,spaces(b-1,2:3));
        elseif DATA.usestdtemplates
            details.bestcl = cluster(objs{1},DATA.TemplateScores(:,DATA.tmplspace(3,nd)));
            details.space = [6 7];
        else
            details.bestcl = cluster(objs{1},DATA.TemplateScores(:,DATA.tmplspace(1,:)));
            details.space = [6 4];
        end
        details.cluster = C;
        if max(distance) > details.gmdprime * 1.5 %1D cut is very poor
            details.usegmcluster = 1;
        else
            details.usegmcluster = 0;
        end
        if DATA.logfid > 2
            fprintf(DATA.logfid,'P%d, T%d at %s\r\n',DATA.probe(1),DATA.cluster.probe,datestr(now));
        end
        
function DATA = TemplateGMFits(DATA)        
    p= DATA.probe(1);
    uid = DATA.uid;
    %when doing this from scratch (remaking template from other cut), use
    %several starts to make sure the PC cut is best. 
    %Now uses GMfit that set a sensible start point.
    
    tmspace = DATA.tmplspace(2,:);
    tmspace = tmspace(tmspace>0);
    pcs = DATA.TemplateScores(uid,tmspace);
    ncells = 2:DATA.ncelltotry;
    for j = 1:length(ncells)
        [P, distance(j), a] = GMfit(pcs,ncells(j),1);
        details.pctook = a.took;
        objs{j} = P;
        [d, dd] = gmdprime(P);
        cid = cluster(P, pcs);% mahal distance is unsigned, not so useful
        if ncells(j) == DATA.ncelltotry(2)
            DATA.gmcid = cid;
            DATA.usegmcid = 1;
            distances = dd.d;
        end
        details.Converged(j) = P.Converged;
        ll(j) = P.NlogL;
    end
    GetFigure(DATA.tag.dips);
    mysubplot(1,4,1:3);
    imagesc(distances);
    caxis([0 5]);
    mysubplot(1,4,4);
    plot(ll,1+[1:length(ll)],'o-');
    set(DATA.toplevel,'UserData',DATA);

    
function [distance, obj, xy, details] = BestSpace(DATA, varargin)

    newtemplate = 0;
    nloops = 0; %to test multiple fits for consistency
    pconly = 0;
    nr = 1;
    ntr = 1;
    nc = 2;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'newtemplate',7);
            newtemplate = 1;
            nr=5;
        elseif strncmpi(varargin{j},'template',5)
            ntr=5; %redoing a template fit. Make sure its right.
        elseif strncmpi(varargin{j},'pconly',5)
            pconly = 1; %redoing a template fit. Make sure its right.
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
    details.pctook = a.took;
    objs{1} = P;
    details.Converged(1) = P.Converged;
    ll(1) = P.NlogL;
    if pconly
        best = 1;
    obj = objs{best};
   [xy, details.cid] = ProjectND(DATA, best, objs{best});
   e(1) = mean(DATA.energy(1,details.cid == 1));
   e(2) = mean(DATA.energy(1,details.cid == 2));
   if e(1) > e(2) && nc < 3
       details.cid = 3 - details.cid;
   end
   return;
    end
%The Var/Energy plot is particularly sensitive to electdoe movement
%artifacts, so exclude spikes after trial end when making the cut
   AllV = GetAllV(DATA);
    iid = setdiff(uid, DATA.postevents);
    xy(:,1) = DATA.energy(1,iid);
    xy(:,2) = DATA.spkvar(p,iid)./DATA.energy(1,iid);
    
    [E, distance(2)] = GMfit(xy,nc,1);
    details.Converged(2) = E.Converged;
    objs{2} = E;
    distance(2) = gmdistance(E);
    ll(2) = E.NlogL;
    
    try
    [V, distance(3), a] = GMfit(squeeze(AllV(p,DATA.vspace,uid))',nc,1);
    catch
        fprintf('ERROR!!! E%dP%d uid %d-%d of %d, clusters %s\n',DATA.exptno,p,min(uid),max(uid),size(AllV,3),nc);
    end
    objs{3} = V;
    details.Converged(3) = V.Converged;

    ll(3) = V.NlogL;
    quick = 1;
    if newtemplate || ~isfield(DATA,'TemplateScores')
        [d, best] = max(distance);
        [xy, details.cid] = ProjectND(DATA, best, objs{best});
%        [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
        %could use gmdistribution/cluster to assign to groups here.
        if quick & length(unique(details.cid)) > 1
            DATA.clid = find(details.cid  == 2);
            DATA.nid = find(details.cid  == 1);
            DATA.clst = details.cid;
        else
            [cluster.crit, a] = GMDip(xy(uid,:),DATA.energy(1,uid),'label',DATA.idstr);
            cluster.sign = a.sign;
            if a.sign >= 0
                DATA.clid = find(xy(:,1) > cluster.crit(1));
                DATA.nid = find(xy(:,1) <= cluster.crit(1));
            else
                DATA.clid = find(xy(:,1) < cluster.crit(1));
                DATA.nid = find(xy(:,1) >= cluster.crit(1));
            end
            DATA.clst(DATA.clid) = 2;
            DATA.clst(DATA.nid) = 1;
        end
        DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
        if DATA.currentcluster == 1
            DATA.cluster.MeanSpike = DATA.MeanSpike;
            DATA.cluster.spts = DATA.spts;
        else
            DATA.cluster.next{DATA.currentcluster-1}.MeanSpike = DATA.MeanSpike;
        end
        TemplatePlot(DATA,'nodip','usemean');
        DATA = get(DATA.toplevel,'UserData');
%        imean = DATA.MeanSpike;  seems unused
        ntr = 1;
    end
    if isfield(DATA,'TemplateScores')

        [objs{4}, distance(4), a] = GMfit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,ntr);
%        objs{5} = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),3,'Options',statset('MaxIter',1000));
        distance(4) = gmdistance(objs{4});
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
    if ismember(best,[2 3]) 
        safeid = [1 4];
        [d, best] = max(distance(safeid));
        best = safeid(best);
        err = sprintf('Using Space %d (%.1f, %.1f) not Space 3 (%.1f. %.1f)\r\n',best,distance(best),ll(best),distance(3),ll(3));
        fprintf('%s',err);
        if DATA.logfid > 2
            fprintf(DATA.logfid,'%s',err);
        end
    end
    details.bestd = d;
    details.besti = best;
    if best == 4
        DATA.cluster.space = [6 4];
        DATA.cluster.shape = 2;
        if ~isfield(DATA,'clid') || isempty(DATA.clid)
        xy = ProjectND(DATA, best, objs{best});
%        [cluster.crit, details] =
%        FindDip(xy(:,1),DATA.energy(1,:));
        [cluster.crit, details] = GMDip(xy(uid,:),DATA.energy(1,uid),'label',DATA.idstr);
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
        distance(4) = gmdistance(objs{4});
        DATA = get(DATA.toplevel,'UserData');
    end
    
    
    
    details.TemplateUsed = DATA.TemplateUsed;
    
    obj = objs{best};
   [xy, details.cid] = ProjectND(DATA, best, objs{best});
   e(1) = mean(DATA.energy(1,details.cid == 1));
   e(2) = mean(DATA.energy(1,details.cid == 2));
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
    ts = now;
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
       all.fail(j) = 0;
        catch
            fprintf('GM Fit fail\n');
            all.Gs{j} = S;
            all.Gs{j}.Converged = -1;
            all.Gs{j}.NlogL = NaN;
            all.fail(j) = 1;
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
            all.fail(j) = 1;
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
            all.fail(j) = 1;
        end
    end
    if all.ds(j) > max(all.ds(1:j-1)) * 1.1
        [a,b] = max(all.ds(1:j-1));
           fprintf('Best manual start %.2f vs %.2f (%d of %d)\n',all.ds(j),a,b,j-1);
    end
    [D,b] = max(all.ds);
    G = all.Gs{b};
    all.took = mytoc(ts);


function [G, xy] = IterateTemplateFit(DATA, G)
j = 1;
xc = [0 0];
DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
T{1} = DATA.MeanSpike.ms;
nc = size(G.mu,1);
    while j < 10 && xc(1,2) < 0.5
    xy = ProjectND(DATA, 4, G);
%    [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
    [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
    cluster.sign = details.sign;
    if details.sign >= 0
        DATA.clid = find(xy(:,1) > cluster.crit(1));
        DATA.nid = find(xy(:,1) <= cluster.crit(1));
    else
        DATA.clid = find(xy(:,1) < cluster.crit(1));
        DATA.nid = find(xy(:,1) >= cluster.crit(1));
    end
    DATA.cluster.neednewtemplate = 1;
    DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
    TemplatePlot(DATA,'nodip','noplot');
    DATA = get(DATA.toplevel,'UserData');
    T{j+1} = DATA.MeanSpike.ms;
    [G, distance(j),a] = GMfit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),nc,1);
%    G = gmdistribution.fit(DATA.TemplateScores(:,DATA.tmplspace(1,:)),2,'Options',statset('MaxIter',1000));
%    D = mahal(G,G.mu);
%    distance(j) = sqrt(2./(1./D(1,2)+1./D(2,1)));
    xc = corrcoef(T{j+1}(:),T{j}(:));
    xcs(j) = xc(1,2);
    j = j+1;
    end
    DATA.xy{1} = xy;
    DATA.ndxy = xy;
    set(DATA.toplevel,'UserData',DATA);
        
function [xy, cid, details] = ProjectND(DATA, best, obj)
    pnorm = 0;
    details.err = 0;
    if size(obj.mu,1) > 2
        v = range(obj.mu);
    else
        v = diff(obj.mu);
    end
    if best ==1 && size(v,2) ~= 4
        details.err = 1;
        xy(:,1) = DATA.pcs(:,1);
        xy(:,2) = DATA.pcs(:,2);
        cid = [];
        DATA = AddErr(DATA,'ERROR! Dimension mismatch\n');
        return;
    end
        if best == 1
            sd = std(DATA.pcs(:,[1 2 3 4]));
            v = v(1:4)./sd;
        if pnorm
            S(:,1) = DATA.pcs(:,1) ./ sd(1);
            S(:,2) = DATA.pcs(:,2) ./ sd(2);
            S(:,3) = DATA.pcs(:,3) ./ sd(3);
            S(:,4) = DATA.pcs(:,4) ./ sd(4);
        else
            S = DATA.pcs(:,1:4);
        end

    elseif best == 2

    elseif best == 3
        if isfield(DATA,'AllV')
        S = squeeze(DATA.AllV(DATA.probe(1),DATA.vspace,:))';
        else
        AllV = GetAllV(DATA);
        S = squeeze(AllV(DATA.probe(1),DATA.vspace,:))';
        end
    elseif ismember(best, [4 7])
        nd = find(DATA.tmplspace(best-3,:) > 0);
        if pnorm
        sd = std(DATA.TemplateScores(:,DATA.tmplspace(best-3,nd)));
        v = v./sd;
        for j = 1:nd
            S(:,j) = DATA.TemplateScores(:,DATA.tmplspace(1,j)) ./ sd(j);
        end
        else
            S = DATA.TemplateScores(:,DATA.tmplspace(best-3,nd));
        end
%        xy(:,1) = v([1 5]) * S(:,[1 5])'; %test
        end
        
        if best == 2
            xy(:,1) = DATA.energy(1,:);
            xy(:,2) = DATA.spkvar(DATA.probe(1),:)./DATA.energy(1,:);
            S = xy;
        elseif size(v,2) == size(S,2)
            xy(:,1) = v * S';
        else
            fprintf('ERROR! Dimension mismatch\n');
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
    if size(ov,2) == size(S,2)
    xy(:,2) = ov * S';
%apply an arbitray sign convetion, so that this isnt' just set byt
%the order of teh two means in the fit
fixsign = 0;
    if mean(xy(:,1)) < 0 && fixsign
        xy(:,1) = -xy(:,1);
        xy(:,2) = -xy(:,2);
    end
    else
        xy = S(:,1:2);
    end
    if size(obj.mu,2) ~= size(S,2) 
            DATA = AddErr(DATA,'ERRROR! Dimensions of space don''t match fit\n');
            showerrs = 0;
            if showerrs
            errordlg(sprintf('Fit %d dimensions, data %d',size(obj.mu,2),size(S,2)),'Fit Dimension Error','modal');
            end
            cid = ones(1,size(S,1));
            details.err = 1;
    elseif ~isobject(obj)
            DATA = AddErr(DATA,'ERRROR! Final GM fit Failed\n');
            cid = ones(1,size(S,1));
            details.err = 2;
    else
        cid = cluster(obj, S);% mahal distance is unsigned, not so useful
    end
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
AllV = GetAllV(DATA);
ctime = now;
if usegm
    DATA.cluster.spts = DATA.spts;
    DATA.cluster.probe = DATA.probe;
    if ~isfield(DATA.cluster,'next')
        DATA.cluster.next = {};
    end
    for j = 1:length(DATA.ncomponents)
    [bd,obj{j},xy,c] = BestSpace(DATA,'newtemplate','ncells',DATA.ncomponents(j));
    end
    E.newscores = 1;
    DATA.ndxy = xy;
    [a,b] = max(bd);
    E.bestspace = [c.bestd c.besti];
    E.bestd = c.besti;
    E.bestll = c.ll;
    E.bestcl = c.cid;
    E.bestfitres = c.Converged;
    if isfield(c,'TemplateUsed')
        E.TemplateUsed = c.TemplateUsed;
    end
    E.auto = 1;
    if c.besti == 4
        [gm,d]= max(bd(1:3));
        s = sprintf('GM %.2f for space %d (%.2f for %d)dt%d csd%d',a,c.besti,gm,d,DATA.dvdt,DATA.csd);
    else
        s = sprintf('GM %.2f for space %s dt%d csd%d',a,DATA.SpaceTypes{b},DATA.dvdt,DATA.csd);
    end
    fprintf('%s\n',s);
     PrintMsg(DATA.logfid,'P%d (S%d) %s %d events at %s\r\n',DATA.probe(1),DATA.savespikes,s,DATA.nevents,datestr(now));
%    [dip, details] = FindDip(xy(:,1),DATA.energy(1,:),'gmix');
    [dip, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
   [E.gmfit2d, d] = GMfit(xy,2,1,'idlist',c.cid);
   E.mahal(4) = d;
   E.angle = 0;
   if d > 2  && d > details.gmdprime * 1.1 %% 2D significantly better - check rotation
       [a, gm,  dipres] = BestAngleGM(xy, E.gmfit2d, details);
       if gm > 2 && gm > details.gmdprime * 1.5
%           xy = dipres.xy; %rotated values
           details.dipres = dipres.gmfit;
           details.gmdprime = gm;
           E.angle = a;
           [dip, details] = GMDip(dipres.xy,DATA.energy(1,:),'label',DATA.idstr);
       end
   end
   E.mahal(4) = details.gmdprime;
   details.pctook = c.pctook;
    
    E.space = [6 c.besti];
    E.bestd = bd;
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.gmdip = dip;
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    %in case rotation was applied by BestAngle above
    x = xyrotate([crit crit],[min(xy(:,2)) max(xy(:,2))],-E.angle);
    E.pos = x([1 3 2 4]);
%    E.pos = [crit min(xy(:,2)) crit max(xy(:,2))];
    E.plottype = DATA.plottype;
    E.gmfit = obj{end};
    if length(obj) > 1
        E.gmfits = obj(1:end-1);
    end
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
        [as(j+n),bs(j+n),c] = BestAngle(AllV(p(j,1),p(j,2),:),AllV(p(j,3),p(j,4),:),2);
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
        x = AllV(p(j,1),p(j,2),:);
        y = AllV(p(j,3),p(j,4),:);
        E.plottype = 2;
    else
        p= DATA.pcplots;
        cluster.space = [1 DATA.pcplots(j,:)];
        x = DATA.pcs(:,p(j,1));
        y = DATA.pcs(:,p(j,2));
        E.plottype = 1;
    end
    xy = xyrotate(x,y,cluster.angle);
%    [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
    [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
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
        [Scores, details.TemplateUsed, details.DprimeUsed] = TemplatePlot(DATA);
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
%                [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
                [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
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
%            [cluster.crit, details] = FindDip(xy(:,1),DATA.energy(1,:));
              [cluster.crit, details] = GMDip(xy,DATA.energy(1,:),'label',DATA.idstr);
            cluster.sign = details.sign;
        end
    end
    cluster.auto = 1;
    E.autotook = mytoc(ctime);
    if length(cluster.space) > 1 && cluster.space(1) == 6 && cluster.space(2) == 4
        E.plottype = 3;
    end
    details.newDATA = newDATA;
    E = BoundaryFromCluster(E,cluster, DATA.currentcluster);
    

    
function DATA = AddErr(DATA,varargin)


s = sprintf(varargin{:});
fprintf('%s\n',s);
if isfield(DATA,'errs')
    DATA.nerr = length(DATA.errs)+1;
else
    DATA.nerr = 1;
end
DATA.errs{DATA.nerr} = s;

if isfield(DATA,'logfid') && DATA.logfid > 0
    try
        fprintf(DATA.logfid,'%s\r\n',s);
    catch
        fprintf('fid %d no longer valid\n',DATA.logfid);
    end
end


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

    
   SetFigure(DATA.tag.vare, DATA);
   subplot(2,2,1);
    x = DATA.energy(1,:);
    y = DATA.spkvar(DATA.probe(1),:)./DATA.energy(1,:);
    Cs(1) = CutAndPlot(x,y,DATA.energy(1,:));
    drawnow;

    subplot(2,2,2);
    x = DATA.energy(1,:).^2;
    y = DATA.spkvar(DATA.probe(1),:).^2 ./DATA.energy(1,:).^2;
    Cs(2) = CutAndPlot(x,y,DATA.energy(1,:));
    drawnow;
    C = Cs(1);
    C.sign = 0;

function [DATA, E] = OptimizeBoundary(DATA);
    C = OptimizeClusterBoundary(DATA);
%    DataClusters{DATA.probe(1)}.angle = C.angle;
%    DataClusters{DATA.probe(1)}.crit = C.crit;
    if DATA.currentcluster > 1
        c = DATA.currentcluster-1;
        DATA.cluster.next{c}.angle = C.angle;
        DATA.cluster.next{c}.crit = C.crit;
        DATA.cluster.next{c}.sign = C.sign;
        if isfield(C,'xyr')
            DATA.cluster.next{c}.xyr = C.xyr;
        end
        DATA.cluster.next{c}.manual = 3;
        DATA.cluster.next{c}.y = C.y;
        DATA.cluster.next{c}.fitdprime = C.fitdprime;
    else
        DATA.cluster.angle = C.angle;
        DATA.cluster.crit = C.crit;
        DATA.cluster.sign = C.sign;
        if isfield(C,'xyr')
            DATA.cluster.xyr = C.xyr;
        end
        DATA.cluster.manual = 3;
        DATA.cluster.y = C.y;
        if isfield(C,'fitdprime')
        DATA.cluster.fitdprime = C.fitdprime;
        end
    end
    E = BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);
    [cl, cluster, DATA.xy{1}] = ClassifySpikes(DATA,DATA.cluster,DATA.quickcutmode);
    if isfield(cluster,'sign')
        E.sign = cluster.sign;
    end
    if isfield(cluster,'fitdpparams')
        E.fitdpparams = cluster.fitdpparams;
        E.fitdprime = cluster.fitdprime;
    end
    DATA.cluster = cluster;
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;


    
function [C, details] = OptimizeClusterBoundary(DATA)
    
    details = [];
    
    if DATA.currentcluster > 1
        preC = DATA.cluster.next{DATA.currentcluster-1};
    else
        preC = DATA.cluster;
    end
    
    
    if preC.shape == 0
        [C,details]  = OptimizeEllipse(DATA);
        C.crit = 1;
        C.sign = 0;
        return;
    end
    if preC.shape == 1 || preC.shape == 2
        [C, details]  = OptimizeLine(DATA);
        C.sign = CheckSign(C,details.r,DATA.energy);
        return;
    end
    usegm = 1;
    
    space = preC.space;
    if space(1) == 5
        C = OptimizeVarE(DATA);
        return;
    elseif space(1) == 6  %% cut in > 2 dimensions.
        if DATA.currentcluster == 1 %should always be tru
            x = DATA.ndxy(:,1);
            y = DATA.ndxy(:,2);
        else
            x = DATA.xy{DATA.currentcluster}(:,1);
            y = DATA.xy{DATA.currentcluster}(:,2);
        end
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
        AllV = GetAllV(DATA);
        x = squeeze(AllV(space(2),xi,:));
        y = squeeze(AllV(space(4),yi,:));
    elseif space(1) == 3
        x = DATA.TemplateScores(:,xi);
        y = DATA.TemplateScores(:,yi);
    end
    C.y = minmax(y);
    
    if usegm
        xy = cat(2,x,y);
        [a,b,c] = BestAngleGM(xy, [],[]);
        C.angle = a;
        C.crit = c.crit;
        C.sign = CheckSign(C,c.xy(:,1),DATA.energy);
        return;
    end
    a = -pi/2:pi/36:pi/2; % use this range because this is what atan returns;
    
    for j = 1:length(a)
        xy = xyrotate(x,y,a(j));
        dip(j) = HartigansDipTest(sort(xy(:,1)));
        [aa,bb] = oldFindDip(xy(:,1),DATA.energy(1,:));
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
    

    [crit,b] = oldFindDip(xy(:,1),DATA.energy(1,:),'plot',DATA.tag.dips,'gmix');
    C.crit = crit(1);
    C.crit = mean(crit([5 6])); %based on GM fit
    C.sign = b.sign;
    SetFigure(DATA.tag.dips, DATA)
    hold off;
    plot(a,dip./max(dip));
    hold on;
    plot(a,mydip./max(mydip),'r');
    plot(a,coeff./max(coeff),'g');
    
    
    
    
    
function FullV = SetFullVNames(DATA, FullV)
    
    [a,b] = fileparts(FullV.loadname);
    eid = sscanf(b,'Expt%d');
    if isfield(DATA,'exptname')
        [c,d, e] = fileparts(DATA.exptname);
        FullV.matfile = [a '/' d e];
    elseif isfield(FullV, 'matfile')
        [c,d, e] = fileparts(FullV.matfile);
        FullV.matfile = [a '/' d e];
    end
    

function PlotXY(xy, clst)    

    ptsz = 1;
    if length(xy) < 1000
        ptsz = 6;
    end
  hold off;
  plot(xy(:,1),xy(:,2),'.','markersize',ptsz);
  hold on;
  id = find(clst == 2);
  plot(xy(id,1),xy(id,2),'r.','markersize',ptsz);
  
 function c = BimodalCoeff(x, e)
    e = 1.3;
    c = (1+skewness(x).^2)./((kurtosis(x).^e)+3);
    
    
    
function ChangeCell(a,b,p)
    DATA = GetDataFromFig(a);
    if isfield(DATA,'CellDetails');
        eid = find(DATA.CellDetails.exptids == DATA.exptno);
    else
        eid = DATA.exptno;
    end
    ecells = squeeze(DATA.CellList(eid,:,:));
    [x,y] = find(ecells == p);
    ChangeProbe(a,b,x);
    
function ChangeProbe(a,b,p)
    DATA = GetDataFromFig(a);
    if strcmp(p,'next')
        if DATA.probelist(DATA.probe(1)) < DATA.allnprobes
            p = DATA.probelist(DATA.probe(1))+1;
        else
            p = 0;
        end
    elseif strcmp(p,'prev')
        if DATA.probelist(DATA.probe(1)) > 2
            p = DATA.probelist(DATA.probe(1))-1;
        else
            p = 0;
        end
    elseif sum(strcmp(p,{'save' 'quicksave'}))
        if strcmp(p,'save')
            PCCluster(a,b,25); %save spikes and cluster
        else
            PCCluster(a,b,12); %save spikes and cluster
        end
        if DATA.probe(1) < DATA.allnprobes
            p = DATA.probelist(DATA.probe(1))+1;
        else
            p = 0;
        end
    end
    set(DATA.toplevel,'name',sprintf('Loading probe %d',p));
    drawnow;
    if ~isfield(DATA,'spkrate')
        DATA.spkrate = 100;
    end
    if p > 0
        if strncmp(DATA.DataType,'GridData',8)
            if DATA.loadfromspikes
                newname = DATA.name;
            else
                newname = regexprep(DATA.fullvname,'\.p[0-9]*FullV',sprintf('.p%dFullV',p));
            end
            if DATA.setspkrate > 0
            AllVPcs(newname,'tchan',p,DATA.probeswitchmode,'spkrate',DATA.setspkrate);
            else
            AllVPcs(newname,'tchan',p,DATA.probeswitchmode);
            end
        else
            if DATA.setspkrate > 0
                AllVPcs(DATA.toplevel,'tchannew',p,DATA.probeswitchmode,'spkrate',DATA.setspkrate);
            else
                AllVPcs(DATA.toplevel,'tchannew',p,DATA.probeswitchmode);
            end
        end
    end
    set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));
    if DATA.quickcutmode.plotspikes
        DATA = get(DATA.toplevel,'UserData');
    end
    figure(DATA.toplevel);

    
function OptionMenu(a,b, fcn)

    onoff = {'off' 'on'};
[DATA, F] = GetDataFromFig(a);
qf = fields(DATA.quickcutmode);
autofields = fields(DATA.auto);
fullvf = fields(DATA.fullvswitchmode);
if strmatch(fcn,'plotspikes')
    DATA.quickcutmode.(fcn) = ~DATA.quickcutmode.(fcn);
    set(a,'checked',onoff{1+DATA.quickcutmode.(fcn)});
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,'comment')
    str  = sprintf('E%dP%d',DATA.exptno,ProbeNumber(DATA));
    GetString(DATA.tag.comments,DATA.toplevel,{@AddComment, []},'label',str);
    return;
elseif strmatch(fcn,'keepspikes')
    set(F,'Tag','OldSpikes');    
    SetFigure(DATA.tag.spikes, DATA,'front');
elseif strmatch(fcn,'keepmeanspikes')
    id = unique(DATA.clst);
    for j = 2:length(id)
        ms{j} = PlotMeanSpike(DATA,'recalc','cluster',j-1);
    end
    figure;
    voff = DATA.voffset - DATA.voffset(DATA.probe(1));
    for j = 2:length(id)
        for k = DATA.chspk
        plot(ms{j}.ms(k,:)+voff(k),'color',DATA.colors{j},'linewidth',2);
        hold on;
        end
    end
    dv = mean(diff(voff(DATA.chspk)))/10;
    for k = DATA.chspk
        text(1,voff(k)+dv,sprintf('E%dP%d',DATA.exptno,k),'fontsize',DATA.gui.fontsize(1));
    end
elseif strmatch(fcn,'spoolspikes')
    SetFigure(DATA.tag.spikes, DATA,'front');
    SpoolSpikes(DATA); %ellipse
    figure(DATA.toplevel);
elseif strmatch(fcn,'autoopt')
    flag = get(a,'tag');
    DATA.auto.(flag) = ~DATA.auto.(flag);
    set(a,'checked',onoff{1+DATA.auto.(flag)});
    DATA.checkclusters = DATA.auto.checkcluster;
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,'checkclusters')
    DATA.checkclusters = ~DATA.checkclusters;
    set(a,'checked',onoff{1+DATA.checkclusters});
    set(DATA.toplevel,'UserData',DATA);

elseif strmatch(fcn,qf)
    if ~strcmp(fcn,'quickest')
        DATA.quickcutmode.quickest = 0;
    end
    DATA.quickcutmode.(fcn) = ~DATA.quickcutmode.(fcn);
    if DATA.quickcutmode.quickest == 1
        for j = 1:length(qf)
            DATA.quickcutmode.(qf{j}) = 0;
        end
        DATA.quickcutmode.quickest = 1;
    end
        
    c = get(get(a,'parent'),'children');
    tags = get(c,'Tag');
    for j = 2:length(c)
        k = strmatch(get(c(j),'Tag'),qf);
        set(c(j),'checked',onoff{1+DATA.quickcutmode.(qf{k})});
    end
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,fullvf)
    DATA.fullvswitchmode.(fcn) = ~DATA.fullvswitchmode.(fcn);
    SetMenuChecks(get(a,'parent'),DATA.fullvswitchmode);
    set(DATA.toplevel,'UserData',DATA);
elseif regexp(fcn,'[1-9]probefullv')
    DATA.plotspk.nfullvprobes = sscanf(fcn,'%d');
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,{'plotsummary'})
    DATA = LoadTrigTimes(DATA,1:DATA.nprobes,'savexy');
    SpikeDraw(DATA,[],  'allquickspks');
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,{'nextfullv'})
    AllVPcs(DATA.toplevel, 'newexpt', DATA.exptno+1, DATA.probeswitchmode);
elseif strmatch(fcn,{'prevfullv'})
    AllVPcs(DATA.toplevel, 'newexpt', DATA.exptno-1, DATA.probeswitchmode);
elseif strmatch(fcn,{'newfullv'})
    s = inputdlg({'New Expt #' 'Probe #'},'New Expt to load',1,{num2str(DATA.exptno+1) num2str(DATA.probe(1))});
    if isempty(s)
        return;
    end
    eid = str2num(s{1});
    pid = str2num(s{2});
    args{1} = 'tchan';
    args{2} = pid;
    args{3} = DATA.probeswitchmode;
    if DATA.fullvswitchmode.summary
        args = {args{:} 'summary'};
    end   
    AllVPcs(DATA.toplevel, 'newexpt', eid, args{:});
elseif strmatch(fcn,{'profiling'})
    DATA.profiling = ~DATA.profiling;
    set(a,'checked',onoff{1+DATA.profiling});
    set(DATA.toplevel,'UserData',DATA);
elseif strmatch(fcn,{'retrigger'})
    args = RetriggerDialog(a,b,'popup');
    if ~isempty(args)
        AllVPcs(DATA.toplevel,'tchan',DATA.probe, args{:});
    end
elseif strmatch(fcn,{'showfullv' 'includeprepost'})
    DATA.plotspk.(fcn) = ~DATA.plotspk.(fcn);
    set(a,'checked',onoff{DATA.plotspk.(fcn)+1});
    set(DATA.toplevel,'UserData',DATA);
end
figure(DATA.toplevel);



    function pos = PlaceUi(a, b ,str)
        space = 0.01;
        pos = get(b,'Position');
        x = get(a,'position');
        if strcmp(str,'left');
            pos(1) = x(1)+x(3);
            pos(2) = x(2);
        elseif strcmp(str,'down');
            pos(1) = x(1);
            pos(2) = x(2)-x(4);
        elseif strcmp(str,'up');
            pos(1) = x(1);
            pos(2) = x(2)+x(4);
        end
        set(b,'position', pos);
        
 function args = RetriggerDialog(a,b, fcn)
        DATA = GetDataFromFig(a);
        persistent trigargs;
        persistent trigstate;
        
        p = ProbeNumber(DATA);
        if strcmp(fcn,'popup')
        F = dialog('Name','Retrigger','resize','on');
        a =uicontrol(F,'style','pushbutton','string','ReTrigger','callback',{@RetriggerDialog, 'go'});
        b =uicontrol(F,'style','pushbutton','string','Cancel','callback',{@RetriggerDialog, 'cancel'});
        PlaceUi(a,b,'up');
        b =uicontrol(F,'style','pushbutton','string','SetSpkrate','callback',{@RetriggerDialog, 'setspkrate'});
        PlaceUi(a,b,'left');
        e = uicontrol(F,'style','pop','string','V|dVdt|accel|energy|Template|Templatedt','tag','trigtype');
        PlaceUi(b,e,'left');
        c = uicontrol(F,'style','text','string','Threshold','tag','labelthrehsold');
        PlaceUi(e,c,'left');
        s = num2str(DATA.cluster.Trigger);
        trigstate.Trigger = str2num(s);
        d = uicontrol(F,'style','edit','string',s,'tag','threshold');
        PlaceUi(c,d,'left');

        aa =uicontrol(F,'style','text','string','Smooth','tag','labelthrehsold');
        PlaceUi(c,aa,'up');
        s = num2str(DATA.triggersmooth);
        trigstate.smooth = str2num(s);
        bb = uicontrol(F,'style','edit','string',s,'tag','smooth');
        PlaceUi(aa,bb,'left');

        e = uicontrol(F,'style','text','string','Spkrate','tag','labelspkrate');
        PlaceUi(d,e,'left');
        if isfield(DATA.cluster,'eventrate')
            trigstate.spkrate = round(DATA.cluster.eventrate .* DATA.samplerate );
        else
            trigstate.spkrate = round(DATA.cluster.nspks./DATA.duration);
        end
        if trigstate.Trigger < 0
            trigstate.spkrate = -trigstate.spkrate;
        end
        f = uicontrol(F,'style','edit','string',num2str(trigstate.spkrate),'tag','spkrate');
        PlaceUi(e,f,'left');
        FitWindow(F);
        set(F,'UserData',DATA.toplevel);
        uiwait(F);
        args = trigargs;
        elseif strcmp(fcn,'cancel')
            args = {};
            close(gcf);
        elseif strcmp(fcn,'setspkrate')
            F = get(a,'parent');
            DATA.setspkrate = str2num(get(findobj(F,'tag','spkrate'),'string'));
            set(DATA.toplevel,'UserData',DATA);
        elseif strcmp(fcn,'go')
            F = get(a,'parent');
            trigtypes = {'threshold','dtthr','dtthr2','dtthr3','triggertemplate','triggertemplatedt'};
            trigtype = get(findobj(F,'tag','trigtype'),'value');
            a = str2num(get(findobj(F,'tag','threshold'),'string'));
            b = str2num(get(findobj(F,'tag','spkrate'),'string'));
            if ismember(trigtype,[5 6])
                if isfield(DATA.cluster,'MeanSpike')
                    Ms = DATA.cluster.MeanSpike.ms;
                else
                    DataClusters = getappdata(DATA.toplevel,'Clusters');
                    Ms = DataClusters{p}.MeanSpike.ms;
                end
                npts = size(Ms,2);
                spts = 16-npts:15;
                trigargs  = {trigtypes{trigtype} Ms 'spts' spts};
            else
            trigargs  =  {trigtypes{trigtype }};
            end
% by default use spkrate to set trigger - hard to guess what you need for
% other spaces.
            if  a ~= trigstate.Trigger
                trigargs = {trigargs{:} a};
            else
%                if trigtype == 1
%                    trigargs = {};
%                end
                if b > 0
                    trigargs = {trigargs{:} 'th+' 'spkrate', b};
                else
                    trigargs = {trigargs{:} 'th-' 'spkrate', -b};
                end
            end
            a = str2num(get(findobj(F,'tag','smooth'),'string'));
            trigargs = {trigargs{:} 'smooth' a};
            trigargs = {trigargs{:} 'reapply'};
            close(gcf);
        end
        
function FitWindow(F)
c = get(F,'children');
set(c,'units','pixels');
b = get(c,'position');
 x = cat(1,b{:});
 w = max(x(:,1)+x(:,3));
 h = max(x(:,2)+x(:,4));
 b  = get(F,'position');
 set(F,'position',[b(1) b(2) w+min(x(1,:)) h+min(x(1,:))]);
set(c,'units','normalized');
 
 
function SetMenuChecks(hm, S)
%Set checked on/off for menu items in hm whose tags
%match fields in structre S
qf = fields(S);
onoff = {'off' 'on'};

c = get(hm,'children');
tags = get(c,'Tag');
for j = 1:length(c)
    k = strmatch(get(c(j),'Tag'),qf);
    if length(k) == 1
    set(c(j),'checked',onoff{1+S.(qf{k})});
    end
end


function ProbeMenu(a,b, fcn)

    onoff = {'off' 'on'};
[DATA, F] = GetDataFromFig(a);
if strmatch(fcn,{'usecluster' 'reclassify' 'autocut' 'simple' 'reapply'})
    DATA.probeswitchmode = fcn;
    SetMenuCheck(a,'exclusive');
    set(DATA.toplevel,'UserData',DATA);
elseif strcmp(fcn,'select')
    s = inputdlg('New Probe #','Change Probe',1,{num2str(ProbeNumber(DATA)+1)});
    p = str2num(s{1});
    args = {};
    if DATA.fullvswitchmode.summary
        args = {args{:} 'summary'};
    end
    ChangeProbe(DATA,[],p);
elseif ismember(fcn,[1 2 3])
    [C, DATA.Evec, DATA.pcs, DATA.dipvals, DATA.chspk] = CalcPCs(DATA,AllV,fcn-1);
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
elseif fcn == 4
    DATA.dvdt = ~DATA.dvdt;
    set(a,'Checked',onoff{DATA.plotdvdt+1});
end
    
    
function res = PlotClusters(a,b,fcn)
[DATA, F] = GetDataFromFig(a);
SetFigure('Clusters', DATA);
DataClusters = getappdata(DATA.toplevel,'Clusters');
res = [];
if fcn == 1
    PlotSpikeTimes(DataClusters);
elseif fcn == 2
    res = PlotSpikeTimes(DataClusters,'xcorr');
    for j = 1:length(DataClusters)
        DataClusters{j}.synci = res.synci(j,:);
    end
    DATA.xcorrs = res;
    set(DATA.toplevel,'UserData',DATA);
elseif fcn == 3
    PlotSpikeTimes(DataClusters,'probequality');
end
setappdata(DATA.toplevel,'Clusters',DataClusters);

function SetCellFromLine(a,b, cluster, cell)
    DATA = GetDataFromFig(a);
    DATA = SetCellEntry(DATA, DATA.cluster, DATA.exptno,DATA.probe(1),cluster, cell);
    if max(cell) > length(DATA.comparecell)
        DATA.comparecell(1+length(DATA.comparecell):cell) = 0;
    end
    set(DATA.toplevel,'UserData',DATA);
    PlotCellList(DATA,'showfig');
    
    
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
    f = SetFigure(DATA.tag.celllist, DATA);
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
        set(DATA.toplevel,'UserData',DATA);
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
eid = find(DATA.CellDetails.exptids == DATA.exptno);
PlotCellIm(DATA.CellList, DATA.CellDetails, DATA.nclusters);
PlotCellIm('mark',DATA.tagged,'y');
PlotCellIm('marklist',[ [DATA.Comments.ex]; [DATA.Comments.p];],'g');
text(ProbeNumber(DATA),0,sprintf('P%d',ProbeNumber(DATA)),'color','r','fontsize',DATA.gui.fontsize(1));
text(-1,eid,sprintf('E%.0f',DATA.exptno),'color','r','fontsize',DATA.gui.fontsize(1));
return;
offset = 0;
nc = max(DATA.nclusters); %max # clusters in any one expt, for each probe
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
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)),'fontsize',DATA.gui.fontsize(1));
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
    h = text(y(b),x(b)-1,sprintf('%d',cells(j)),'fontsize',DATA.gui.fontsize(1));
    set(h,'color',colors{j});
end
iscellim = sum(DATA.CellList,3) > 0;

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
    oldcell = DATA.CellList(e,p,c);
    id = find(DATA.CellList(e,p,:) == cellid);
    DATA.CellList(e,p,id) = 0; %remove this cell from other clusters this expt
    DATA.CellList(e,p,c) = cellid;
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
    if CellIsEmpty(DATA.CellDetails,'MeanSpike', cellid)
        DATA.CellDetails.MeanSpike{cellid}.ms = C.MeanSpike.ms(DATA.chspk,:);
        DATA.CellDetails.MeanSpike{cellid}.chspk = DATA.chspk;
    end
    if savelist
        SaveCellList(DATA);
    end

function SaveCellList(DATA)
    cellfile = [DATA.name '/CellList.mat'];
    CellList = DATA.CellList;
    CellDetails = DATA.CellDetails;
    CellChanges = DATA.CellChanges;
    save(cellfile, 'CellList','CellDetails','CellChanges');

function DATA = AddAxisContextMenu(DATA, ax)
    
    if ~isfield(DATA,'CellList')
        return;
    end
    cmenu = uicontextmenu;
    
    cellid = unique(DATA.CellList);
    cellid = cellid(cellid>0);

    for j = 1:length(cellid)
        h = uimenu(cmenu,'label',sprintf('Cell %d',cellid(j)),'Callback',{@SetCellCompare, cellid(j)});
        if DATA.comparecell(cellid(j))
            set(h,'checked','on');
        end
    end
    set(ax,'UIContextMenu',cmenu);
    
    function cmenu = AddLineContextMenu(DATA, h)


        if ~isfield(DATA,'CellList')
            cmenu = [];
            return;
        end
        
        h = h(h>0);
        p = ProbeNumber(DATA);
        if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'exptids')
            eid = find(DATA.CellDetails.exptids == DATA.exptno);
        else
            eid = DATA.exptno;
        end
        
        ecells = squeeze(DATA.CellList(eid,:,:));
        pcells = squeeze(DATA.CellList(:,p,:));
        pid = find(sum(pcells,2) > 0);  %probes with defined cells
        if isempty(pid) || isempty(eid)
            nearestex = 0;
            nearestcell = 0;
        else
            [a,b] = min(abs(eid-pid));
            nearestex = pid(b); %nearest expt with this probe defined as cell
            nearestcell = DATA.CellList(nearestex,p,:);
        end
    for j = 2:length(h)
        cmenu = uicontextmenu;
        if j > size(DATA.CellList,3)+1 || DATA.exptno > size(DATA.CellList,1)
            k = 0;
        else
            k = DATA.CellList(DATA.exptno,DATA.probe(1),j-1);
        end
        if k > 0
            uimenu(cmenu,'label',sprintf('Spike%d Cell%d',j-1,k),'foregroundcolor',get(h(j),'color'));
            c = uimenu(cmenu,'label','clear','Callback',{@DeleteCellFromLine, j-1,  k});
        else
            uimenu(cmenu,'label',sprintf('Spike%d',j-1),'foregroundcolor',get(h(j),'color'));
        end
        ccolor = get(h(j),'color');
        if size(ecells,2) >= j-1 && size(ecells,1) >= p
        ccell = ecells(p,j-1);
        else
            ccell = 0;
        end
        for k = 1:20
            c = uimenu(cmenu,'label',sprintf('Cell %d',k),'Callback',{@SetCellFromLine, j-1,  k});
            if k == 1
                set(c,'separator','on');
            end
            if k == ccell
                set(c,'foregroundcolor',ccolor)
            elseif sum(ecells(:) == k) %this cell already defined
                set(c,'foregroundcolor',[0 0.8 0])
            elseif ismember(k,nearestcell)
                cl = find(nearestcell == k);
                set(c,'label',sprintf('Cell %d (E%d)',k,nearestex),'foregroundcolor',DATA.colors{cl+1});
            end
                
        end
        set(h(j),'uicontextmenu',cmenu);
    end

function M = CalcDistanceMatrices(DATA, nc, varargin)
%overfit with GMs, and calculate distance Matrix
    p= DATA.probe(1);
    AllV = GetAllV(DATA);
    uid = DATA.uid;   
    ntr = 1;
    M.errs = {};
    M.errstates = {};
    nerr = 0;
    try
    objs{1} = gmdistribution.fit(DATA.pcs(uid,DATA.pcspace),nc,'Options',statset('MaxIter',1000));
    catch ME
        objs{1} = [];
        M = AddErr(M,'GMFit ERROR!!!A E%dP%d uid %d-%d of %d, clusters %d (%s)\n',DATA.exptno,p,min(uid),max(uid),size(AllV,3),nc,ME.message);
        nerr = length(M.errs);
        M.errstates{nerr} = ME;
    end
    try
        objs{2} = gmdistribution.fit(squeeze(AllV(p,DATA.vspace,uid))',nc,'Options',statset('MaxIter',1000));
    catch ME
        objs{2} = [];
        M = AddErr(M,'GMFit ERROR!!!B E%dP%d uid %d-%d of %d, clusters %d (%s)\n',DATA.exptno,p,min(uid),max(uid),size(AllV,3),nc,ME.message);
        nerr = length(M.errs);
        fprintf('AllV %s, Vpcs %s p%d',sprintf('%d ',size(AllV)),sprintf('%d ',DATA.vspace),p);
        M.errstates{nerr} = ME;
    end
    DATA = SetTemplateData(DATA,1, 'force');
    try
    objs{3} = gmdistribution.fit(DATA.TemplateScores(uid,DATA.tmplspace(1,:)),nc,'Options',statset('MaxIter',1000));
    catch ME
        M = AddErr(M,'Error Fitting  Gms C (%s)\n',ME.message);
        nerr = length(M.errs);
        M.errstates{nerr} = ME;
    end
    nf = 0;
    goodfit = zeros(1,3); %in fits all crash, length(objs) is 0
    for j = 1:length(objs)
        if ~isempty(objs{j})
            nf = nf +1;
            M.converged(nf) = objs{j}.Converged;
            [a,b] = gmdprime(objs{j});
            M.D(:,:,nf) = b.d;
            goodfit(j) = 1;
        end
    end
    best = find(goodfit);
    if ~isempty(best)
        best = best(1);
        [a,b(1)] = gmdprime(objs{best});
        f = fieldnames(objs{best});
        if isprop(objs{1},'mu')
            M.cid = cluster(objs{best}, DATA.pcs(uid,DATA.pcspace));
            M.gmfit = objs{best}; %save the fit
        else
            M.cid = zeros(size(DATA.uid));
        end
    end
    [a, M.isacell] = isacell(DATA, DATA.exptno, p);
    if isfield(DATA.cluster,'mahal')
     M.mahal(1) = DATA.cluster.mahal(4);
    end
    for j = 1:length(DATA.cluster.next)
        if isfield(DATA.cluster.next{j},'mahal')
            M.mahal(j+1) = DATA.cluster.next{j}.mahal(4);
        end
        if isfield(DATA.cluster.next{j},'fitdprime')
            M.fitdprime(j+1) = DATA.cluster.next{j}.fitdprime(1);
        end
    end
    if isfield(DATA.cluster,'fitdprime')
     M.fitdprime(1) = DATA.cluster.fitdprime(1);
    end
    

function cid = AssignCluster(DATA, G)
    cid = cluster(G, DATA.pcs(DATA.uid,DATA.pcspace));

function SetPlot(a,b, fcn)
    [DATA, F] = GetDataFromFig(a);
    if strcmp(fcn,'testscore')
        DATA.plottype = 12;
    elseif strcmp(fcn,'allspace')
        DATA.plottype = 11;
    end
    ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    
function PCCluster(a,b, fcn)

[DATA, F] = GetDataFromFig(a);
DataClusters = getappdata(DATA.toplevel,'Clusters');
onoff = {'off' 'on'};
plotpcs = 0;
finish = 0;

if strcmp(fcn,'all6cells')
    b = CalcDistanceMatrices(DATA,6);
    uid = DATA.uid;
    for j = 1:11
        objs{j} = gmdistribution.fit(DATA.pcs(uid,DATA.pcspace),6,'Options',statset('MaxIter',1000));
        [a,d] = gmdprime(objs{j});
        [D, scores] = squishDistanceMatrix(d.d);
        bestd(j) = scores(2);
        fprintf('Fit %d D = %.2f\n',j,scores(2));
    end
    [a,c] = sort(bestd);
    besti = c(6);
%    [a, besti] = min(bestd);
    DATA.gmcid = AssignCluster(DATA, objs{besti});
    [a,D] = gmdprime(objs{besti});
    b.D(:,:,1) = D.d;

    [nspk, a] = Counts(DATA.gmcid);
    DATA.usegmcid = 1;
    GetFigure('DistanceMatrix');
    [E,V] = eig(squeeze(b.D(:,:,1)));
    [c,d] = sort(E(:,1));
    D = MatrixPermute(b.D(:,:,1),d);
    [D, scores, dlist] = squishDistanceMatrix(squeeze(D(:,:,1)));
    imagesc(D);
    title(sprintf('Score %.2f %.2f',scores(1),scores(2)));
    caxis([0 5]);
    for j = 1:6
        text(j,1,num2str(nspk(dlist(j))));
    end
    ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
end

if ~ismember(fcn,[4 7])  %if intereactive, look at plots when change cuts
    DATA.watchplots = 1;
    DATA.watcharg = {'front'};
    DATA.interactive = 1;
end
classargs = {};

if DATA.profiling
    profile on;
end
if fcn == 1
    if isempty(a)
        return;
    elseif isstruct(a)
        E = b;
         if E.boundarytype == 2
            E.LineA = BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);
        end
    else
        DATA.clustermode = 0;
        DATA = SetEllipseDrawing(DATA, 0); %called from gui
        return;
    end
    DATA.usegmcid = 0;
    classargs = {classargs{:} DATA.quickcutmode};
%    plotpcs = 1;
elseif strncmp(fcn,'Clear2',5)
    c = sscanf(fcn(6:end),'%d');
    if isfield(DATA.cluster,'next')
        DATA.cluster.next{c-1} = [];
    end
    DATA.clst(DATA.clst ==c+1) = 1;
    DATA.currentcluster = 1;
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif strncmp(fcn,'Clear2',5)
    c = sscanf(fcn(6:end),'%d');
    if isfield(DATA.cluster,'next')
        DATA.cluster.next{c} = [];
    end
    DATA.clst(DATA.clst ==c+1) = 1;
    DATA.currentcluster = 1;
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif strncmp(fcn,'plotcelllist',5)

elseif strncmp(fcn,'iteratefit',9)
    IterateFit(DATA, 1);
    DATA = get(DATA.toplevel,'UserData');
    finish = 1;
elseif strncmp(fcn,'mark',4)
    DATA.cluster.marked = find(strcmp(fcn,{'marknone' 'markgood' 'markbadprobe' 'markmugood'}));
    finish = 1;
elseif strncmp(fcn,'savedef',7)
    outname = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    if DATA.auto.saveref
        SaveClusters(DATA, outname,'quick');
    end
    outname = regexprep(outname,'Expt[0-9]*ClusterTimes','RefClusters');
    SaveClusters(DATA, outname,'quick');
    if DATA.auto.advanceprobe
        ChangeProbe(DATA,[],'next');
        return;
    end
    finish = 1;
elseif strncmp(fcn,'plotother',7)
    DATA.plottype = 9;
    ReplotPCs(DATA,[]);
    return;
elseif strncmp(fcn,'dvdtpc',6)
    DATA.dvdt = ~DATA.dvdt;
    set(a,'Checked',onoff{1+DATA.dvdt});
    DATA = SetPCs(DATA, 1);
    ReplotPCs(DATA,[]);
elseif strncmp(fcn,'plotica',7)
    DATA.plottype = 10;
    if 1 || ~isfield(DATA,'icas')
        DATA = CalcICA(DATA,5);
    end
    set(DATA.toplevel,'UserData',DATA);
    ReplotPCs(DATA,[]);
    return;
elseif strcmp(fcn,'quantifyquick')
    QuantifyQuickClusters(DATA,[1:DATA.nprobes],'savespikes');
    return;
elseif strncmp(fcn,'Ellipse',7)
    DATA.usegmcid = 0;
    if isstruct(a)
        E = b;
    else
        c = sscanf(fcn(8:end),'%d');
        if c == 0 %use line for cluster 1
            c = 1;
            clustertype = 1;
        else
            clustertype = 0;
        end
        if c ~= DATA.currentcluster && DATA.auto.replotcluster
            replot = 1;
        else
            replot = 0;
        end
        DATA = SetEllipseDrawing(DATA, clustertype,'cluster',c); %called from gui
%        DATA.clustericon = SetClusterIcon(DATA);
        if replot
            TemplatePlot(DATA,'recalc');
            DATA = get(DATA.toplevel,'UserData');
            ReplotPCs(DATA,[], 'autospace');
        end
        set(DATA.toplevel,'UserData',DATA);
        return;
    end
    manualcut = 1;
    plotpcs = 1;
elseif strncmp(fcn,'flipsign',7)
    [cl, cluster] = ClassifySpikes(DATA,DATA.cluster,'sign',DATA.cluster.sign*-1,'quick', 1);
    DATA = ReplotPCs(DATA,[]);
    DATA.cluster = rmfield(cluster,'r');
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    DATA.cluster.ctime = now;
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 2
    it = SetFigure(DATA.tag.vare,DATA, DATA.watcharg{:});
    DATA = SetEllipseDrawing(DATA, 0,'cluster',1,'figure',DATA.tag.vare); %called from gui

    return;
    
elseif fcn == 3 %replot
    PlotMeanSpike(DATA);
    DATA = ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 4 %toggle density plot
    DATA.clplot = ~DATA.clplot;
    if isstruct(a)
        it = findobj(DATA.toplevel,'tag','Density');
        set(it,'Checked',onoff{DATA.clplot+1});
    else
        set(a,'Checked',onoff{DATA.clplot+1});
    end
    ReplotPCs(DATA,[]);
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif ismember(fcn, [5 32])
    if fcn == 5
        TemplatePlot(DATA,'recalc');
    else
        TemplatePlot(DATA,'stdtemplate','nodip');
    end
    DATA = get(DATA.toplevel,'UserData');
    DATA.plottype = 3;
%    DATA = ReplotPCs(DATA,[]); %done in templateplot
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 6
    it = SetFigure(DATA.tag.tmplscore,DATA, DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    E.pcplot = get(gca,'UserData');
    PlotHistogram(DATA,E);
    return;
elseif fcn == 7
    it = SetFigure(DATA.tag.vare,DATA, DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    E.space = 5;
    E.pcplot = [];
    DATA.cboundary = PlotHistogram(DATA,E);
elseif fcn == 8
    DATA = SetEllipseDrawing(DATA, 1,'cluster',1); %called from gui

    return;
elseif fcn == 11
    DATA = SetEllipseDrawing(DATA, 1,'cluster',1,'boundarytype',2); %called from gui
            
    return;
elseif ismember(fcn, [8 11])
    DATA.watchplots = 1;
    it = SetFigure(DATA.tag.top,DATA,DATA.watcharg{:});
    E = AddEllipse(it,'wait','color','r','line');
    if isempty(E)
        return;
    end
    E.cluster = DATA.currentcluster;
    E.pcplot = get(gca,'UserData');
    classargs = {classargs{:} DATA.quickcutmode};
    plotpcs = 1;
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
    outname = ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);
    if fcn == 25 %force save Spikes
        DATA.savespikes = 2;
    end
    stime = now;
    if DATA.profiling
        profile on;
    end
    args = {};
    if fcn == 12 
        args = {args{:} 'quick'};
    end
    [DATA, id] = SaveClusters(DATA,outname,args{:});
    if fcn == 25 %save Spikes
        SaveSpikes(DATA,id);
    end
        fprintf('Save took %.2f\n',mytoc(stime));
    if DATA.profiling
        profile viewer;
    end
    set(DATA.toplevel,'UserData', DATA);
    if DATA.auto.advanceprobe
        ChangeProbe(DATA,[],'next');
    elseif DATA.auto.advanceexpt
        AllVPcs(DATA.toplevel, 'newexpt', DATA.exptno+1, DATA.probeswitchmode);
    end
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
%        DATA = ReplotPCs(DATA,BoundaryFromCluster([],DATA.cluster, DATA.currentcluster));
        DATA = ReplotPCs(DATA,[]);
        SetFigure(DATA.tag.tmplscore, DATA);
        subplot(1,1,1);
        hold off;
        plot(DATA.TemplateScores(:,2)./DATA.spkvar(DATA.probe(1),:)',DATA.TemplateScores(:,2),'.');
%        plot(DATA.spkvar(DATA.probe(1),:),DATA.TemplateScores(:,2),'.');
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
    DataClusters{DATA.probe(1)} = [];
    setappdata(DATA.toplevel,'Clusters',DataClusters);
    return;
elseif fcn == 18  %optimnize boundary line in current plot
    set(DATA.toplevel,'Name','Optimizing');
    drawnow;
    [DATA, E] = OptimizeBoundary(DATA);
    DATA.cluster.ctime = now;
    PlotHistogram(DATA,E);
    set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 19  %make Automatic cut
    set(DATA.toplevel,'Name','Thinking');
    drawnow;
    DATA.watchplots = 1; 
    [E, scores, dips, xy, details] = AutoCut(DATA,'refine');
    if details.newDATA
        DATA = get(DATA.toplevel,'UserData');
    end
    if E.angle ~= 0
        DATA.xy{1} = xyrotate(xy(:,1),xy(:,2),E.angle);
    end
    DATA.ndxy = xy;
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
    DATA.plottype = 8;
    set(DATA.toplevel,'UserData',DATA);
    DATA = ReplotPCs(DATA,[]);
    return;
elseif fcn == 21
    if isfield(DATA,'Clusters') && length(DataClusters) >= DATA.probe(1) && ~isempty(DataClusters{DATA.probe(1)})
        DATA.cluster = DataClusters{DATA.probe(1)};
        DATA = CheckTemplates(DATA,DATA.cluster);
        if DATA.cluster.space(1) == 6
                   [DATA.ndxy, DATA.gmcid] = ProjectND(DATA, DATA.cluster.space(2), DATA.cluster.gmfit);
                   DATA.cluster.bestcl = DATA.gmcid;
        end
        DATA.watchplots = 1;  %must be in gui
        [cl, cluster] = ClassifySpikes(DATA,DATA.cluster);
        DATA.cluster.MeanSpike = cl.MeanSpike;
        DATA.clid = cl.id;
        DATA.clst = cl.clst;
        DATA.cluster.ctime = cluster.ctime;
        PlotHistogram(DATA,[]);
    end
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 22 %evaluate Gaussian means fit in the cluster space
    a = FitGaussMeans(DATA.xy{1},2,'verbose');
    DataClusters{DATA.probe(1)}.mahal = [a.mahal a.dprime];
    SetFigure(DATA.tag.hist, DATA);
    subplot(2,1,2);
    hold on;
    ezcontour(@(x,y)pdf(a.obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
    set(DATA.toplevel,'UserData',DATA);
    setappdata(DATA.toplevel,'Clusters',DataClusters);
    return;
elseif sum(strcmp(fcn, {'2cells' '3cells' '4cells'})) | ...
        ismember(fcn,[23 26 27 31]) %Check spacee out with guassian mixture model.
    usegm = 0;
    ntest = 1; % set high to test that minima are reliably found
 
    nc = find(strcmp(fcn, {'2cells' '3cells' '4cells'}));
    if usegm
        [E, scores, dips, DATA.ndxy] = AutoCut(DATA, 'usegm');
        if E.angle ~= 0
            DATA.xy{1} = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),E.angle);
        end

    else
        if nc
            [d, obj, xy, details] = BestSpace(DATA, 'ncells',nc+2);
        elseif fcn == 26
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
        elseif fcn == 31
            [d, obj, xy, details] = BestSpace(DATA, 'pconly');
            E.bestspace = [1 d];
            E.bestd = [d 0 0 0];
        else
            [d, obj, xy, details] = BestSpace(DATA);
        end
    if ~isfield(DATA,'TemplateScores')
        DATA = get(DATA.toplevel,'UserData');
    end
    DATA.xy{1} = xy;
    DATA.ndxy = xy;
    E.bestcl = details.cid;
    E.gmfit = obj;
    [a,b] = max(d);
    if b == 4
        [c,d]= max(d(1:3));
        fprintf('GM %.2f for space %d (%.2f for %d)dt%d csd%d\n',a,b,c,d,DATA.dvdt,DATA.csd);
    else
        fprintf('GM %.2f for space %d dt%d csd%d\n',a,b,DATA.dvdt,DATA.csd);
    end

%    [dip, details] = FindDip(DATA.xy{1}(:,1),DATA.energy(1,:),'gmix');
    [dip, details] = GMDip(DATA.xy{1},DATA.energy(1,:),'gmix','label',DATA.idstr);
    E.space = [6 b];
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    E.pos = [crit min(DATA.xy{1}(:,2)) crit max(DATA.xy{1}(:,2))];
        
    end
    DATA.gmcid = E.bestcl;
    DATA.cboundary = E;
    DATA.currentcluster = 1;
    [cl, cluster] = ClassifySpikes(DATA,E);
    if fcn == 31
%        [P, d, a] = GMfit(DATA.pcs(DATA.uid,DATA.pcspace),2,1,'idlist',cl.clst);
    end
    E.mahal = cluster.mahal;
    PlotHistogram(DATA,E,'plotgm');
    DATA.clid = cl.id;
    DATA.nid = cl.nid;
    DATA.clst = cl.clst;
    if fcn == 27
        DATA.clst = E.bestcl+1;
    end
    cluster.MeanSpike = cl.MeanSpike;
    DATA.cluster = cluster;
    DATA.cluster.gmfit = obj;
    DATA.MeanSpike = cl.MeanSpike;
    if DATA.gmtypes(b)
        DATA.plottype = DATA.gmtypes(b);
        DATA = ReplotPCs(DATA,E);
    else
        SetFigure(DATA.tag.vare, DATA,DATA.watcharg{:});
        subplot(1,1,1);
        PlotVarE(DATA);
        hold on;
        ezcontour(@(x,y)pdf(obj,[x y]),get(gca,'xlim'),get(gca,'ylim'));
    end
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 24
    SetFigure('Cluster', DATA,DATA.watcharg{:});
    xy = xyrotate(DATA.xy{1}(:,1),DATA.xy{1}(:,2),-DATA.cluster.angle);
    hold off;
    x = mean(xy);
    if DATA.plotspk.muscale < 1
    a = prctile(abs(DATA.xy{1}(:,1) - DATA.cluster.crit),2);
    id = find(abs(DATA.xy{1}(:,1) - DATA.cluster.crit) < a);
    r = rand(size(id)).*a;
    xid = find(abs(DATA.xy{1}(id,1) - DATA.cluster.crit) > r);
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
elseif strcmp(fcn, 'PCmultiple')
    for j = 1:3
        [d, obj, xy, details] = BestSpace(DATA, 'pconly','ncells',j+1);
        ll(j) = obj.NlogL;
        E.bestspace = [1 d];
        E.bestd = [d 0 0 0];
        [d, dd] = gmdprime(obj);
        distance(j,1) = min(dd.d(dd.d > 0));
    end
    DATA.xy{1} = xy;
    [d, dd] = gmdprime(obj);
    E.bestcl = details.cid;
    DATA.clst = details.cid;
    E.gmfit = obj;
    [dip, details] = GMDip(DATA.xy{1},DATA.energy(1,:),'gmix','label',DATA.idstr);
    E.space = [6 1];
    E.pcplot = [2 8 10 11 12];
    crit = dip(1);
    E.xyr = [dip(1) 0];
    E.shape = 2;
    E.sign = details.sign;
    E.pos = [crit min(DATA.xy{1}(:,2)) crit max(DATA.xy{1}(:,2))];
    GetFigure(DATA.tag.dips);
    hold off;
    mysubplot(1,4,1:3);
    imagesc(dd.d);
    mysubplot(1,4,4);
    plot(ll,1+[1:length(ll)],'o-');
    DATA.plottype = 1;
    ReplotPCs(DATA,E);
    set(DATA.toplevel,'UserData',DATA);
    return;

elseif strcmp(fcn, 'NCellTemplate')
    DATA = TemplateGMFits(DATA);
    ReplotPCs(DATA,[]);
    return;
elseif strcmp(fcn, 'NDTemplate') || strcmp(fcn, 'NDStdTemplate') %Classify in template space
    if strcmp(fcn,'NDStdTemplate') && (DATA.usestdtemplates == 0 || ~isfield(DATA,'TemplateScores'))
        TemplatePlot(DATA,'stdtemplate','nodip');
        DATA = get(DATA.toplevel,'UserData');
    end
        [a,b, DATA.xy{1}, details] = TemplateSpace(DATA,'template','recalc');
        if details.space(1) == 6
            DATA.ndxy = DATA.xy{1};
        end
        DATA.cluster.gmfit = b;
        DATA.cluster.gmdprime = details.gmdprime;
        DATA.cluster.sign = details.cluster.sign;
        DATA.cluster.crit = details.cluster.crit;
        if DATA.usestdtemplates
            fprintf('Template GM %.2f 4D %.2f 2d, %.2f1D\n',a(1),a(2),details.gmdprime)
        else
            fprintf('Template GM %.2f 6D %.2f 2d, %.2f1D\n',a(1),a(2),details.gmdprime)
        end
        DATA.cluster.space = details.space;
        if details.space(1) == 6
            DATA.plottype = details.space(2);
            DATA.cluster.shape = 2;
        else
            DATA.plottype = details.space(1);
            DATA.cluster.shape = 1;
        end
        DATA.cluster.angle = 0;
        DATA.cluster.templatesrc = DATA.cluster.probe;
        DATA.cluster.bestcl = details.bestcl;
        DATA.cluster.usegmcluster = details.usegmcluster;
elseif fcn == 33
    if DATA.hidecluster == 0
        DATA = NextPCs(DATA);
        set(a,'label','Use All Spikes')
    else
        DATA = UseAllEvents(DATA);
        set(a,'label','PCs Without Cell')
        ReplotPCs(DATA,[]);
    end
    set(DATA.toplevel,'UserData',DATA);
    return;
elseif fcn == 34
    ts = now;
    DATA.recluster = 0;
    ClassifyAll(DATA,0);
    mytoc(ts);
    return;
else
    for j = 1:size(DATA.pcplots,1);
        yl = prctile(DATA.pcs(:,DATA.pcplots(j,2)),[100-fcn fcn]);
        xl = prctile(DATA.pcs(:,DATA.pcplots(j,1)),[100-fcn fcn]);
        SetFigure(DATA.tag.top,DATA);
        subplot(2,4,j);
        set(gca,'xlim',xl,'ylim',yl);
    end
    return;
end

if finish
    set(DATA.toplevel,'UserData',DATA);
    figure(DATA.toplevel);
    return;
end

DATA.recluster = 0;
DATA.cluster.auto = 0;
ts = now;
if ~exist('E','var')
    E = BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);
else
    f = fields(DATA.cluster);
    for j = 1:length(f)
        if strmatch(f{j},{'sign'}) 
        elseif ~isfield(E,f{j})
            E.(f{j}) = DATA.cluster.(f{j});
        end
    end
end
if isfield(E,'bestcl')
    DATA.gmcid = E.bestcl;
end
set(DATA.toplevel,'Name','Classifying Spikes...');
drawnow;
[cl, cluster, DATA.xy{DATA.currentcluster}] = ClassifySpikes(DATA,E,classargs{:});
if DATA.checkclusters
    CheckClusters(DataClusters, 'CheckNexts','Clasified');
end
if isfield(cluster,'r')
%    cluster = rmfield(cluster,'r');
end
DATA.cluster = cluster;
DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(E);
if isfield(cl,'MeanSpike')
    if DATA.currentcluster == 1
    DATA.cluster.MeanSpike = cl.MeanSpike;
    else
    DATA.cluster.next{DATA.currentcluster-1}.MeanSpike = cl.MeanSpike;
    end
    DATA.MeanSpike = cl.MeanSpike;
end

     if DATA.plot.expt
         DATA.Expt = PlotExptCounts(DATA);
     end
if fcn == 11 
    if isfield(DATA.cluster,'distance')
    SetFigure(DATA.tag.tmplscore, DATA);
    plot(DATA.cluster.distance,r,'.','markersize',1);
    end
end
DATA.clid = cl.id;
DATA.nid = cl.nid;
    DATA.clst = cl.clst;

if ismember(fcn, [1 2 6 7 8 11 30])%cluster cut called from GUI, so replot
    if isfield(cluster,'gmfit2dman') && cluster.gmfit2dman.Converged > -1
        E.gmfit2dman = cluster.gmfit2dman;
        E.mahal = cluster.mahal;
    end
    if DATA.watchplots == 0 || plotpcs%otherwise done in ClassifySpikes
    DATA = ReplotPCs(DATA,DATA.cluster);
    end
    f = {'r' 'gmfit1d' 'fitdprime' 'fitdpparams'};
    C = GetSubCluster(DATA.cluster, DATA.currentcluster);
    for j = 1:length(f)
        if isfield(C,f{j})
            E.(f{j}) = C.(f{j});
        end
    end
    PlotHistogram(DATA, E,classargs{:});  %? add a 'quick' option here too. But most often want the GM disp number
end
mytoc(ts);
if DATA.plot.expt
    SetFigure(DATA.tag.Expt, DATA);
    DATA.Expt = PlotExptCounts(DATA);
end
if DATA.profiling
    profile viewer;
end
set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));
set(DATA.toplevel,'UserData',DATA);
figure(DATA.toplevel);

function C = GetSubCluster(Cluster, c)
%    return struct with details for a given cluster
if c ==1
    C = Cluster;
elseif c <= length(Cluster.next)+1
    C = Cluster.next{c-1};
else
    C = [];
end

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

function [Expt, matfile] = LoadExpt(DATA, ei)
    aargs = {};
    if DATA.usealltrials
        aargs = {aargs{:} 'usealltrials'};
    end
        ei = floor(ei); %plain .mat files are not split just because FullV files are
        if regexp(DATA.name,'Expt[0-9]*Spikes')
            monk = GetMonkeyName(DATA.name);
            [a,b] = fileparts(DATA.name);
            [c,d] = fileparts(a);
            smrname = regexprep(DATA.name,'Expt([0-9]*)Spikes',[monk d]);
        elseif regexp(DATA.name,'lem/M[0-9]*')
            smrname = regexprep(DATA.name,'lem/M([0-9]*).*','$0/lemM$1');
            smrname = regexprep(smrname,'online/lemM([0-9]*).*','$0/lemM$1');
        elseif regexp(DATA.name,'dae/M[0-9]*')
            smrname = regexprep(DATA.name,'dae/M([0-9]*).*','$0/daeM$1');
            smrname = regexprep(smrname,'online/lemM([0-9]*).*','$0/daeM$1');
        else %for file moved to unusual directories
            monk = GetMonkeyName(DATA.name);
            smrname = regexprep(DATA.name,'/M([0-9]*).*',['$0/' monk 'M$1']);
        end
        exfile = [smrname '.' num2str(ei) 'idx.mat'];
        matfile = [smrname '.' num2str(ei) '.mat'];
        if exist(exfile,'file') || exist(matfile,'file') %either will do
            ts = now;
            if exist(exfile,'file')
                fprintf('Loading %s (%s)\n',exfile,IDStr(DATA));
                [Trials, Expts] = APlaySpkFile(exfile,'noerrs','nospikes', aargs{:});
            else
                fprintf('Loading %s (%s)\n',matfile,IDStr(DATA));
                [Trials, Expts] = APlaySpkFile(matfile,'noerrs','nospikes', aargs{:});
            end
            if isempty(Expts)
                Expt = [];
                return;
            elseif length(Expts) > 1
                fprintf('WARNING: %d Expts in %s\nTrials:',length(Expts),exfile);
                for j = 1:length(Expts)
                    nt(j) = length(Expts{j}.Trials);
                end
                fprintf([sprintf(' %d',nt) '\n']);
                [a,b] = max(nt);
                Expt = Expts{b};
                if isfield(DATA,'toplevel') && DATA.interactive >= 0
                    setappdata(DATA.toplevel,'Expts',Expts);
                end
            else
            Expt = Expts{1};
            end
            if ~isfield(Expt.Header,'expname')
                Expt.Header.expname = Expt2Name(Expt);
            end
            
            Expt.Header.trialdur = sum([Expt.Trials.dur]);
            Expt = FillTrials(Expt,'ed');
            Expt = FillTrials(Expt,'st');
            [a,b] = fileparts(Expt.Header.Name);
            Expt.Header.title = [b '.' Expt.Header.expname];
            Expt.Header.preperiod = DATA.preperiod * 10000;
            Expt.Header.postperiod = DATA.postperiod * 10000;
            mytoc(ts);
            if ~isfield(Expt.Header,'ReadMethod')
                Expt.Header.ReadMethod = -1;
            end
            fprintf('Expt%d: %s %d trials\n',ei,Expt.Header.title,length(Expt.Trials));
        else
            Expt = [];
        end
        
function Expt = LoadExptA(DATA, exfile, ei)
    Expt = [];
    matfile = [];
    aargs = {};
    if DATA.usealltrials
        aargs = {aargs{:} 'usealltrials'};
    end
    %Loads expt from stadard smr file where one sme has all/many expts
    ei = floor(ei); %plain .mat files are not split just because FullV files are
    if exist(exfile,'file')
        ts = now;
        fprintf('LoadingA %s (%s)\n',exfile,IDStr(DATA));
        if ei == 0
            load(exfile);
        else
        [Trials, Expts] = APlaySpkFile(exfile,'noerrs','nospikes',aargs{:});
        if length(Expts) == 1
            Expt = Expts{1};
        elseif length(Expts) < ei 
            return;
        else
            Expt = Expts{ei};
        end
        end
        if ~isfield(Expt.Header,'expname')
            Expt.Header.expname = Expt2Name(Expt);
        end

        Expt.Header.trialdur = sum([Expt.Trials.dur]);
        Expt = FillTrials(Expt,'ed');
        Expt = FillTrials(Expt,'st');
        [a,b] = fileparts(Expt.Header.Name);
        Expt.Header.title = [b '.' Expt.Header.expname];
        fprintf('Expt%d: %s %d trials\n',ei,Expt.Header.title,length(Expt.Trials));
        mytoc(ts);
        if ei > 0
        if isfield(Trials,'DataType')
            Expt.Header.DataType = Trials.DataType;
        end
        Expt.Header.preperiod = DATA.preperiod * 10000;
        Expt.Header.postperiod = DATA.postperiod * 10000;
        if sum([Trials.ExptList.result] ==2) >= 1 && isfield(DATA,'toplevel')
            for j = 1:length(Expts)
            Expts{j}.Header.trialdur = sum([Expts{j}.Trials.dur]);
            Expts{j}.Header.DataType = Trials.DataType;
            Expts{j} = FillTrials(Expts{j},'ed');
            Expts{j} = FillTrials(Expts{j},'st');
            Expts{j}.Header.title = [b '.' Expt.Header.expname];
            end
            if isfield(DATA,'toplevel') && DATA.interactive >=0
            setappdata(DATA.toplevel,'Expts',Expts);
            end
        end
        end
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
            if length(nid) > 1 && length(id) > 1
                S.mu(1,:) = mean(X(id,:),1);
                S.mu(2,:) = mean(X(nid,:),1);
                S.Sigma(:,:,1) = cov(X(id,:));
                S.Sigma(:,:,2) = cov(X(nid,:));
                S.PComponents(1) = length(id)./size(X,1);
                S.PComponents(2) = length(nid)./size(X,1);
            end
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
    distance = gmdistance(G);
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
    
    
    
    
function [d, details]  = gmdprime(G, varargin)
%calcualte drpime between two Gaussians in gmdistribution fit        
    if ~isobject(G)
        d = NaN;
        details.d = [NaN NaN NaN];
        return;
    end
    nc =size(G.mu,1);
if size(G.mu,1) == 3
    %find three distances. one is allowed to be smaltt (splitting hash into
    %two. But one cluster must be distant from both of these. So take
    %lowest of top two = middle value
distance = mahal(G,G.mu);
    d(1) = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
    d(2) = sqrt(2./((1./distance(1,3))+(1./distance(3,1))));
    d(3) = sqrt(2./((1./distance(2,3))+(1./distance(3,2))));
    details.d =d;
    ds = sort(d);
    d  = d(2);  
elseif nc > 3
    distance = mahal(G,G.mu);
    for j = 1:nc
        for k = 1:j-1
            d(j,k) = sqrt(2./((1./distance(j,k))+(1./distance(k,j))));
            d(k,j) = d(j,k);
        end
    end
%need to return a scalar for typical use.
%put full matix in details
    details.d = d;
    d  = mean(d(d>0));  
%    d = sqrt(2./((1./distance(2,1))+(1./distance(1,2))));
else
    details.d = gmdistance(G);
    d = details.d;
end
    
function [x,y] = GetClusterXYData(DATA, p)
    C = GetClusterDef(DATA.cluster,DATA.currentcluster);
    if isfield(C,'space')
        usespace = C.space(1);
    else
        usespace = DATA.plottype;
    end
    if isempty(p) %Var/E plot
        y = DATA.spkvar(DATA.probe(1),DATA.uid)'./DATA.energy(1,DATA.uid)';
        x = DATA.energy(1,DATA.uid)';
    elseif usespace == 2 %voltage pairs
        AllV = GetDataStruct(DATA, 'AllV');
        x = AllV(p(1),p(2),DATA.uid);
        y = AllV(p(3),p(4),DATA.uid);
    elseif ismember(usespace, [3 4 7])
        x = DATA.TemplateScores(DATA.uid,p(1));
        y = DATA.TemplateScores(DATA.uid,p(2));
    else
        x = DATA.pcs(DATA.uid,p(1));
        y = DATA.pcs(DATA.uid,p(2));
    end
x = squeeze(x);
y = squeeze(y);


function x = Rprime(r)
    x = sqrt(r);
        
function Cut = PlotHistogram(DATA, E, varargin)
    
    plotdip = 0;
    checkdprimes = 0;
    plotgm = 0;
   replotpts = 0;
   quickmode.quick = 0;
    j = 1;
    while j <= length(varargin)
        if isfield(varargin{j}, 'quickest')
            quickmode = varargin{j};
            quickmode.quick = 1;
        elseif strncmpi(varargin{j},'plotdip',5)
            plotdip = 1;
        elseif strncmpi(varargin{j},'plotgmdetails',10)
            plotgm = 2;
        elseif strncmpi(varargin{j},'plotgm',5)
            plotgm = 1;
        elseif strncmpi(varargin{j},'quick',5)
            quickmode.quick = 1;
            plotgm = 0;
        elseif strncmpi(varargin{j},'fit1cut',6)
            quick = 2;
            plotgm = 0;
        end
        j = j+1;
    end
% E.shape == 2 means that DATA.xy{1} has already been rotated, so don't rotate here.
% if classify spike hasn't been called yet, this doesn't work.  Need to
% check for this....
if isempty(E)
    E = BoundaryFromCluster([],DATA.cluster, DATA.currentcluster);
end
if DATA.currentcluster == 1
    C = DATA.cluster;
else
    C = DATA.cluster.next{DATA.currentcluster-1};
end

if isfield(E,'pcplot')
    pcplot = E.pcplot;
elseif isfield(E,'space')
    pcplot = E.space(2:end);
end

if E.shape(1) == 2 || E.space(1) == 6
%    xy = DATA.xy{1};
    angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
    exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);
    crit = mean(exy(:,1));
    xy = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),angle);
else
    pos = E.pos;
[allx, ally] = GetClusterXYData(DATA, pcplot);
    if E.shape == 1
    angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
    xy = xyrotate(allx,ally,angle);
    else
        xy = cat(2,allx,ally);
        angle = 0;
    end
exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);
crit = mean(exy(:,1));
end

if DATA.interactive < 0
    [dip, mydip] = GMDip(xy, DATA.energy(1,:),'crit',crit,'label',DATA.idstr);
    Cut.gmdprime = mydip.gmdprime;
    return;
end


if ~isfield(E,'sign')
    E.sign = 0;
end
if E.shape(1) == 0
    clid = find(DATA.clst(DATA.uid) == DATA.currentcluster);
    nid = find(DATA.clst(DATA.uid) ~= DATA.currentcluster);
elseif E.sign > 0
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
cfig = SetFigure(DATA.tag.hist, DATA,DATA.watcharg{:});
subplot(2,1,2);
xid = [];
showgmfit = 0;
if isfield(E,'bestcl') && length(E.bestcl) == length(xy) && DATA.usegmcid
    fprintf('XY using GM clustering\n');
    nid = find(E.bestcl == 1);
    clid = find(E.bestcl == 2);
    xid=find(E.bestcl ==3);
    diffid = find(E.bestcl ==1 & xy(:,1) .* E.sign > crit .* E.sign);
    diffcid = find(E.bestcl ==2 & xy(:,1) .* E.sign < crit .* E.sign);
elseif E.shape == 0
    fprintf('XY using Ellipse\n');
    diffcid = [];
elseif isfield(E,'gmfit2dman') && DATA.usegmcid
    showgmfit = 1;
    fprintf('XY using GM 2D clustering\n');
    if E.shape == 0
        E.bestcl = cluster(E.gmfit2dman,DATA.xy{1});
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
   if isfield(E,'bestspace') && isfield(E,'mahal')
       str = sprintf(' Mahal %.2f (1D %.2f, 2D %.2f)',E.bestspace(2),E.mahal(4),E.mahal(1));
   elseif isfield(E,'bestspace')
       str = sprintf(' Mahal %.2f (%.2f)',a.mahal,E.bestspace(1));
   elseif isfield(E,'mahal')
       str = sprintf(' Mahal 1D %.2f, 2D %.2f',E.mahal(4),E.mahal(1));
   else
       str = sprintf(' Mahal %.2f (%.2f)',a.mahal);
   end
elseif isfield(E,'bestspace') 
    str = sprintf('Mahal %.2f 1D %.2f 2D %.2f',E.bestspace(1),E.mahal(4),E.mahal(1));
elseif isfield(E,'mahal')
    str = sprintf('Mahal %.2f.%.2f',E.mahal(1),E.mahal(4));
else
    str = [];
end

if E.shape == 2 || E.space(1) == 6
    title(sprintf('ND: %s%s',DATA.gmtypelabels{E.space(2)},str));
elseif isempty(pcplot)
    title(sprintf('Var-E %s',str));
else
    title(sprintf('%d: %dvs%d%s',DATA.plottype,pcplot(1),pcplot(2),str));
end



hdat = get(gcf,'UserData');
hdat.elmousept.pos(1) = crit;
hdat.elmousept.pos(3) = crit;
hdat.elmousept.pos(2) = exy(1,2);
hdat.elmousept.pos(4) = exy(2,2);
hdat.elmousept.shape = E.shape;
if isfield(E,'angle')
hdat.elmousept.angle = E.angle;
else
hdat.elmousept.angle = 0;
end
if isfield(E,'cluster')
hdat.elmousept.color = DATA.colors{E.cluster+1};
else
hdat.elmousept.color = 'r';
end
hdat.elmousept.down = 0;
%hdat.elmousept.h = DrawEllipse(hdat.elmousept);
if E.shape ~= 0 
    plot(exy(:,1),exy(:,2),'r-');
else
    hdat.elmousept.h = DrawEllipse(hdat.elmousept);
end

dp = (mean(xy(clid,1))-mean(xy(nid,1)))./sqrt(mean([var(xy(clid,1)) var(xy(nid,1))]));
hold off;
subplot(2,1,1);
hold off; 
if E.shape == 0 && isfield(E,'r')
[a,x] = hist(Rprime(E.r),500);
bar(x,a,1);
axis('tight');
hold on;
plot([1 1],get(gca,'ylim'),'r-');
else
[a,x] = hist(xy(:,1),500);
bar(x,a,1);
axis('tight');
hold on;
plot(exy(:,1),get(gca,'ylim'),'r-');
end
area = trapz(x,a);
if ~isfield(E,'quick')
    E.quick = 0;
end
%[dip, mydip] = FindDip(xy(:,1),DATA.energy(1,:),'eval',crit,'plot','gmix');
if plotgm == 2
    [dip, mydip] = GMDip(xy, DATA.energy(1,:),'plot','crit',crit,'label',DATA.idstr);
    if mydip.converged(1) == 0
        oldFindDip(xy(:,1),DATA.energy(1,:),'eval',crit,'plot',DATA.tag.dips,'gmix');
    end
elseif isfield(E,'gmfit1d') & strmatch('mu',fieldnames(E.gmfit1d),'exact') & quickmode.quick == 0
    if E.shape == 0 && isfield(E,'r')
        [dip, mydip] = GMDip(Rprime(E.r), E.gmfit1d,'crit',1,'label',DATA.idstr);
    elseif E.quick == 0  || DATA.quickcutmode.fit1cut; %if its a new cut and gmfit1d is old, don't show
        [dip, mydip] = GMDip(xy, E.gmfit1d,'crit',crit,'label',DATA.idstr);
    else
        dip(1:4) = NaN;
        mydip.type = 0;
        mydip.dipsize = 0;
        mydip.gmdprime = 0;
        mydip.cdipsize = 0;
    end
    mydip.sign = E.sign;
elseif quickmode.quick == 0 
    [dip, mydip] = GMDip(xy, DATA.energy(1,:),'crit',crit,'label',DATA.idstr);
else
    if E.shape == 1 %manual line, at least cal GM fit starting with crit
        [dip, mydip] = GMDip(xy, DATA.energy(1,:),'critonly',crit,'label',DATA.idstr);
    end
    dip(1:4) = 0;
    mydip.sign = 0;
    mydip.type = 0;
    mydip.dipsize = 0;
    mydip.gmdprime = 0;
    mydip.cdipsize = 0;
end
if E.sign == 0
    E.sign  = mydip.sign;
end
mycrit = dip;
if isfield(mydip,'gxy') %pdf of GM fit in 1D
    plot(mydip.gxy(:,1),sum(mydip.gxy(:,[2 3]),2).*area,'r');
    plot(mydip.gxy(:,1),mydip.gxy(:,2).*area,'g');
    plot(mydip.gxy(:,1),mydip.gxy(:,3).*area,'g');
end
if isfield(C,'fitdpparams')
    dpxy(:,1) = FitGauss(x,C.fitdpparams(1,:),'eval');
    dpxy(:,2) = FitGauss(x,C.fitdpparams(2,:),'eval');
    if length(C.fitdprime) > 3
        pscale = mean(diff(x))./C.fitdprime(4);
    else
    pscale = area./trapz(x,sum(dpxy,2));
    end
    
    plot(x,dpxy(:,1).*pscale,'m');
    plot(x,dpxy(:,2).*pscale,'m');
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
p = ProbeNumber(DATA);
if quickmode.quick == 0
    dip = HartigansDipTest(sort(xy(:,1))).*100;
    bii = BimodalCoeff(xy(:,1),1.5);
    t = sprintf('P%d%s Dip %.1f(%.1f,%.2f gm%.2f)',p,DATA.probelabel,dip,mydip.dipsize(1),bii,mydip.gmdprime);
else
    t = sprintf('P%d%s Dip (%.1f gm%.2f)',p,DATA.probelabel,mydip.dipsize(1),mydip.gmdprime);
end
if isfield(E,'fitdprime')
    t = [t sprintf(' G2 %.2f',E.fitdprime(1))];
end
title(t);
DATA.clid = clid;
DATA.nid = nid;
Cut = E;
Cut.area = area;
Cut.dip = mycrit;
Cut.angle = angle;
Cut.crit = [crit mycrit];
Cut.hdip = dip;
Cut.gmdprime = mydip.gmdprime;
Cut.mydip = [mydip.cdipsize mydip.dipsize];
Cut.space = [DATA.plottype pcplot];
Cut.shape = E.shape;
Cut.y = exy(:,2);
DATA.cluster.shape = 1;
hdat.cluster = Cut;
C = ClusterFromBoundary(E, Cut);
set(gcf,'UserData',hdat);
%DATA = ReplotPCs(DATA,E);
if C.shape == 1
newE = [];
newE = BoundaryFromCluster(newE, C, DATA.currentcluster);
%DrawLine(newE);
end

function ApplyLayout(DATA,varargin)
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'print',4)
            fprintf('Applying %s\n',DATA.layoutfile);
        end
        j = j+1;
    end
    
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

function [F, isnew] = SetFigure(lb, varargin)

    onoff = {'off' 'on'};
    
   PCDATA = [];
   toptag = 'PCs';
   DATA.interactive = 0;
   if length(varargin) && isstruct(varargin{1})
       DATA = varargin{1};
       toptag = DATA.tag.top;
       if length(varargin) > 1
           varargin = varargin(2:end);
       else
           varargin = {};
       end
   end
   if DATA.interactive < 0
       F = 0;
       isnew = 0;
       return;
   end
   if DATA.interactive == 0
       [F, isnew] = GetFigure(lb,'noforce',varargin{:});
   else
       [F, isnew] = GetFigure(lb,varargin{:});
   end
if isnew
    if strcmp(lb,toptag)
        DATA.toplevel = F;
        if exist(DATA.layoutfile) && DATA.readlayout
            load(DATA.layoutfile);
            setappdata(F,'Figpos',Figpos);
            if isfield(Figpos,lb)
                set(F,'Position',Figpos.(lb));
            end
        elseif DATA.readlayout
            fprintf('Cant read %s\n', DATA.layoutfile)
        end
        hm = uimenu(F,'Label','&Cluster','Tag','ClusterMenu');
        sm = uimenu(hm,'Label','Draw &Ellipse','Tag','EllipseMenu');
        uimenu(sm,'Label','Cluster &1','foregroundcolor',DATA.colors{2},'Callback',{@PCCluster, 'Ellipse1'});
        uimenu(sm,'Label','Cluster &2','foregroundcolor',DATA.colors{3},'Callback',{@PCCluster, 'Ellipse2'});
        uimenu(sm,'Label','Cluster &3','foregroundcolor',DATA.colors{4},'Callback',{@PCCluster, 'Ellipse3'});
        uimenu(sm,'Label','Cluster &4','foregroundcolor',DATA.colors{5},'Callback',{@PCCluster, 'Ellipse4'});
        sm = uimenu(hm,'Label','&Clear Cluster','Tag','EllipseMenu');
        uimenu(sm,'Label','Clear Cluster &2','foregroundcolor',DATA.colors{3},'Callback',{@PCCluster, 'Clear2'});
        b = uimenu(sm,'Label','Clear Cluster &3','foregroundcolor',DATA.colors{4},'Callback',{@PCCluster, 'Clear3'});
        uimenu(sm,'Label','Clear Cluster &4','Callback',{@PCCluster, 'Clear4'});
        uimenu(hm,'Label','&Use This Cluster','Callback',{@PCCluster, 34});
%        uimenu(hm,'Label','VarE Ellipse','Callback',{@PCCluster, 2});
        uimenu(hm,'Label','&Density','Callback',{@PCCluster, 4},'Tag','Density');
        uimenu(hm,'Label','Scale Density','Callback',{@MiscMenu,'scaledensity'},'Tag','ScaleDensity');
        sm = uimenu(hm,'Label','&Mark','Tag','MarkMenu');
        uimenu(sm,'Label','&Good','Callback',{@PCCluster, 'markgood'},'Tag','MarkGood');
        uimenu(sm,'Label','&Bad Probe','Callback',{@PCCluster, 'markbadprobe'},'Tag','MarkGood');
        uimenu(sm,'Label','Good &MU','Callback',{@PCCluster, 'markmugood'},'Tag','MarkGood');
%        uimenu(sm,'Label','&Good','Callback',{@PCCluster, 'markgood'},'Tag','MarkGood');
        uimenu(hm,'Label','&ReDraw','Callback',{@PCCluster, 3});
        uimenu(hm,'Label','&Template','Callback',{@PCCluster, 5});
        uimenu(hm,'Label','StdTemplate','Callback',{@PCCluster, 32});
        uimenu(hm,'Label','TemplateLine','Callback',{@PCCluster, 6});
 %       uimenu(hm,'Label','Var E Line','Callback',{@PCCluster, 7});
        uimenu(hm,'Label','PC &Line','Callback',{@PCCluster, 'Ellipse0'});
        sm = uimenu(hm,'Label','E&xtras','Tag','ClusterMenu');
        uimenu(sm,'Label','&Flip Sign','Callback',{@PCCluster, 'flipsign'});
        uimenu(sm,'Label','PC &Line2','Callback',{@PCCluster, 11});
        uimenu(sm,'Label','&Optimize Cluster','Callback',{@PCCluster, 18});
        uimenu(sm,'Label','&Iterate Fit','Callback',{@PCCluster, 'iteratefit'});
        uimenu(sm,'Label','&PlotExpt (once)','Callback',{@PlotResult, 'plotexpt'});
        uimenu(sm,'Label','&AutoMatic Cut','Callback',{@PCCluster, 19});
        uimenu(sm,'Label','&GaussMeans','Callback',{@PCCluster, 22});
        uimenu(sm,'Label','&Xcorr (clusters)','Callback',{@CalcXcorr, 'clusters'});
        uimenu(sm,'Label','Xcorr &Probes','Callback',{@CalcXcorr, 'probes'});
        uimenu(sm,'Label','PC on &dvdt','Callback',{@PCCluster, 'dvdtpc'});
        uimenu(sm,'Label','&Quick AutoCut','Callback',{@QuickAutoCut});
        sm = uimenu(hm,'Label','N-D','Tag','ClusterMenu');
        uimenu(sm,'Label','BestSpace','Callback',{@PCCluster, 23});
        uimenu(sm,'Label','BestSpace (scratch)','Callback',{@PCCluster, 26});
        uimenu(sm,'Label','BestSpace (2 cells)','Callback',{@PCCluster, '2cells'});
        uimenu(sm,'Label','BestSpace (3 cells)','Callback',{@PCCluster, '3cells'});
        uimenu(sm,'Label','BestSpace (4 cells)','Callback',{@PCCluster, '4cells'});
        uimenu(sm,'Label','AllSpaces (6 clusters)','CallBack',{@PCCluster, 'all6cells'});
        uimenu(sm,'Label','K means','Callback',{@PCCluster, 28});
        uimenu(sm,'Label','Use ND clid','Callback',{@PCCluster, 29});
        uimenu(sm,'Label','TemplateSpace','Callback',{@PCCluster, 'NDTemplate'});
        uimenu(sm,'Label','TemplateSpace N cell','Callback',{@PCCluster, 'NCellTemplate'});
        uimenu(sm,'Label','StdTemplateSpace','Callback',{@PCCluster, 'NDStdTemplate'});
        uimenu(sm,'Label','PCSpace','Callback',{@PCCluster, 31});
        uimenu(sm,'Label','PCSpace Multiple Cells','Callback',{@PCCluster, 'PCmultiple'});
        sm = uimenu(hm,'Label','PCs','Tag','ClusterMenu');
        uimenu(sm,'Label','Recalc without cell','Callback',{@PCCluster, 33});
        uimenu(sm,'Label','csd','Callback',{@SetOption},'tag','csd');
        uimenu(sm,'Label','dvdy','Callback',{@SetOption},'tag','dvdy');
        uimenu(sm,'Label','dvdy','Callback',{@SetOption},'tag','dvdt');
        
        sm = uimenu(hm,'Label','&Plot','Tag','ClusterMenu');
        uimenu(sm,'Label','&PCS','Callback',{@PCCluster, 9});
        uimenu(sm,'Label','&ADC','Callback',{@PCCluster, 10});
        uimenu(sm,'Label','&Templates','Callback',{@PCCluster, 14});
        uimenu(sm,'Label','Templates&2','Callback',{@PCCluster, 15});
        uimenu(sm,'Label','&dvdt','Callback',{@PCCluster, 20});
        uimenu(sm,'Label','&Other Bits','Callback',{@PCCluster, 'plotother'});
        uimenu(sm,'Label','Test Scores','Callback',{@SetPlot, 'testscore'});
        uimenu(sm,'Label','All Clusters','Callback',{@SetPlot, 'allspace'});
        if ~isempty(strfind(path,'FastICA'))
            uimenu(sm,'Label','&ICA','Callback',{@PCCluster, 'plotica'});
        end
        uimenu(sm,'Label','Cell &List','Callback',{@MiscMenu, 'plotcelllist'});
        uimenu(sm,'Label','&Cluster','Callback',{@PCCluster, 24});
        uimenu(sm,'Label','&Histogram','Callback',{@PCCluster, 26});
        uimenu(sm,'Label','use gmcids','Callback',{@SetOption, 'usegmcid'},'tag','usegmcid');
        uimenu(sm,'Label','&ISI Hist','Callback',{@PlotISI, 1});
        uimenu(sm,'Label','&ISI Sequence','Callback',{@PlotISI, 2});
        uimenu(sm,'Label','&XY Sequence','Callback',{@PlotMenu, 'xyseq'});
        uimenu(sm,'Label','&Rate Sequence','Callback',{@PlotMenu, 'rateseq'});
        uimenu(sm,'Label','&Expt','Callback',{@PlotResult, 1});
        uimenu(sm,'Label','&CrossCorr','Callback',{@PlotMenu, 'xcorrprobes'});
        sm = uimenu(hm,'Label','Axes:tight','Tag','ClusterMenu');
        uimenu(sm,'Label','100%','Callback',{@PCCluster, 100});
        uimenu(sm,'Label','99%','Callback',{@PCCluster, 99});
        sm = uimenu(hm,'Label','PlotClusters','Tag','PlotClusters');
        uimenu(sm,'Label','SpkTimes','Callback',{@PlotClusters, 1});
        uimenu(sm,'Label','xCorr','Callback',{@PlotClusters, 2});
        uimenu(sm,'Label','Quality','Callback',{@PlotClusters, 3});
        uimenu(hm,'Label','Delete','Callback',{@PCCluster, 17});
        uimenu(hm,'Label','Revert','Callback',{@PCCluster, 21});
        uimenu(hm,'Label','Save (&quick)','Callback',{@PCCluster, 12});
        uimenu(hm,'Label','&Save Spikes','Callback',{@PCCluster, 25});
        uimenu(hm,'Label','Save To Def&initions','Callback',{@PCCluster, 'savedef'});
        uimenu(hm,'Label','Quantify Quick Clusters','Callback',{@PCCluster, 'quantifyquick'});
        sm = uimenu(hm,'Label','&Window Layout/Setting','Tag','LayoutMenu');
        uimenu(sm,'Label','Windows to &Front','Callback',{@MiscMenu, 'tofront'});
        uimenu(sm,'Label','&Save Layout','Callback',{@MiscMenu, 'savelayout'});
        uimenu(sm,'Label','&Save Default Layout','Callback',{@MiscMenu, 'savedefaultlayout'});
        uimenu(sm,'Label','&Load','Callback',{@MiscMenu, 'loadlayout'});
        uimenu(sm,'Label','Load &Layout','Callback',{@MiscMenu, 'loadlayout'});
        uimenu(sm,'Label','Save &Config','Callback',{@MiscMenu, 'saveconfig'});
        uimenu(sm,'Label','Save Config as &Default','Callback',{@MiscMenu, 'savedefaultconfig'});
        uimenu(sm,'Label','Load &Configuration','Callback',{@MiscMenu, 'loadconfig'});
        
        
        hm = uimenu(F,'Label','&Probes','Tag','ProbeMenu');
        uimenu(hm,'Label','&Next','Tag','NextButton','Callback',{@ChangeProbe, 'next'});
        uimenu(hm,'Label','&QuickSave+Next','Tag','SaveNextButton','Callback',{@ChangeProbe, 'quicksave'});
        uimenu(hm,'Label','&Prev','Tag','NextButton','Callback',{@ChangeProbe, 'prev'});
        uimenu(hm,'Label','&Save+Next','Tag','SaveNextButton','Callback',{@ChangeProbe, 'save'});
        uimenu(hm,'Label','1','Callback',{@ProbeMenu, 1});
        uimenu(hm,'Label','3','Callback',{@ProbeMenu, 2});
        uimenu(hm,'Label','5','Callback',{@ProbeMenu, 3});
        uimenu(hm,'Label','dvdt','Callback',{@ProbeMenu, 4});
        uimenu(hm,'Label','csd','Callback',{@ProbeMenu, 5});
        uimenu(hm,'Label','dvdy','Callback',{@ProbeMenu, 'dvdy'});
        uimenu(hm,'Label','select','Callback',{@ProbeMenu, 'select'});

        sm =  uimenu(hm,'Label','S&witch To','Tag','ProbeSwitchMenu');
        pchars = ['1':'9' '0' 'a':'z'];
        for j = 1:DATA.allnprobes
            if j < 10
            uimenu(sm,'Label',['&' num2str(j)] ,'Callback',{@ChangeProbe, j});
            elseif j <= length(pchars)
                uimenu(sm,'Label',[num2str(j) ' (&'  pchars(j) ')'] ,'Callback',{@ChangeProbe, j});
            else
                uimenu(sm,'Label', num2str(j) ,'Callback',{@ChangeProbe, j});
            end
        end
        sm =  uimenu(hm,'Label','Switch &Mode','Tag','ProbeModeMenu');
        a(1) = uimenu(sm,'Label','usecluster (quicker)','CallBack',{@ProbeMenu, 'usecluster'},'tag','usecluster');
        a(2) = uimenu(sm,'Label','Reclassify','CallBack',{@ProbeMenu, 'reclassify'},'tag','reclassify');
        a(3) = uimenu(sm,'Label','Reapply','CallBack',{@ProbeMenu, 'reapply'},'tag','reapply');
        a(4) = uimenu(sm,'Label','AutoCut','CallBack',{@ProbeMenu, 'autocut'},'tag','autocut');  
        a(5) = uimenu(sm,'Label','Simple cut (quickest)','CallBack',{@ProbeMenu, 'simple'},'tag','simple');

        j = strmatch(DATA.probeswitchmode,{'usecluster' 'reclassify' 'reapply' 'autocut' 'simple'});
        if length(j) == 1
            set(a(j),'Checked','on');
        end
        hm = uimenu(F,'Label','&Options','Tag','OptionMenu');
        sm = uimenu(hm,'Label','When &Cutting','Tag','CutOptionMenu');
        uimenu(sm,'Label','&Quick','Tag','NextButton','Callback',{@OptionMenu, 'quickest'},'Tag','quickest');
        uimenu(sm,'Label','+&GM for 1d','Tag','NextButton','Callback',{@OptionMenu, 'fit1cut'},'Tag','fit1cut');
        uimenu(sm,'Label','+&2Gauss Fit','Tag','NextButton','Callback',{@OptionMenu, 'fit2gauss'},'Tag','fit2gauss');
        uimenu(sm,'Label','+&Drop Index','Tag','NextButton','Callback',{@OptionMenu, 'dropi'},'Tag','dropi');
        uimenu(sm,'Label','+&Clac Mean','Tag','NextButton','Callback',{@OptionMenu, 'calcmean'},'Tag','calcmean');
        uimenu(sm,'Label','Do all fits','Tag','NextButton','Callback',{@OptionMenu, 'fitallcut'},'Tag','fitallcut');
        uimenu(sm,'Label','Plot Spikes','Tag','NextButton','Callback',{@OptionMenu, 'plotspikes'},'Tag','plotspikes');
        uimenu(hm,'Label','Ne&w FullV','Tag','NextButton','Callback',{@OptionMenu, 'newfullv'},'Tag','NewFullV');
        uimenu(hm,'Label','&Next FullV','Tag','NextButton','Callback',{@OptionMenu, 'nextfullv'},'Tag','NextFullV');
        uimenu(hm,'Label','&Prev FullV','Tag','NextButton','Callback',{@OptionMenu, 'prevfullv'},'Tag','PrevFullV');
        sm = uimenu(hm,'Label','When->New FullV','Tag','SwitchOptionMenu');
        uimenu(sm,'Label','Summary Plot','Callback',{@OptionMenu, 'summary'},...
            'Tag','summary','checked',onoff{DATA.fullvswitchmode.summary+1});
        uimenu(sm,'Label','Check Clusters','Callback',{@OptionMenu, 'checkclusters'},...
            'Tag','checkclusters','checked',onoff{DATA.checkclusters+1});
        sm = uimenu(hm,'Label','Auto','Tag','AutoMenau');
        uimenu(sm,'Label','Advance Probe When Save','Callback',{@OptionMenu, 'autoopt'},'Tag','advanceprobe',...
            'checked',onoff{1+DATA.auto.advanceprobe});
        uimenu(sm,'Label','Advance Expt When Save','Callback',{@OptionMenu, 'autoopt'},'Tag','advanceexpt');
        uimenu(sm,'Label','Save when define ref','Callback',{@OptionMenu, 'autoopt'},'Tag','saveref');
        uimenu(sm,'Label','Xcorr for multiple clusters','Callback',{@OptionMenu, 'autoopt'},'Tag','checkxcorr');
        uimenu(sm,'Label','Re-use last  cluster for new expt','Callback',{@OptionMenu, 'autoopt'},'Tag','uselastcluster');
        uimenu(sm,'Label','Replot when change cluster','Callback',{@OptionMenu, 'autoopt'},'Tag','replotcluster');
        uimenu(sm,'Label','Replot when change cluster','Callback',{@OptionMenu, 'autoopt'},'Tag','checkcluster');

        uimenu(hm,'Label','S&ummary Plot','Callback',{@OptionMenu, 'plotsummary'});
        uimenu(hm,'Label','&Retrigger','Tag','NextButton','Callback',{@OptionMenu, 'retrigger'},'Tag','Retrigger');
        uimenu(hm,'Label','&Spool Spikes','Tag','SpoolSpikes','Callback',{@OptionMenu, 'spoolspikes'});
        uimenu(hm,'Label','Co&mment','Tag','Commnent','Callback',{@OptionMenu, 'comment'});
        uimenu(hm,'Label','Profiling','Tag','Profiling','Callback',{@OptionMenu, 'profiling'});
        
        set(F, 'KeyPressFcn',@PCKeyPressed);
        set(F,'CloseRequestFcn',@exitallv);
        PCDATA = DATA;
        set(F,'UserDATA',DATA);
    else
       it = findobj('type','figure','tag',toptag);
       PCDATA = get(it,'UserData');
    end    
    if strcmp(lb,'FullV')
%            set(F, 'WindowScrollWheelFcn',@ScrollV);
        hm = uimenu(F,'Label','Plot','Tag','ClusterMenu');
        uimenu(hm,'Label','1 Probe','Callback',{@OptionMenu, '1probefullv'});
        uimenu(hm,'Label','3 probes','Callback',{@OptionMenu,'3probefullv'});
        uimenu(hm,'Label','5 probes','Callback',{@OptionMenu, '5probefullv'});
        set(F,'UserData',DATA.toplevel);
    elseif strcmp(lb,'Hist')
        DATA.parentfigtag = toptag;
        set(F,'UserData',DATA);
        set(F, 'KeyPressFcn',@HistKeyPressed);
        set(F, 'WindowButtonDownFcn',@HistButtonPressed);
        set(F, 'WindowButtonMotionFcn',@HistButtonDragged);
        set(F, 'WindowButtonUpFcn',@HistButtonReleased);
        hm = uimenu(F,'Label','Plot','Tag','ClusterMenu');
        uimenu(hm,'Label','Dip Criterion','Callback',{@HistMenu, 1});
        uimenu(hm,'Label','mahal/angle','Callback',{@HistMenu, 2});
        uimenu(hm,'Label','Hist Replot','Callback',{@HistMenu, 3});
        uimenu(hm,'Label','Flip Criterion','Callback',{@HistMenu, 4});
    elseif strcmp(lb,'ExptFig')
        hm = uimenu(F,'Label','Fit','Tag','ExptFigMenu');
        uimenu(hm,'Label','Fit','Callback',{@ExptFigMenu, 'fit'});
        X.toplevel = PCDATA.toplevel;
        set(F,'UserData',X);
    elseif strcmp(lb,'Clusters')
        DATA.parentfigtag =toptag;
        hm = uimenu(F,'Label','Plot','Tag','ClusterMenu');
        set(F,'UserData',DATA);
        uimenu(hm,'Label','Times','Callback',{@PlotCluster, 1});
        uimenu(hm,'Label','xcorr','Callback',{@PlotCluster, 2});
        uimenu(hm,'Label','xcorr adjacent','Callback',{@PlotCluster, 3});
        uimenu(hm,'Label','xcorr 2sep','Callback',{@PlotCluster, 4});
        uimenu(hm,'Label','Quality Scatter','Callback',{@PlotCluster, 5});
        uimenu(hm,'Label','QualityProbes','Callback',{@PlotCluster, 6});
    elseif strcmp(lb,'Spikes')
        DATA.parentfigtag = toptag;
        set(F,'UserData',DATA);
        set(F, 'WindowScrollWheelFcn',@ScrollSpikes);
        set(F, 'KeyPressFcn',@KeyPressed);

        hm = uimenu(F,'Label','Scroll','Tag','ClusterMenu','ButtonDownFcn',{@SpikeDraw, 'menu'});
        uimenu(hm,'Label','Spool Spikes','Callback',{@SpikeDraw, 3});
        uimenu(hm,'Label','Spool by Trial','Callback',{@SpikeDraw, 10},'Tag','PlotByTrial');
        uimenu(hm,'Label','More Spikes','Callback',{@SpikeDraw, 1});
        uimenu(hm,'Label','Fewer Spikes','Callback',{@SpikeDraw, 2});
        uimenu(hm,'Label','Spool All Probes','Callback',{@SpikeDraw, 13});
        uimenu(hm,'Label','QuickSpk All Probes','Callback',{@SpikeDraw, 'allquickspks'});
        uimenu(hm,'Label','Mean All Probes','Callback',{@SpikeDraw, 'allmeans'});
        uimenu(hm,'Label','dVdt','Callback',{@SpikeDraw, 4});
        uimenu(hm,'Label','csd','Callback',{@SpikeDraw, 5});
        uimenu(hm,'Label','dvdy','Callback',{@SpikeDraw, 'dvdy'});
        uimenu(hm,'Label','Subtact Trigger','Callback',{@SpikeDraw, 'subtrigger'});
        uimenu(hm,'Label','Subtact mean','Callback',{@SpikeDraw, 6});
        uimenu(hm,'Label','Subtract peak','Callback',{@SpikeDraw, 7});
        uimenu(hm,'Label','Subtract Min','Callback',{@SpikeDraw, 8});
        uimenu(hm,'Label','Main Probe Only','Callback',{@SpikeDraw, 9});        
        uimenu(hm,'Label','Exclude Later Spikes','Callback',{@SpikeDraw, 11});        
        uimenu(hm,'Label','Exclude Earlier Spikes','Callback',{@SpikeDraw, 12});        
        uimenu(hm,'Label','Exclude Selected Trials','Callback',{@SpikeDraw, 'excludetrials'});        
        uimenu(hm,'Label','Use All Spikes','Callback',{@SpikeDraw, 14});        
        uimenu(hm,'Label','xCorr Selected','Callback',{@SpikeDraw, 15});        
        uimenu(hm,'Label','xCorr Adjacent','Callback',{@SpikeDraw, 'xcorradj'});        
        uimenu(hm,'Label','xCorr Aall','Callback',{@SpikeDraw, 'xcorrall'});        
        uimenu(hm,'Label','Spoolwith mean','Callback',{@SpikeDraw, 'spoolwithmean'});        

        hm = uimenu(F,'Label','Options','Tag','OptionMenu');
%        sm = uimenu(hm,'Label','When Cutting','Tag','CutOptionMenu');
        sm = hm;
        uimenu(sm,'Label','Show Full V','Tag','FullVToggle','Callback',{@OptionMenu, 'showfullv'});
        uimenu(sm,'Label','Include Pre/Post period','Tag','PrePostToggle','Callback',{@OptionMenu, 'includeprepost'});
        hm= uimenu(sm,'Label','Define ADC sample','Tag','ADCSetButton');
        uimenu(hm,'Label','#1 (Main probe)','Tag','ADCSetButton1' , 'Callback', {@SetADC ,1});
        uimenu(hm,'Label','#2 (Main probe)','Tag','ADCSetButton2', 'Callback', {@SetADC ,2});
        uimenu(hm,'Label','#3 (Main probe)','Tag','ADCSetButton3', 'Callback', {@SetADC ,3});
        uimenu(hm,'Label','#1 (Main probe)','Tag','ADCSetButton1' , 'Callback', {@SetADC ,4});
        uimenu(hm,'Label','#2 (Main probe)','Tag','ADCSetButton2', 'Callback', {@SetADC ,5});
        uimenu(sm,'Label','Keep','Tag','KeepButton','Callback',{@OptionMenu, 'keepspikes'});
        uimenu(sm,'Label','Keep Means','Tag','KeepMeans','Callback',{@OptionMenu, 'keepmeanspikes'});
        
        
        bp = [0.80 0.95 0.1 0.05];
        uicontrol(F,'style','check','string','stop','Units','Normalized','Position',bp,'Tag','StopSpool',...
            'value',0);
        bp = [0.9 0.05 0.1 0.95];

        if isfield(PCDATA,'Expt') && isfield(PCDATA.Expt,'Trials')
        h = uicontrol(gcf,'Style', 'list',...
        'String', num2str([PCDATA.Expt.Trials.Trial]'), 'Tag', 'ChooseTrial','Units','norm', 'Position', bp,'value',1,...
        'Max',3,'Min',1,'Callback',@SelectTrial);
        else
            disp('Missing PCS Data');
        end
%        uicontrol(F,'style','pop','string','1|2|3','Units','Normalized','Position',bp,'Tag','ChooseTrial',...
%            'Callback', @SelectTrial);
    elseif strcmp(lb,'TemplateScores')
        DATA.parentfigtag = toptag;
        set(F,'UserData',DATA);
        hm = uimenu(F,'Label','Scroll','Tag','Classify');
        uimenu(hm,'Label','Line','Callback',{@PCCluster, 6});
    elseif isfield(PCDATA,'toplevel');
        set(F,'UserData',PCDATA.toplevel);
    end
    if ~isempty(PCDATA)
        SetFigPos(PCDATA,lb);
    end
end

function DATA = LoadCellFile(DATA)

    cellfile = [DATA.name '/CellList.mat'];
    if exist(cellfile,'file')
        load(cellfile);
        DATA.CellList = CellList;
        DATA.CellDetails = CellDetails;
        DATA.CellChanges = CellChanges;
        DATA = LoadComments(DATA);
        if isfield(DATA,'tagged')
            [a,b] = find(DATA.tagged > 0);
            id = find(DATA.CellDetails.exptids(a) == DATA.exptno);
            if ~isempty(id)
                DATA.TaggedProbes = DATA.tagged(a(id(1)),:);
                for j = 1:length(id)
                    fprintf('Expt %d Probe %d Tagged %d\n',DATA.exptno,b(id),DATA.TaggedProbes(b(id(j))));
                end
                ShowTaggedProbes(DATA);
            end
        end
        nc = max(DATA.CellList(:));
        n = length(DATA.comparecell);
        if nc > n
            DATA.comparecell(n+1,nc) = 0;
        end
    end
        

    
    
function [true, cellid] = isacell(DATA, ei, p)
% [true, cellid] = isacell(DATA, ei, p)
% true is any clusters on this expt(ei_ and probe(p) is in the cell list
    cellid = 0;
    if ~isfield(DATA,'CellList') || ei > size(DATA.CellList,1)
        true = 0;
        return;
    end
    if ~isfield(DATA.CellDetails,'exptids')
        fprintf('CellList missing exptids\n');
        true = 0;
        return;
    end
    row  = find(DATA.CellDetails.exptids == ei);

    true = sum(DATA.CellList(row,p,:),3) > 0;
    if true
        cellid = squeeze(DATA.CellList(row,p,:));
  %      cellid = cellid(cellid>0);
    end

    
function AddCellMenu(DATA)
    if DATA.interactive < 0
        return;
    end
    hm = findobj(DATA.toplevel,'Tag','ProbeMenu');
    if ~isfield(DATA,'CellList')
        return;
    end
    sm = findobj(hm, 'Tag','CellSwitchMenu');
    delete(sm);
    sm =  uimenu(hm,'Label','&Cell','Tag','CellSwitchMenu');
    cellids = unique(DATA.CellList(:));
    cellids = cellids(cellids > 0);
    cellids = reshape(cellids,1,length(cellids));
    pchars = ['1':'9' '0' 'a':'z'];
    for j = cellids
        if j < 10
            uimenu(sm,'Label',['&' num2str(j)] ,'Callback',{@ChangeCell, j});
        elseif j <= length(pchars)
            uimenu(sm,'Label',[num2str(j) ' (&'  pchars(j) ')'] ,'Callback',{@ChangeCell, j});
        else
            uimenu(sm,'Label', num2str(j) ,'Callback',{@ChangeCell, j});
        end
    end

    
    
function exitallv(src, evnt)
    DATA = get(src,'UserData');
    if isfield(DATA,'tag')
        f = fields(DATA.tag);
        for j = 1:length(f)
            if ~strcmp(f{j},'top')
                CloseTag(DATA.tag.(f{j}));
            end
        end
    end
    try
        fclose(DATA.logfid);
    end
    delete(src);

    
function vpts = SetVsamples(vsmps, probe, np, nv)
%DATA.vsmps = [20 6 15 11 30 20];

    vpts = [0 vsmps(1) 0 vsmps(2); ...
            0 vsmps(1) 0 vsmps(3); ...
            0 vsmps(1) 0 vsmps(4); ...
            0 vsmps(1) 1 vsmps(5); ...
            0 vsmps(1) -1 vsmps(3); ...
            0 vsmps(1)  -1 vsmps(5); ...
            0 vsmps(4) 0 vsmps(5); ...
            0 vsmps(1) 1 vsmps(3)];
        
vpts(:,1) = vpts(:,1) + probe(1);
vpts(:,3) = vpts(:,3) + probe(1);
id = find(vpts(:,1) < 1);
vpts(id,1) = 1;
id = find(vpts(:,1) > np);
vpts(id,1) = np;
id = find(vpts(:,3) < 1);
vpts(id,3) = 1;
id = find(vpts(:,3) > np);
vpts(id,3) = np;
id = find(vpts(:,2) > nv);
vpts(id,2) = nv;
id = find(vpts(:,4) > nv);
vpts(id,4) = nv;


    
function SetADC(a,b,fcn)
    DATA = GetDataFromFig(a);
    DATA.adcmousept.down = 0;
    DATA.adcmousept.h = -1;
        
    if fcn == 10
        DATA.vsmps(DATA.setadcpos) = DATA.adcmousept.x;
        DATA.vpts = SetVsamples(DATA.vsmps,DATA.probe, DATA.nprobes, DATA.nvpts);
        set(gcf,'ButtonDownFcn',{@ShowADCPos, 0});
        set(gca,'ButtonDownFcn',{@ShowADCPos, 0});
        set(gcf,'WindowButtonMotionFcn',{@ShowADCPos, 0});
        set(gcf,'WindowButtonUpFcn',{@ShowADCPos, 0});
        DATA.plottype = 2;
        ReplotPCs(DATA,[]);
    else
        DATA.setadcpos = fcn;
        set(gcf,'ButtonDownFcn',{@ShowADCPos, 1});
        set(gca,'ButtonDownFcn',{@ShowADCPos, 1});
        set(gcf,'WindowButtonMotionFcn',{@ShowADCPos, 2});
        set(gcf,'WindowButtonUpFcn',{@ShowADCPos, 3});
    end
    set(DATA.toplevel,'UserData',DATA);
    
function ShowADCPos(src, data, type)
    
 if type == 0
     return;
 end
DATA = GetDataFromFig(src);

DATA.ts = now;
start = get(gca,'CurrentPoint');
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
    x = round(start(1,1));
    if type ==3 
        if  DATA.adcmousept.down &&  ishandle(DATA.adcmousept.h)
            delete(DATA.adcmousept.h);
            SetADC(src, data, 10);
        end
    else
    hold on;
    if type == 2 && ishandle(DATA.adcmousept.h)  &&  DATA.adcmousept.down == 1
        set(DATA.adcmousept.h,'xdata',[x x]);
        DATA.adcmousept.x = x;
        set(DATA.toplevel,'UserData',DATA);
    elseif type == 1
        DATA.adcmousept.down = 1;
        DATA.adcmousept.h = plot([x x],get(gca,'ylim'),'k:');
        DATA.adcmousept.x = x;
        set(DATA.toplevel,'UserData',DATA);
    end

    end
function SetFigPos(DATA, tag)
    Figpos = getappdata(DATA.toplevel,'Figpos');
    if isfield(Figpos,tag)
        it = findobj('type','figure','Tag',tag);
        if length(it) == 1
            set(it,'position',Figpos.(tag));
        end
    end
        
        
function CalcXcorr(a,b,fcn)
DATA = GetDataFromFig(a);
DataClusters = getappdata(DATA.toplevel,'Clusters');
colors = mycolors;

if strcmp(fcn,'clusters')
ta = DATA.t(DATA.clst ==2);        
tb = DATA.t(DATA.clst ==3);
xc = xcorrtimes(ta,tb);
SetFigure(DATA.tag.xcorr, DATA,'front');
subplot(1,1,1);
plot([-200:200],xc);
elseif strcmp(fcn,'probes')
    p = ProbeNumber(DATA);
    SetFigure(DATA.tag.xcorr, DATA,'front');
    chspk = setdiff(DATA.chspk,DATA.probe(1));
    ta = DATA.t(DATA.clst ==DATA.currentcluster+1);
    hold off;
    for j = 1:length(chspk)
        tb = DataClusters{chspk(j)}.times;
        xc = xcorrtimes(ta,tb);
        plot([-200:200],xc./mean(xc),'color',colors{j});
        hold on;
    end
    legend(num2str(chspk'));
end

function FullVKeyPressed(src, ks)

if sum(strcmp(ks.Key,{'leftarrow' 'rightarrow'})) 
    DATA = GetDataFromFig(src);
    VT = get(gcf,'UserData');
    if strcmp(ks.Key,'leftarrow')
        VT.spkid = VT.spkid-1;
    elseif strcmp(ks.Key,'rightarrow')
        VT.spkid = VT.spkid+1;
    end
    if VT.spkid > 0 && VT.spkid <= length(VT.ids{2})
    t = VT.t(VT.ids{2}(VT.spkid));
    set(gca,'xlim',[t-0.001 t + 0.001]);
    end
    set(gcf,'UserData',VT);
    if VT.spkid > length(VT.ids{2})
        DATA = get(VT.toplevel,'UserData');
        [DATA.currenttrial, DATA.spklst] = PlotTrialSpikes(DATA, DATA.currenttrial+1);
        set(DATA.toplevel,'UserData',DATA);
    end
end

function PCKeyPressed(src, ks)

DATA = GetDataFromFig(src);
spk = [0 0];
changed = 1;


if strmatch(ks.Key,'p')
    if ~strcmp('alt',ks.Modifier)
    PCCluster(DATA,[],9); %plot PCS
    end
elseif strmatch(ks.Key,{'1' '2' '3' '4'})
    PCCluster(DATA.toplevel,[],['Ellipse' ks.Key]); %set ellipse
elseif strmatch(ks.Key,'a') 
    PCCluster(DATA,[],10); %plot ADCs
elseif strmatch(ks.Key,'d') 
    PCCluster(DATA,[],4); %plot Density
elseif strmatch(ks.Key,'e') 
    PCCluster(DATA,[],1); %ellipse
elseif strmatch(ks.Key,'h') %replot histogram
    PlotHistogram(DATA,[],'plotgmdetails');
    figure(DATA.toplevel);
elseif strmatch(ks.Key,'l') 
    PCCluster(DATA,[],8); %cut with line
elseif strmatch(ks.Key,'n')
    ChangeProbe(DATA,[],'next')
%    PCCluster(DATA,[],13); %Next Probe
    if strmatch('shift',ks.Modifier)
    end
elseif strmatch(ks.Key,'q')
    PCCluster(DATA,[],12); %quick save
elseif strmatch(ks.Key,'r')
    ReplotPCs(DATA,[]);
elseif strmatch(ks.Key,'s') 
    if strmatch('shift',ks.Modifier)
        PCCluster(DATA,[],12); %save cluster
    else
        PCCluster(DATA,[],25); %save spikes+cluster
    end
elseif strmatch(ks.Key,'t') 
    PCCluster(DATA,[],5); %plot Templates
elseif strcmp(ks.Character,'V') 
    OptionMenu(src, ks,  'newfullv');
elseif strmatch(ks.Key,'v') 
    SetFigure(DATA.tag.spikes, DATA,DATA.watcharg{:},'front');
    SpoolSpikes(DATA); %ellipse
    figure(DATA.toplevel);
elseif strmatch(ks.Key,'w') %sWitch probe
    ProbeMenu(src, ks,'select');
elseif strmatch(ks.Key,'x') 
    PCCluster(DATA,[],20); %plot extra adc
elseif strcmp(ks.Key,'leftarrow')
    RotateCluster(DATA, 0.05);
elseif strcmp(ks.Key,'rightarrow')
    RotateCluster(DATA, -0.05);
else
    changed = 0;
end

function SelectTrial(src, b)
DATA = GetDataFromFig(src);
id = get(src,'value');
PlotOneTrial(DATA,id)

function PlotOneTrial(DATA,id)
[DATA.currenttrial, DATA.spklst] = PlotTrialSpikes(DATA,id,'showall'); 
T = DATA.Expt.Trials(DATA.currenttrial);
fprintf('Trial %d: %.2f-%.2f\n',DATA.currenttrial,T.Start(1)./10000,T.End(end)./10000);
set(DATA.toplevel,'UserData',DATA);

    
function SetTrialList(DATA)
    if DATA.interactive < 0
        return;
    end
    it = findobj('Tag','ChooseTrial');
    if ~isempty(it) & isfield(DATA.Expt,'Trials')
        tlist = [DATA.Expt.Trials.Trial];
        idlist = [DATA.Expt.Trials.id];
        id = find(ismember(idlist,DATA.excludetrialids));
        tlist(id) = tlist(id) .* -1;
        set(it,'string',sprintf('%d|',tlist),'value',1);
    end
        
        

function HistKeyPressed(src, ks)

DATA = GetDataFromFig(src);
spk = [0 0];
changed = 1;

if ~isfield(DATA,'cboundary')
    return;
end
if strmatch(ks.Key,'rightarrow')
    if DATA.cboundary.shape == 2
        DATA.xy{1} = xyrotate(DATA.xy{1}(:,1),DATA.xy{1}(:,2),5 * pi/180);
    else
        DATA.cboundary.pos = RotateLine(DATA.cboundary.pos,5 * pi/180);
        if isfield(DATA.cboundary,'axis')
            axes(DATA.cboundary.axis);
        end
    end
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plotgmdetails');
    SetFigure(DATA.tag.hist, DATA,DATA.watcharg{:});
elseif strmatch(ks.Key,'leftarrow')
    if DATA.cboundary.shape == 2
        DATA.xy{1} = xyrotate(DATA.xy{1}(:,1),DATA.xy{1}(:,2),-5 * pi/180);
    else
        DATA.cboundary.pos = RotateLine(DATA.cboundary.pos,-5 * pi/180);
        if isfield(DATA.cboundary,'axis')
            axis(DATA.cboundary.axis);
        end
    end
%    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plotgmdetails');
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary);
    SetFigure(DATA.tag.hist, DATA,DATA.watcharg{:});
elseif strmatch(ks.Key,'d')
    DATA.cboundary = PlotHistogram(DATA,DATA.cboundary,'plotgmdetails');
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

if ~isfield(DATA,'currenttrial') %keypress before ready.
    return;
end
if isempty(DATA.spklst)
    DATA.spklst = 1:100;
end
w = DATA.spksperview;
nt = DATA.currenttrial;

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
if DATA.plotspk.bytrial || DATA.plotspk.allprobes
    DATA.currenttrial = PlotTrialSpikes(DATA,nt,'setcontext');
    FinishSpikePlot(DATA);
    set(DATA.toplevel,'UserData',DATA);
    return;
end
if spk(2) > 0 && spk(1) < DATA.nevents
    spklst = spk(1):spk(2);
    DATA.spklst = spklst;
    PlotSpikes(DATA,DATA.spklst,'fixy');
    set(DATA.toplevel,'UserData',DATA);
end


function PlotQuickSpikes(DATA, nspk, varargin)
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    npts = length(DATA.spts);
    Vall = getappdata(DATA.toplevel,'Vall');
    probes = 1:size(Vall.V,1);
    DATA.handles = ones(DATA.nprobes,2).*-1;
    subplot('position',[0.01 0.01 0.98 0.98]);
    hold off;
    for p = probes
        step = ceil(length(DATA.trigtimes{p})./nspk);
        spklst = 1:step:length(DATA.trigtimes{p});
        xoff = floor((p-1)/6) .* npts;
        yoff = (rem(p-1,6)) * DATA.vsep;
        PlotProbeSpikes(DATA, Vall, p, DATA.trigtimes{p}(spklst),[npts 8],[xoff yoff]);
        text(xoff,yoff+1,sprintf('%d:%.1f',p,DataClusters{p}.mahal(1)),'fontsize',DATA.gui.fontsize(1));
    end
    set(gca,'ylim',[-2 5 .* DATA.vsep+2],'xtick',[],'ytick',[]);
    set(gca,'xlim',[-8 length(DATA.spts)-8+xoff]);
    set(gca,'ButtonDownFcn',@SpikeButtonPressed);


function [nt, spklst] = PlotTrialSpikes(DATA, nt, varargin)

    
    
spklst = [];
if length(nt) > 1
    for j = 1:length(nt)
        if j == length(nt)
            [k, spklst] = PlotTrialSpikes(DATA, nt(j), varargin{:},'setcontext');
        else
            [k, spklst] = PlotTrialSpikes(DATA, nt(j), varargin{:});
        end
        drawnow;
    end
    nt = k;
    FinishSpikePlot(DATA);
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

args = {};
probes = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'probes',5)
        j = j+1;
        probes = varargin{j};
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

if DATA.plotspk.allprobes
    npts = length(DATA.spts);
    Vall = getappdata(DATA.toplevel,'Vall');
    if isempty(probes)
        probes = 1:size(Vall.V,1);
    end
    T = DATA.Expt.Trials(nt);
    for p = probes
        spklst = find(Vall.t(DATA.trigtimes{p}) .*10000 > DATA.Expt.Trials(nt).Start(1) & Vall.t(DATA.trigtimes{p}) .*10000 < DATA.Expt.Trials(nt).End(end));
        xoff = floor((p-1)/6) .* npts;
        yoff = (rem(p-1,6)) * DATA.vsep;
        PlotProbeSpikes(DATA, Vall, p, DATA.trigtimes{p}(spklst),[npts 8],[xoff yoff]);
    end
    title(sprintf('Trial %d (du%.2f): ed%.2f,st%.0f',DATA.Expt.Trials(nt).Trial, T.dur./10000,T.ed,T.st));
    drawnow;
else

    t(1) = DATA.Expt.Trials(nt).Start(1);
    t(2) = DATA.Expt.Trials(nt).End(end);
    if DATA.plotspk.includeprepost
        t(1) = t(1) - DATA.preperiod.*10000;
        t(2) = t(2) + DATA.postperiod.*10000;
    end
    spklst = find(DATA.t .*10000 > t(1) & DATA.t.*10000 < t(2));
    DATA.spklst = spklst;
    DATA.currenttrial = nt;
    PlotSpikes(DATA,DATA.spklst,'fixy',args{:});
    D = get(gcf,'UserData');
    D.plottype = 'TrialSpikes';
    set(gcf,'UserData',D);

end
DATA.currenttrial = nt;
% set(DATA.toplevel,'UserData',DATA);


function PlotCluster(a,b, mode)
DATA = GetDataFromFig(a);
DataClusters = getappdata(DATA.toplevel,'Clusters');

if mode == 1
    PlotSpikeTimes(DataClusters);
elseif mode == 2
    PlotSpikeTimes(DataClusters,'xcorr');
elseif mode == 3
    PlotSpikeTimes(DataClusters,'xcorrstep',1);
elseif mode == 4
    PlotSpikeTimes(DataClusters,'xcorrstep',2);
elseif mode == 5
    PlotSpikeTimes(DataClusters,'dips');
elseif mode == 6
    PlotSpikeTimes(DataClusters,'probequality');
end

function HistMenu(a,b, mode)
DATA = GetDataFromFig(a);
if mode == 1
%    dip = FindDip(DATA.xy{1}(:,1),DATA.energy(1,:),'plot');
    dip = GMDip(DATA.xy{1},DATA.energy(1,:),'plot','label',DATA.idstr);
elseif mode == 2
    [a,b,c] = BestGMAngle(DATA.xy{1}(:,1),DATA.xy{1}(:,2));
    SetFigure(DATA.tag.dips, DATA,DATA.watcharg{:});
    hold off; 
    plot(c.angles,c.mahal);
    title(sprintf('Mahal 1D %.3f at %.2fdeg , 2D %.3f',b, a .*180/pi,c.mahal2d));
elseif mode == 3
    H = PlotHistogram(DATA, [],'plotgmdetails');
    dip = MyDip(DATA.xy{1}(:,1));
    fprintf('Dip size%s\n',sprintf(' %.2f',dip.dipsize));
    sc = H.area./trapz(dip.x,dip.finey);
    plot(dip.x,dip.finey.*sc,'c');
    plot([dip.x(dip.dip(3)) dip.x(dip.dip(3))],get(gca,'ylim'),'c');
elseif mode == 4
    DATA.cluster.crit(1) = -DATA.cluster.crit(1);
%    PlotHistogram(DATA, []);
    set(DATA.toplevel,'UserData',DATA);
    PCCluster(a,b,34);
end
    
function ExptFigMenu(a,b, mode)
DATA = GetDataFromFig(a);
onoff = {'off' 'on'};

if strcmp(mode,'fit')
%    dip = FindDip(DATA.xy{1}(:,1),DATA.energy(1,:),'plot');
    DATA.plot.exptfit = ~DATA.plot.exptfit;
    set(a,'Checked',onoff{DATA.plot.exptfit+1});
    if DATA.plot.exptfit
        [Expt, res] = PlotExptCounts(DATA);
    end
    set(DATA.toplevel,'UserData',DATA);
elseif mode == 2
    [a,b,c] = BestGMAngle(DATA.xy{1}(:,1),DATA.xy{1}(:,2));
    SetFigure(DATA.tag.dips, DATA,DATA.watcharg{:});
    hold off; 
    plot(c.angles,c.mahal);
    title(sprintf('Mahal 1D %.3f at %.2fdeg , 2D %.3f',b, a .*180/pi,c.mahal2d));
elseif mode == 3
    PlotHistogram(DATA, [],'plotgmdetails');
elseif mode == 4
    DATA.cluster.crit(1) = -DATA.cluster.crit(1);
%    PlotHistogram(DATA, []);
    set(DATA.toplevel,'UserData',DATA);
    PCCluster(a,b,34);
end



function PlotXcorr(a,b, pa, pb)
DATA = GetDataFromFig(a);
if strcmp(pa,'zoom')
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

DataClusters = getappdata(DATA.toplevel,'Clusters');
mysubplot(2,2,2);
P = DataClusters{pa};
Q = DataClusters{pb};
[xc, details]  = xcorrtimes(P.times,Q.times);
plot(details.xpts,xc,'k','linewidth',2);
axis('tight');
yl = get(gca,'ylim');
xl = get(gca,'xlim');
[a,b] = max(xc);
text(xl(2),yl(2),sprintf('%.0fms (%d+-%.1f) %d %d',details.xpts(b).*1000,prctile(xc,50),std(xc),pa,pb),'horizontalalignment','right','verticalalignment','top');
line([0 0],yl,'linestyle','--');
set(gca,'buttondownfcn', {@PlotXcorr, 'zoom', 2});

function SpikeDraw(a,b, mode)
DATA = GetDataFromFig(a);
DataClusters = getappdata(DATA.toplevel,'Clusters');

onoff = {'off' 'on'};
redraw = 0;
if mode == 2
    DATA.spksperview = round(DATA.spksperview/2);
elseif mode == 1
    DATA.spksperview = DATA.spksperview *2;
elseif mode == 3
   DATA.plotspk.allprobes = 0;
   SpoolSpikes(DATA);
elseif mode == 4
    DATA.plotdvdt = ~DATA.plotdvdt;
    set(a,'Checked',onoff{DATA.plotdvdt+1});
    PlotMeanSpike(DATA);
elseif mode == 5
    if DATA.plotcsd == 1
        DATA.plotcsd = 0;
    else
        DATA.plotcsd = 1;
    end
    set(a,'Checked',onoff{DATA.plotcsd+1});
    PlotMeanSpike(DATA);
elseif strcmp(mode,'dvdy')
    DATA.plotcsd =  ~(DATA.plotcsd == 2)*2;
    set(a,'Checked',onoff{DATA.plotcsd/2+1});
    redraw = 1;
    PlotMeanSpike(DATA);
elseif strcmp(mode,'subtrigger')
    DATA.plotspk.subtrigger = ~DATA.plotspk.subtrigger;
    set(a,'Checked',onoff{DATA.plotspk.subtrigger+1});
    redraw = 1;
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
    set(a,'Checked',onoff{DATA.plotspk.submin+1});
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
    DATA = SetPCs(DATA,1);
elseif strcmp(mode,'spoolwithmean')
    DATA.plotspk.showmean = ~DATA.plotspk.showmean;
    set(a,'Checked',onoff{DATA.plotspk.showmean+1});
    redraw = 1;
elseif strcmp(mode,'allmeans')
    PlotAllMeans(DATA);
elseif strcmp(mode,'menu') % hit menu itself
    bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})

elseif strcmp(mode,'excludetrials')
    f = findobj('type','figure','Tag',DATA.tag.spikes);
    it = findobj(f, 'Tag','ChooseTrial');
    id = get(it,'value');
    DATA = ExcludeTrials(DATA,id, 1);
    DATA = SetPCs(DATA,1);
elseif mode == 12 %restrict time range
    if ~isempty(DATA.restricttimerange)
        t(2)= DATA.restricttimerange(2);
    else
        t(2) = DATA.t(end)+0.01;
    end
    t(1)= DATA.t(DATA.spklst(1))-0.01;
    DATA = RestrictTimeRange(DATA,t);
    DATA = SetPCs(DATA,1);
elseif mode == 13 %Spool All Probes
    DATA.plotspk.allprobes = 1;
    SpoolAllSpikes(DATA);
elseif mode == 14 %Ues all spikes again
    DATA = UseAllEvents(DATA);
    set(DATA.toplevel,'UserData',DATA);
elseif mode == 15 %Xcorr for selected probes
    SetFigure(DATA.tag.xcorr, DATA,'front');
    ps = find(DATA.selectprobe);
    np = length(ps);
    for j = 1:np
        for k = 1:j
            P = DataClusters{ps(j)};
            Q = DataClusters{ps(k)};
            xc = xcorrtimes(P.times,Q.times);
            if k == j
                xc(201) = 0;
            end
            subplot(np,np,k+(j-1)*length(ps));
            bar(-200:200,xc);
            axis('tight');
            if k == j
                title(sprintf('P%d',ps(j)));
            end
        end
    end
elseif strcmp(mode,'allquickspks') %QuickSpks All probes
    DATA.plotspk.allprobes = 1;
    SetFigure(DATA.tag.allspikes, DATA);
    PlotQuickSpikes(DATA,100);
    SetFigure(DATA.tag.allxy, DATA);
    PlotAllProbes(DATA,'xy');
elseif strcmp(mode,'xcorradj')  %Xcorr for all adjacent probes
    xpts = linspace(DATA.spts(1),DATA.spts(end),401);
    for j = 1:DATA.nprobes-1
        xoff = floor((j-1)/6) .* length(DATA.spts);
        yoff = rem(j-1,6) .* DATA.vsep;
        P = DataClusters{j};
        Q = DataClusters{j+1};
        xc = xcorrtimes(P.times,Q.times);
        plot(xpts+xoff,yoff+(xc.*DATA.vsep./max(xc)),'k','linewidth',2);
        drawnow;
    end
elseif strcmp(mode,'xcorrall') %Xcorr for all adjacent probes
    GetFigure(DATA.tag.xcorr);
    ClearPlot;
    xpts = linspace(DATA.spts(1),DATA.spts(end),401);
    for j = 1:DATA.nprobes-1
    for k = 1:j-1
        mysubplot(DATA.nprobes,DATA.nprobes,(j-1) * DATA.nprobes+k);
        xoff = floor((j-1)/6) .* length(DATA.spts);
        yoff = rem(j-1,6) .* DATA.vsep;
        P = DataClusters{j};
        Q = DataClusters{k};
        [xc, details] = xcorrtimes(P.times,Q.times);
        h = plot(xpts+xoff,yoff+(xc.*DATA.vsep./max(xc)),'k','linewidth',2);
        [a,b]= max(xc);
        if a > 40 && a > mean(xc)+4*std(xc) && details.xpts(b) ~= 0 && abs(details.xpts(b)) < 5
            set(h,'color','r');
        end
        set(gca,'xtick',[],'ytick',[],'buttondownfcn',{@PlotXcorr, j,k});
        set(h,'buttondownfcn',{@PlotXcorr, j,k});
        axis('tight');
        drawnow;
    end
    end
end
if ismember(mode,[1 2 4 5 6 7 8 9]) | redraw
    DATA.spklst = DATA.spklst(1):DATA.spklst(1)+DATA.spksperview-1;
    PlotSpikes(DATA,DATA.spklst);
    set(DATA.toplevel,'UserData',DATA);
end

 function PlotClusterXY(DATA, C)

     cid = unique(C.clst);
     cid = cid(cid > 0);
     for j = 1:length(cid)
         c = cid(j);
         id = find(C.clst == c);
         plot(C.xy(id,1),C.xy(id,2),'.','markersize',1,'color',DATA.colors{c});
         hold on;
     end
     axis('tight');
     set(gca,'xtick',[],'ytick',[]);
     
function h = AddCellLabel(DATA,e,p)
    [iscell, cellid] = isacell(DATA, e, p);
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
    if iscell
        cellid = cellid(cellid > 0);
        yo=0;
        h =text(xl(2),yl(1)+yo,sprintf('%d:Cell %d',p,cellid(1)),'color',DATA.colors{2},'horizontalalignment','right','verticalalignment','bottom','fontsize',DATA.gui.fontsize(1));
    else
        [a,b] = find(squeeze(DATA.CellList(:,p,:)) > 0);
        if length(a)
            [c, id] = min(abs(a - e));
            cellid = DATA.CellList(a(id),p,b(id));
        h = text(xl(2),yl(1),sprintf('%d:(cell %d E%d)',p,cellid,a(id)),'color','r','horizontalalignment','right','verticalalignment','bottom','fontsize',DATA.gui.fontsize(1));
        else
        h = text(xl(2),yl(1),sprintf(' %d',p),'color','k','horizontalalignment','right','verticalalignment','bottom','fontsize',DATA.gui.fontsize(1));
        end
    end
    
    
function PlotAllProbes(DATA,type, varargin)
    ClusterDetails = getappdata(DATA.toplevel,'ClusterDetails');
    pid = 1:length(ClusterDetails);
    linewidth = 0.1;
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'linewidth',5)
            j = j+1;
            linewidth = varargin{j};
        elseif strncmpi(varargin{j},'probes',5)
            j = j+1;
            pid = varargin{j};
        end
        j = j+1;
    end
    [nr,nc] = Nsubplots(length(ClusterDetails));
    nr=4;
    nc=6;
    SetFigure(DATA.tag.allxy, DATA);
    for j = pid
        k = ceil(j/6) + (5-rem(j-1,6)).*4;
        mysubplot(6,4,k);
        hold off;
        PlotClusterXY(DATA, ClusterDetails{j});
        h = AddCellLabel(DATA, DATA.exptno,j);
        set(h,'ButtonDownFcn',{@SummaryHit, j});
        set(gca,'linewidth',linewidth);
    end
    
function SummaryHit(a,b, p)
    DATA = GetDataFromFig(a);
    AllVPcs(DATA.toplevel,'tchannew',p,'usecluster');
    figure(DATA.toplevel);

    
function PlotAllMeans(DATA)
    type = 'lineonly';

    Clusters = getappdata(DATA.toplevel,'Clusters');
    [nr,nc] = Nsubplots(length(Clusters));
    SetFigure(DATA.tag.allxy, DATA);
    for j = 1:length(Clusters)
        mysubplot(nr,nc,j);
        hold off; 
        PlotMeanSpikes(Clusters{j},0,1,type);
        set(gca,'Xtick',[],'Ytick',[]);
        h = get(gca,'title');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        a = get(h,'position');
        a(2) = yl(2);
        a(1) = mean(xl);
        set(h,'position',a,'VerticalAlignment','top');
        set(gca,'ButtonDownFcn',{@HitXYPlot, j});
        if DATA.selectprobe(j)
            set(gca,'xcolor','r','ycolor','r','linewidth',3);
        else
            set(gca,'xcolor','k','ycolor','k','linewidth',1);
        end
    end

function PlotMeanSpikes(C, p, cluster, varargin)
    addstr = [];
    plots = [ 1 1];
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'addtitile',5)
            j = j+1;
            addstr = varargin{j};
        elseif strncmpi(varargin{j},'lineonly',5)
            plots = [0 1];
        elseif strncmpi(varargin{j},'imageonly',7)
            plots = [1 0];
        end
        j = j+1;
    end
        
    if sum(plots) > 1
    subplot(1,2,1);
    end
    if plots(1)
    hold off;
    imagesc(C.MeanSpike.ms);
    if p <= 0 && isfield(C,'probe');
        p = C.probe(1);
    end
    line([0 5],[p p],'color','r');
    title(sprintf('P%d Ex %d Gm %.2f (%.2f) %s',p,C.exptno,C.mahal(1),C.mahal(2),addstr));
    if sum(plots) > 1
    subplot(1,2,2);
    end
    end
    if plots(2)
    hold off;
    v = std(C.MeanSpike.ms');
    id = find(v > max(v)/2);
    
    for j = id
        plot(C.MeanSpike.ms(j,:),'r');
        hold on;
        if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j
            plot(C.MeanSpike.dp(j,:),'g');
        end
        plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);
    end
    end

    function bad = BadCluster(C)
    bad = 0;
    if C.mahal(1) < 2
        bad = 1;
    end
    
    
 function DATA = UseAllEvents(DATA)
     if min(DATA.clst) < 0  %cluster was removed - need to recalc PCS
        recalc = 1;
    else
        recalc = 0;
    end
    if length(DATA.clst) < DATA.nevents
        clst = DATA.clst;
        DATA.clst = ones(DATA.nevents,1);
        DATA.clst(DATA.uid(clst)) = clst;
    end
    DATA.clst = abs(DATA.clst);
    DATA.clst(DATA.clst == 0) = 1;
    DATA.uid = 1:DATA.nevents;
    DATA.excludetrialids = [];
    DATA.restricttimerange = [];
    DATA.cluster = rmfields(DATA.cluster,'restricttimerange');
    DATA.excludetrialids = [];
    DATA.hidecluster = 0;
    DATA.currentcluster = 1;
    if recalc
        DATA = SetPCs(DATA,1);
    end
    SetTrialList(DATA);
    
function HitXYPlot(a,b, p)
DATA = GetDataFromFig(a);        
DATA.selectprobe(p) = ~DATA.selectprobe(p);
if DATA.selectprobe(p)
set(gca,'xcolor','r','ycolor','r','linewidth',3);
else
set(gca,'xcolor','k','ycolor','k','linewidth',1);
end
set(DATA.toplevel,'UserData',DATA);
    
function SpikeButtonPressed(a,b)
DATA = GetDataFromFig(a);        
DataClusters = getappdata(DATA.toplevel,'Clusters');
xy = get(a,'currentpoint');
npts = length(DATA.spts);
col = round((xy(1,1)+DATA.spts(1))./npts);
row = round(xy(1,2)./DATA.vsep)+1;
if col >= 0 && col < DATA.nprobes && row > 0 && row <= 6;
    p = row + col * 6;
    x = DATA.spts + col * npts;
    y = [DATA.vsep/2 -DATA.vsep/2] + (row-1) * DATA.vsep;
    imagesc(x,y,DataClusters{p}.MeanSpike.ms,'ButtonDownFcn',{@HitImage,p});
    plot(x([1 end]),[y(1)-(p-1) * DATA.vsep/(DATA.nprobes-1) y(1)-(p-1) * DATA.vsep/(DATA.nprobes-1)],'w--');
    DATA.selectprobe(p) = 1;
end
set(DATA.toplevel,'UserData',DATA);

function HitImage(a,b,p)
DATA = GetDataFromFig(a);        
DATA.selectprobe(p) = 0;        
delete(a);
%PlotTrialSpikes(DATA,DATA.currenttrial,'probes',p);

function SpoolAllSpikes(DATA, varargin)
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'front',3)
            SetFigure(DATA.tag.spikes, DATA,'front');
        end
            
        j = j+1;
    end
SetTrialList(DATA);
DataClusters = getappdata(DATA.toplevel,'Clusters');

DATA.vsep = 4;
np = DATA.nprobes;
if isempty(DATA.trigtimes) || sum(DATA.trigcheck) <np
    [DATA, DataClusters] = LoadTrigTimes(DATA, 1:np);
end
DATA.selectprobe = zeros(1,1:DATA.nprobes);
hold off;
plot(0,0,'w.');
hold on;
set(gca,'Position',[0.05 0.05 0.9 0.9]);
set(gca,'ButtonDownFcn',@SpikeButtonPressed);
npts = length(DATA.spts);
nc = 2;
DATA.handles = zeros(DATA.nprobes,2);
for p = 1:DATA.nprobes
    %nc = length(unique(DataClusters{p}.clst));
    DATA.handles(p,1) = plot([-8 32],[0 0], 'color',[0.5 0.5 0.5]);
    
    DATA.handles(p,2) = plot([-8 32],[0 0], 'color','r');
    if ~isfield(DataClusters{p},'auto')
        DataClusters{p}.auto = 1;
    end
    if DataClusters{p}.auto == 1 & BadCluster(DataClusters{p})
        set(DATA.handles(p,2),'color','b');
    elseif DataClusters{p}.auto ==1
        set(DATA.handles(p,2),'color',[0.8 0 0]);
    end
    if isfield(DataClusters{p},'xtimes')
        nx = length(DataClusters{p}.xtimes)
        for j = 1:nx
            DATA.handles(p,j+2) = plot([-8 32],[0 0], 'color','g');
        end
    end
    xoff = floor((p-1)/6) .* npts;
    yoff = (rem(p-1,6)) * DATA.vsep;
    text(xoff,yoff+1,sprintf('%d:%.1f',p,DataClusters{p}.mahal(1)));
end
set(gca,'xlim',[-8 (npts) * 4-8],'ylim',[-DATA.vsep DATA.vsep*6]); 
th =  title('Trial');
ts = now;
stopobj = findobj(gcf,'tag','StopSpool');
for nt = 1:length(DATA.Expt.Trials)
    PlotTrialSpikes(DATA,nt);
    spstop = get(stopobj,'value');
    if spstop
        break;
    end
%    set(th,'string',sprintf('Trial %d',nt));
end
set(stopobj,'value',0);



fprintf('Took %.2f\n',mytoc(ts));
    DATA.currenttrial = nt;
    set(DATA.toplevel,'UserData',DATA);

function PlotProbeSpikes(DATA, Vall, p, spklist,npts,offset)
    j = 1;
    xoff = offset(1);
    yoff = offset(2);
  hold on;
  if isempty(spklist)
      for j = 1:size(DATA.handles,2)
          if ishandle(DATA.handles(p,j)) & DATA.handles(p,j) > 0
              set(DATA.handles(p,j),'Ydata',[0 0]+yoff,'Xdata',[-npts(2) npts(1)-npts(2)]+xoff);
          end
      end
      return;
  end
DataClusters = getappdata(DATA.toplevel,'Clusters');

  x = [1:npts(1)]-npts(2);
  X = repmat(x,length(spklist),1)';
  id = repmat(spklist,npts(1),1) + X;
  nV = Vall.V(p,id);
  nV(npts(1)+1:npts(1):end) = NaN;
  X(end,:) = NaN;
  if ishandle(DATA.handles(p,1))
  set(DATA.handles(p,1),'Ydata',nV(:)+yoff,'Xdata',X(:)+xoff);
  else
      DATA.handles(p,1) = plot(X(:)+xoff,nV(:)+yoff,'color',DATA.colors{1});
  end
  tid = find(ismember(Vall.t(spklist),DataClusters{p}.times));
  if length(tid)
      X = repmat(x,length(tid),1)';
      id = repmat(spklist(tid),npts(1),1) + X;
      nV = Vall.V(p,id);
      nV(npts(1)+1:npts(1):end) = NaN;
      X(end,:) = NaN;
  end
  if length(tid) &&  ishandle(DATA.handles(p,2))
      set(DATA.handles(p,2),'Ydata',nV(:)+yoff,'Xdata',X(:)+xoff);
  elseif ishandle(DATA.handles(p,2))
      set(DATA.handles(p,2),'Ydata',[0 0]+yoff,'Xdata',[1 1]+xoff);
  else
      DATA.handles(p,2) = plot(X(:)+xoff,nV(:)+yoff,'color',DATA.colors{2});
  end
  for j = 1:length(DataClusters{p}.next)
      if isfield(DataClusters{p}.next{j},'times')
      tid = find(ismember(Vall.t(spklist),DataClusters{p}.next{j}.times));
      else
          tid = [];
      end
      if length(tid)
          X = repmat(x,length(tid),1)';
          id = repmat(spklist(tid),npts(1),1) + X;
          nV = Vall.V(p,id);
          nV(npts(1)+1:npts(1):end) = NaN;
          X(end,:) = NaN;
          if size(DATA.handles,2) < j+2 || ~ishandle(DATA.handles(p,j+2))
              DATA.handles(p,j+2) = plot(X(:)+xoff,nV(:)+yoff,'color',DATA.colors{j+2});
          else
          set(DATA.handles(p,j+2),'Ydata',nV(:)+yoff,'Xdata',X(:)+xoff);
          end
      end
  end

function SpoolSpikes(DATA, varargin)
    
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'front',3)
            SetFigure(DATA.tag.spikes, DATA,'front');
        end
            
        j = j+1;
    end
    
    if ~isfield(DATA,'plotcsd');
        DATA.plotcsd = 0;
    end
DATA.stopspool = 0;
SetTrialList(DATA);
AllV = GetAllV(DATA);
spklist = 1:DATA.spksperview;
%only do this for probes that will be used - it can take a long time
if DATA.plotspk.showfullv;
    usedprobes = union(DATA.plotspk.probes, UseProbeList(DATA,DATA.plotspk.nfullvprobes));
else
    usedprobes = DATA.plotspk.probes;
end
for j = usedprobes
if DATA.plotcsd
    AllCSD = GetCSD(DATA,DATA.plotcsd);
maxV(j,:,:) = max(max(AllCSD(j,:,:),[],3),[],2);
minV(j,:,:) = min(min(AllCSD(j,:,:),[],3),[],2);
else
maxV(j,:,:) = max(max(AllV(j,:,:),[],3),[],2);
minV(j,:,:) = min(min(AllV(j,:,:),[],3),[],2);
end
end
voffset = cumsum(max(abs([maxV minV]),[],2));
DATA.voffset(1:length(voffset)) = voffset;
DATA.spkyrange(1) = min(DATA.voffset(DATA.plotspk.probes)-DATA.voffset(DATA.probe(1))) + minV(min(DATA.plotspk.probes));
DATA.spkyrange(2) = max(DATA.voffset(DATA.plotspk.probes)-DATA.voffset(DATA.probe(1))) + maxV(max(DATA.plotspk.probes));
set(DATA.toplevel,'UserData',DATA);

go = 1;
nt = 0;
if DATA.plotspk.bytrial
    id = find([DATA.Expt.Trials.TrialStart] > DATA.t(DATA.uid(1)).*10000);
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
stopobj = findobj(gcf,'tag','StopSpool');
h = NaN;
while (isempty(spklist) || spklist(end) < DATA.nevents) && go
    if DATA.plotspk.bytrial
        nt = nt+1;
        spklist = find(DATA.t .*10000 > DATA.Expt.Trials(nt).Start(1) & DATA.t.*10000 < DATA.Expt.Trials(nt).End(end));
        if nt >= length(DATA.Expt.Trials) || DATA.Expt.Trials(nt).Start(1) > maxtime
            go = 0;
        elseif DATA.Expt.Trials(nt).Start(1) < mintime
            go = 2;
        else
            go = 1;
        end
    else
        spklist = spklist + DATA.spksperview;
    end
    if go == 1
    DATA.currenttrial = nt;
    h = PlotSpikes(DATA,spklist,'fix',DATA.spkyrange);
    drawnow;
    DATA = get(DATA.toplevel,'UserDATA');
    end
    spstop = get(stopobj,'value');
%    set(th
    if DATA.stopspool || spstop
        go = 0;
        DATA.spklst = spklist;
    end
end
set(stopobj,'value',0);
D = get(gcf,'UserData');
D.plottype = 'SpoolSpikes';
set(gcf,'UserData',D);

DATA.currenttrial = nt;
DATA.stopspool = 0;
if ishandle(h)
AddLineContextMenu(DATA, h);
DATA = AddAxisContextMenu(DATA, gca);
end
set(DATA.toplevel,'UserData',DATA);

function csd = GetCSD(DATA, ndiff)
    csd = getappdata(DATA.toplevel,'AllCSD');
    currdiff = getappdata(DATA.toplevel,'plotcsd');
    if isempty(csd) || ndiff ~= currdiff
        AllV = GetAllV(DATA);
        csd = diff(AllV,3-ndiff,1);
        setappdata(DATA.toplevel,'AllCSD',csd);
    end
    setappdata(DATA.toplevel,'plotcsd',ndiff);

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
    DATA.currenttrial = PlotTrialSpikes(DATA,nt);
    set(DATA.toplevel,'UserData',DATA);
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

if spk(2) > 0 && spk(1) < DATA.nevents
DATA.spklst = spk(1):spk(2);
PlotSpikes(DATA,DATA.spklst);
end
set(DATA.toplevel,'UserData',DATA);
    
        

function PlotPCs(pcs, a,b, type, id, colors, C,varargin)
ptsz = [1 1];
fixrange  = 0;
nid = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fixrange',5)
        fixrange = 1;
    elseif strncmpi(varargin{j},'ptsz',4)
        j = j+1;
        ptsz = varargin{j};
    elseif strncmpi(varargin{j},'setnid',5)
        j = j+1;
        nid = varargin{j};
        ptsz =  [1 8];
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
            plot(pcs(nid,a),pcs(nid,b),'.','markersize',ptsz(2),'color',colors{nc(j)});
            hold on;
        end
        else
            plot(pcs(id,a),pcs(id,b),'r.','markersize',ptsz(1));
        end
    elseif length(id)
        if isempty(nid)
        nid = setdiff(1:size(pcs,1),id);
        end
        if length(nid)
        plot(pcs(nid,a),pcs(nid,b),'.','markersize',ptsz(2));
        hold on;
        end
        plot(pcs(id,a),pcs(id,b),'r.','markersize',ptsz(1));
        hold on;
    else
        plot(pcs(:,a),pcs(:,b),'.','markersize',ptsz(2));
    end
    if fixrange
        set(gca,'xlim',xl);
        set(gca,'ylim',yl);
    end
else
    hold off;
    if length(id) == size(pcs,1)
        id = find(id > 0) ; %dont include spikes from excluted trials
        [a, D] = DensityPlot(pcs(id,a),pcs(id,b),'sd',[2 2],'ynormal');
    else
    [a,D] = DensityPlot(pcs(:,a),pcs(:,b),'sd',[2 2],'ynormal');
    end
    if type == 2
        r = CalcRadius(C,[D.x(:) D.y(:)]);
        id = find(r < 1);
        cmax = max(D.z(id));
        caxis([0 cmax]);
    end
end

function r = CalcRadius(E,xy)
    
    rx = E.xyr(3);
    ry = E.xyr(4);
    xys = xyrotate(xy(:,1)-E.xyr(1),xy(:,2)-E.xyr(2),E.angle);
    r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;

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
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    if length(id) == DATA.nevents
        nc = unique(id);
        nc = nc(nc > 0);
        for j = 1:length(nc)
            ids{j} = find(id == nc(j));
        end
    elseif length(id)
        ids{2} = id;
        ids{1} =setdiff(1:DATA.nevents,id);
    else 
        ids{1} = 1:DATA.nevents;
    end
    if DATA.plotcsd
        AllCSD = GetCSD(DATA,DATA.plotcsd);
        for j = 1:length(ids)
            plot(squeeze(AllCSD(a(1),a(2),ids{j})),squeeze(AllCSD(b(1),b(2),ids{j})),'.','markersize',ptsz(1),'color',colors{nc(j)});
            hold on;
        end
    elseif DATA.plottype == 9
        for j = 1:length(ids)
            plot(squeeze(DATA.energy(a(1),ids{j})),squeeze(DATA.energy(b(1),ids{j})),'.','markersize',ptsz(1),'color',colors{j});
            hold on;
        end        
    elseif DATA.plottype == 6
        for j = 1:length(ids)
            plot(squeeze(DATA.dV(a(1),a(2),ids{j})),squeeze(DATA.dV(b(1),b(2),ids{j})),'.','markersize',ptsz(1),'color',colors{j});
            hold on;
        end        
    else
        AllV = GetAllV(DATA);
        for j = 1:length(ids)
        plot(squeeze(AllV(a(1),a(2),ids{j})),squeeze(AllV(b(1),b(2),ids{j})),'.','markersize',ptsz(1),'color',colors{j});
        hold on;
        end
    end
    axis('tight');
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');

    if fixrange
        set(gca,'xlim',xl);
        set(gca,'ylim',yl);
    end
    title(sprintf('%d:%d vs %d:%d',a(1),a(2),b(1),b(2)),'position',[mean(xl) yl(2)],'VerticalAlignment','Top');
else
    AllV = GetAllV(DATA);
    DensityPlot(AllV(a(1),a(2),:),AllV(b(1),b(2),:),'sd',[2 2],'ynormal');
end

function DATA = ClassifyAll(DATA, force,varargin)
    
    quickcutmode.quick = 0;
    cargs = {};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'quick',5)
            quickcutmode = DATA.quickcutmode;
            quickcutmode.quick = 1;
            cargs = {cargs{:} 'quick'};
        elseif strncmpi(varargin{j},'recluster',5)
            cargs = {cargs{:} varargin{j}};
        end
        j = j+1;
    end
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    if DATA.interactive >= 0
    oldname = get(DATA.toplevel,'Name');
    set(DATA.toplevel,'Name','Classifying...');
    drawnow;
    end

    ctimes(1) = DATA.cluster.ctime;
    iscluster(1) = 1;
    e = DATA.cluster.exptno;
    for j = 1:length(DATA.cluster.next)
        if isfield(DATA.cluster.next{j},'ctime')
            ctimes(j+1) = DATA.cluster.next{j}.ctime;
        else
            if isempty(DATA.cluster.next{j})
                ctimes(j+1) = 0;
            else
                errordlg(sprintf('Missing ctime Cluster %d E%.1f',j+1,e),'Cluster Cut Time','modal');
                ctimes(j+1) = j+1;
            end
        end
        if isfield(DATA.cluster.next{j},'space')
            iscluster(j+1) = 1;
        else
            iscluster(j+1) = 0;
        end
    end
    cls = 1:(1+length(DATA.cluster.next));
    cls = cls(ctimes>0 & iscluster > 0);
    ctimes = ctimes(ctimes > 0 & iscluster > 0);
    [a, tlist] = sort(ctimes);
    tlist = cls(tlist);
    cc = DATA.currentcluster;
    if ~isfield(DATA.cluster,'chspk')
        DATA.cluster.chspk = DATA.chspk;
    end

    autorefine = 0;
    for j = 1:length(tlist)
        needcl = NeedClusterData(DATA.cluster,tlist(j));
        DATA.currentcluster = tlist(j);
        if force
            DATA = SetTemplateData(DATA,tlist(j));
        end
        DATA.cluster = CheckScoreScaling(DATA, DATA.cluster);
        if tlist(j) == 1
            DATA.currentcluster = 1;
%needcl == 2 means that already classified, but fits not done. If force is
%on, the cut is going to be redone anyway
            if needcl == 2  && force == 0
                if DATA.interactive > 0
                    set(DATA.toplevel,'Name',sprintf('Quantifying Cluster %d...',tlist(j)));
                    drawnow;
                end
                DATA.cluster = ClassifyFit(DATA, DATA.cluster,1);
            elseif force || needcl
                if DATA.interactive > 0
                    set(DATA.toplevel,'Name','Classifying Cluster 1...');
                    drawnow;
                end
%if autorefine, then will have to claissfy twice if ther are > 1clusters.
%So on first pass, only do the quck classif. But set quick to zero so that
%the second time it treated as new.
                if DATA.autorefine > 0 && isfield(DATA.cluster,'fitdprime') && DATA.cluster.fitdprime(1) < DATA.autorefine
                    [cl, cluster, DATA.xy{1}] = ClassifySpikes(DATA,DATA.cluster,'quick','noplot');
                    cluster.quick = 0;
                    DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(BoundaryFromCluster([],cluster, DATA.currentcluster));
                    DATA.cluster.clst = cl.clst;
                    DATA = OptimizeBoundary(DATA);
                    autorefine = autorefine+1;
                else
                    [cl, cluster, DATA.xy{1}] = ClassifySpikes(DATA,DATA.cluster,cargs{:},'noplot','notrig');
                    DataClusters{DATA.probe(1)}.mahal = cluster.mahal;
                    DATA.clusterboundary{DATA.currentcluster} = CondenseCluster(BoundaryFromCluster([],cluster, DATA.currentcluster));
                end
                DATA.clid = cl.id;
                DATA.nid = cl.nid;
                DATA.clst = cl.clst;

                if length(cl.id) > 1 && DATA.autorefine == 0
                    if quickcutmode.quick == 1 
                        cl.MeanSpike = DATA.cluster.MeanSpike;
                    end
                    DATA.MeanSpike = cl.MeanSpike;
                    DATA.cluster = cluster;
                    DATA.cluster.MeanSpike = cl.MeanSpike;
                    DATA.cluster.minspke = prctile(DATA.energy(1,DATA.clid),1) .* 0.95;
                    DATA.cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),DATA.clid),1) .* 0.95;
                    DATA.cluster.ctime = now;
                elseif length(cl.id) <= 1
                    DATA.cluster = cluster;
                end
                if DATA.cluster.mahal(4) < 0.001 && quickcutmode.quick == 0
                    c = PlotHistogram(DATA, []);
                    DATA.cluster.mahal(4) = c.gmdprime;
                end
            end
        elseif isfield(DATA.cluster.next{tlist(j)-1},'space')
            if needcl == 2 && force == 0
                if DATA.interactive > 0
                set(DATA.toplevel,'Name',sprintf('Quantifying Cluster %d...',tlist(j)));
                drawnow;
                end
                DATA.cluster.next{tlist(j)-1} = ClassifyFit(DATA, DATA.cluster.next{tlist(j)-1},tlist(j));
            elseif force || needcl
                if DATA.interactive > 0
                set(DATA.toplevel,'Name',sprintf('Classifying Cluster %d...',tlist(j)));
                drawnow;
                end

            DATA.currentcluster = tlist(j);
            if DATA.autorefine > 0 && abs(DATA.cluster.next{tlist(j)-1}.fitdprime(1)) > DATA.autorefine
                [xcl, DATA.cluster, DATA.xy{tlist(j)}] = ClassifySpikes(DATA,DATA.cluster,'quick','noplot','notrig');
                DATA = OptimizeBoundary(DATA);
                autorefine = autorefine+1;
            else
                [xcl, DATA.cluster, DATA.xy{tlist(j)}] = ClassifySpikes(DATA,DATA.cluster,'noplot', 'notrig',cargs{:});
            end
            if ~isempty(xcl)
                DATA.clst = xcl.clst;
            end
            if DATA.cluster.next{tlist(j)-1}.mahal(4) < 0.001 && quickcutmode.quick == 0
                c = PlotHistogram(DATA, [], quickcutmode);
                DATA.cluster.next{tlist(j)-1}.mahal(4) = c.gmdprime;
            end
            %
            DATA.cluster.next{tlist(j)-1}.ctime = now;
            if isfield(xcl,'MeanSpike') %%may not be true on first "quick" pass
            DATA.cluster.next{tlist(j)-1}.MeanSpike = xcl.MeanSpike;
            end
            DATA.clusterboundary{tlist(j)} = CondenseCluster(BoundaryFromCluster([],DATA.cluster, tlist(j)));
            end
        end
%if only calculating metrics, don't call this       
        if DATA.iteratefit && DATA.autorefine > 0 && (needcl == 1 || force)
              IterateFit(DATA, DATA.iteratefit);
              DATA = get(DATA.toplevel,'UserData');
        end
        if DATA.watchplots && (needcl || force) && quickcutmode.quick == 0
            PlotHistogram(DATA, []);
            ReplotPCs(DATA,[]);
        end
    end
%Unrefined Cluster N can steal events from refined cluster N-1, so need to appl one more time 
%and need to call Classify again anyway to apply refined bounday.
%Also set cluster.times in each Cluster now for the same reason

    setappdata(DATA.toplevel,'Clusters',DataClusters);
    if autorefine 
        DATA.autorefine = 0;
        DATA = ClassifyAll(DATA,1);
        DATA.autorefine = 1;
    end
    
%    quantify dropi for other clusters
    if quickcutmode.quick == 0 || quickcutmode.dropi
        for j = 1:length(DATA.cluster.next)
            if ~isempty(DATA.cluster.next{j})
                cluster = PlotTriggerHist(DATA,DATA.cluster.next{j},quickcutmode);
                DATA.cluster.next{j}.vhist = cluster.vhist;
                DATA.cluster.next{j}.dropi = cluster.dropi;
                DATA.cluster.next{j}.times = DATA.t(DATA.clst == j+2);
            end
        end
    end
    DATA.cluster = PlotTriggerHist(DATA,DATA.cluster,'showall',quickcutmode);
    if DATA.watchplots > 1
        drawnow;
    end
    DATA.cluster.times = DATA.t(DATA.clst == 2);
    DATA.currentcluster = cc;
    DATA.clid = find(DATA.clst == cc+1);
    if DATA.checkclusters
        CheckClusters(DataClusters, 'CheckNexts','Classify');
        CheckClusters(DataClusters,'CheckFitSpace');
    end
    if DATA.interactive >= 0
    set(DATA.toplevel,'UserData',DATA);
    set(DATA.toplevel,'Name',oldname);
    end
    DATA.endtime = now;

  function need = NeedClusterData(Cluster, ci)
      need = 0;
      if ci > 1
          C = Cluster.next{ci-1};
      else
          C = Cluster;
      end
      % if there is no 'quick' field, should mean this has not been applied
      % yet. quick is set to 0 when its been applied. Useful to
      % differentiae cuts from command line(1) from interactive changes (2)
      if ~isfield(C, 'quick')
          need = 1;
      elseif C.quick > 0
          need = 2;
      end
      
function D = CondenseCluster(C)
    D = rmfields(C,{'gmfit2dman' 'gmfit2d' 'gmfit1d' 'gmfit'});
    
function DATA = IterateFit(DATA, niter)
j = 1;
xc = [0 0];

if DATA.currentcluster == 1
    C = DATA.cluster;
else
    C = DATA.cluster.next{DATA.currentcluster-1};
end
  
ReplotPCs(DATA,[]);
for j = 1:niter
    DATA.MeanSpike = PlotMeanSpike(DATA,'recalc');
    TemplatePlot(DATA,'nodip','noplot');
    DATA = get(DATA.toplevel,'UserData');
    T{j} = DATA.MeanSpike.ms;
    if C.space(1) == 6
        G = DATA.cluster.gmfit;
%        xy = ProjectND(DATA, size(G.mu,2), G);
        OptimizeBoundary(DATA);
    elseif C.space(1) == 3
        DATA = OptimizeBoundary(DATA);
        ReplotPCs(DATA,[]);
    end
end
set(DATA.toplevel,'UserData',DATA);

    
 function DATA = SetTemplateData(DATA,  cl, varargin)
     
     force = 0;
     j = 1;
     while j <= length(varargin)
         if strncmpi(varargin{j},'force',4)
             force = 1;
         end
         j = j+1;
     end
     DataClusters = getappdata(DATA.toplevel,'Clusters');
     if cl ==1 
         C = DATA.cluster;
     else
         C = DATA.cluster.next{cl-1};
         if ~isfield(C,'TemplateUsed') && isfield(DATA.cluster,'TemplateUsed')
             C.TemplateUsed = DATA.cluster.TemplateUsed;
             if isfield(DATA.cluster,'DprimeUsed')
                 C.DprimeUsed = DATA.cluster.DprimeUsed;
             end
             DATA.cluster.next{cl-1} = C;
         end
         if ~isfield(C,'chspk') 
             if isfield(DATA.cluster,'chspk')
                 C.chspk = DATA.cluster.chspk;
             else
                 C.chspk = DATA.chspk;
             end
         end
     end
     p = DATA.cluster.probe(1);
     [needtemplate, plottypes] = NeedTemplateForCluster(C,0);
     
     if needtemplate || force
         recalc = 0;
         if ~isfield(DATA,'TemplateScores')
             recalc = 1;
         elseif isfield(C,'TemplateUsed') && sum(size(DATA.TemplateUsed) == size(C.TemplateUsed)) <2
             recalc = 1;
         elseif sum(abs(DATA.TemplateUsed(:)-C.TemplateUsed(:))) > 0.1 .* std(DATA.TemplateUsed)
             recalc = 1;
         elseif max(DATA.uid) > size(DATA.TemplateScores,1) %retriggered
             recalc = 1;
         end
             
         if isfield(C,'TemplateUsed') && recalc
             
             if cl == 1
                 DATA.cluster.neednewtemplate = 0;
             else
                 DATA.cluster.next{cl-1}.neednewtemplate = 0;
             end
             if DATA.plottype == 7
                 DATA = CalcTemplatesFromMean(DATA,DATA.cluster.TemplateUsed,'stdtemplate');
             else
                 if size(C.TemplateUsed,1) == length(C.chspk)
                     C.MeanSpike.ms(DATA.cluster.chspk,:) = C.TemplateUsed;
                     if isfield(C,'DprimeUsed') && size(C.DprimeUsed,1) == length(DATA.chspk)
                         C.MeanSpike.vdprime(DATA.cluster.chspk,:) = C.DprimeUsed;
                     end
                     if isfield(C,'mumeanUsed') && size(C.mumeanUsed,1) == length(DATA.chspk)
                         C.MeanSpike.mu(DATA.cluster.chspk,:) = C.mumeanUsed;
                     end
                     if isfield(C.MeanSpike,'othermeans')
                         C.MeanSpike.othermeans = C.MeanSpike.othermeans;
                     end
                     if isfield(C,'next') && length(C.next) &&  isfield(C.next{1},'TemplateUsed')
                         if size(C.next{1}.TemplateUsed,1) == length(C.chspk) && isfield(C.next{1},'MeanSpike')
                             C.MeanSpike.othermeans{1} = C.next{1}.MeanSpike.ms;
                             C.MeanSpike.othermeans{1}(C.chspk,:) = C.next{1}.TemplateUsed;
                         else
                             C.MeanSpike.othermeans{1} = C.next{1}.TemplateUsed;
                         end
                     end
                     DATA = CalcTemplatesFromMean(DATA,C.MeanSpike);
                 else
                     MeanSpike = C.MeanSpike;
                     if isempty(C.TemplateUsed)
                         DATA = AddErr(DATA,'Empty TemplateUsed in  C%d',C.probe(1));
                     else
                         MeanSpike.ms = C.TemplateUsed;
                     end
%if calculations used dVdt, size will be 1 shorter                     
                     if isfield(C,'DprimeUsed') && size(C.DprimeUsed,1) == length(DATA.cluster.chspk)
                         MeanSpike.vdprime(DATA.cluster.chspk,1:size(C.DprimeUsed,2)) = C.DprimeUsed;
                     end
                     if isfield(C,'mumeanUsed') && size(C.mumeanUsed,1) == length(DATA.cluster.chspk)
                         MeanSpike.mu(DATA.cluster.chspk,:) = C.mumeanUsed;
                     end
                     if isfield(DATA.cluster.MeanSpike,'othermeans')
                         MeanSpike.othermeans = DATA.cluster.MeanSpike.othermeans;
                     end
                     if cl == 1
                         if isfield(C,'next') && length(C.next) && isfield(C.next{1},'TemplateUsed')
                             if size(C.TemplateUsed,1) == length(C.chspk)
                                 MeanSpike.othermeans{1} = C.next{1}.MeanSpike.ms;
                                 MeanSpike.othermeans{1}(C.chspk,:) = C.next{1}.TemplateUsed;
                             else
                                 MeanSpike.othermeans{1} = C.next{1}.TemplateUsed;
                             end
                         end
                     else
                         MeanSpike.othermeans{1} = DATA.cluster.MeanSpike.ms;
                     end
                     nm=1;
                     if length(DataClusters) >= p
                     if p > 1 && isfield(DataClusters{p-1},'TemplateUsed')
                         nm=nm+1;
                         MeanSpike.othermeans{nm} = DataClusters{p-1}.TemplateUsed;
                     end
                     if p < DATA.nprobes && isfield(DataClusters{p+1},'TemplateUsed')
                         nm=nm+1;
                         MeanSpike.othermeans{nm} = DataClusters{p+1}.TemplateUsed;
                     end
                     end
                     DATA = CalcTemplatesFromMean(DATA,MeanSpike);
                 end
                 if cl > 1
                     DATA.cluster.next{cl-1}.DprimeUsed = DATA.DprimeUsed;
                 else
                     DATA.cluster.DprimeUsed = DATA.DprimeUsed;
                 end
             end
         elseif recalc
             DATA = CalcTemplatesFromMean(DATA,C.MeanSpike);
         end
     end

function DATA = ClassifyAndFit(DATA)

    oldname = get(DATA.toplevel,'Name');
    set(DATA.toplevel,'Name','Classifying...');
    DataClusters = getappdata(DATA.toplevel,'Clusters');

    drawnow;
    if isfield(DATA.cluster,'sign') && abs(DATA.cluster.sign) == 1
        [cl, DATA.cluster, DATA.xy{DATA.currentcluster}] = ClassifySpikes(DATA, DATA.cluster,'sign',DATA.cluster.sign);
    else
        [cl, DATA.cluster, DATA.xy{DATA.currentcluster}] = ClassifySpikes(DATA, DATA.cluster);
    end
    if length(cl.id)
        DATA.MeanSpike = cl.MeanSpike;
    end
    DataClusters{DATA.probe(1)} = DATA.cluster;
    DataClusters{DATA.probe(1)}.times = DATA.t(cl.id);
    setappdata(DATA.toplevel,'Clusters',DataClusters);
    DATA.plottype = WhichPlotType(DATA.cluster, DATA.currentcluster);
    c = PlotHistogram(DATA, []);
    if DATA.cluster.mahal(4) < 0.001
        DATA.cluster.mahal(4) = c.gmdprime;
    end
    if DATA.checkclusters
    CheckClusters(DataClusters, 'CheckNexts','Classify');
    CheckClusters(DataClusters,'CheckFitSpace');
    end
    set(DATA.toplevel,'UserData',DATA);
    set(DATA.toplevel,'Name',oldname);
    drawnow;
    
function cluster = ClassifyFit(DATA, E, cnum)
    xy = DATA.xy{cnum};
    id = find(DATA.clst == cnum+1);
    nid = find(DATA.clst ~= cnum+1 & DATA.clst > 0);
cluster = E; %must always be a cluster, not a "boundary"
cluster.dprime = CalcDprime(xy(id,1),xy(nid,1));
cluster.hdip = HartigansDipTest(sort(xy(DATA.uid,1)));
cluster.clst = DATA.clst;
if ismember(cluster.shape,[1 2]) %cut with line
    [a,b] = GMDip(xy(DATA.uid,:),DATA.energy(1,DATA.uid),'crit',cluster.crit,'label',DATA.idstr);
    [dp,dpfits, dpdetails] = Fit2Gauss(cluster, xy(:,1),DATA);
    dip = MyDip(xy(DATA.uid,1));
else
    if ~isfield(cluster, 'r')
        cluster.r = CalcClusterDistance(cluster, xy);
    end
    [a,b] = GMDip(Rprime(cluster.r(DATA.uid,:)),DATA.energy(1,DATA.uid),'crit',1,'label',DATA.idstr);
    [dp,dpfits, dpdetails] = Fit2Gauss(cluster, cluster.r,DATA);
    dip = MyDip(cluster.r(DATA.uid));
end
cluster.mydip = dip.x(dip.dip);
cluster.mydipsize = dip.dipsize;
cluster.fitdprime = [dp dpdetails.fitpos dpdetails.dx dpdetails.minxpt];
cluster.fitdpparams = cat(1,dpfits{1}.params, dpfits{2}.params);
cluster.autodipsize = b.dipsize;
cluster.dipsize = b.cdipsize;
cluster.gmdprime = b.gmdprime;
cluster.gmfit1d = b.G{b.best};
x = FitGaussMeans(xy(DATA.uid,:),2,'clusterid',find(ismember(DATA.uid,id)));

cluster.bmc = BimodalCoeff(cluster.r,1.5);
a = FitGaussMeans(xy(DATA.uid,:),2);
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
if isfield(cluster,'gmdprime')
    cluster.mahal(4) = cluster.gmdprime;
elseif isfield(E,'gmdprime')
    cluster.mahal(4) = E.gmdprime;
elseif isfield(E,'gmfit1d') && isfield(E.gmfit1d,'mu')
    cluster.mahal(4) = GMdprime(E.gmfit1d);
else
    cluster.mahal(4) = 0;
end
if DATA.watchplots
    PlotHistogram(DATA, cluster);
end
cluster = PlotTriggerHist(DATA, cluster);

cl.MeanSpike = PlotMeanSpike(DATA,'recalc');
cluster.MeanSpike = cl.MeanSpike;

 cluster.quick = 0;        

 
function r = CalcClusterDistance(cluster, xy)
    cx = cluster.xyr(1);
    cy = cluster.xyr(2);
    rx = cluster.xyr(3);
    ry = cluster.xyr(4);

    if isfield(cluster,'aspectratio') & cluster.aspectratio > 0
        xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./cluster.aspectratio,cluster.angle);
        r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./cluster.aspectratio)).^2;
    else
        xys = xyrotate(xy(:,1)-cx,xy(:,2)-cy,cluster.angle);
        r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
    end
    
function X = GetDataStruct(DATA, f)
        
    
    if DATA.interactive >= 0
        X = getappdata(DATA.toplevel,f);
    else
        X = DATA.(f);
    end

function [cl, cluster, xy] = ClassifySpikes(DATA, E, varargin)

    quick = 0;
    forcesign = 0;
    replot = 1;
    recluster = 0;
    j = 1;
%    if quickmode.quick is nonzero, then avoid slow steps that are not
%    always neeeded.  other flags in quickmode detemine what fits are used.
    quickmode.quickest = 0;
    quickmode.quick = 0;
    quickmode.fit1cut = 0;
    quickmode.fit2gauss = 0;
    quickmode.dropi = 0;
    quickmode.mean = 0;
    quickmode.triggerhist = 1;
    
    while j <= length(varargin)
        if isfield(varargin{j}, 'quickest')
            quickmode = varargin{j};
            quickmode.quick = 1;
        elseif strncmpi(varargin{j},'noplot',4)
            replot = 0;
        elseif strncmpi(varargin{j},'notrigger',6)
            quickmode.triggerhist = 0;
        elseif strncmpi(varargin{j},'quick',4)
            quickmode.quick = 1;
        elseif strncmpi(varargin{j},'recluster',4) %Not chnaging
            recluster = 1;
        elseif strncmpi(varargin{j},'sign',4)
            j = j+1;
            forcesign = varargin{j};
        end
        j = j+1;
    end
    cl = [];
    cluster = [];
    xy = [];

    if ~isfield(E,'next') && isfield(DATA.cluster,'next')
        E.next = DATA.cluster.next;
    end
    
    if ~isfield(E,'pcplot') && isfield(E,'space')
        cluster = E;
        if cluster.space(1) == 6
            if cluster.space(2) == 3
                DATA.plottype = 2;
            elseif cluster.space(2) == 4
                DATA.plottype = 3;
            else
                DATA.plottype = cluster.space(2);
            end
        elseif cluster.space(1) == 5  %Var/E plot is another plot
            DATA.plottype = 1;
        else
            DATA.plottype = cluster.space(1);
        end
        if DATA.currentcluster > 1 & isempty(cluster.next{DATA.currentcluster-1})
            return;
        end
        E = BoundaryFromCluster([],cluster, DATA.currentcluster);

    else
        E.boundarytype = 1;
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
    cluster.nspks = DATA.nevents;
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
    cluster.probe = ProbeNumber(DATA);
    cluster.ctype = 0;
    if length(DATA.excludetrialids)
        cluster.excludetrialids = DATA.excludetrialids;
    end
    clnum = DATA.currentcluster+1;
    
    if isfield(DATA.Expt,'exptno')
        cluster.exptno = DATA.Expt.exptno;
    end
    f = {'bestspace' 'bestd' 'auto' 'autotook' 'firstspace' 'firstbmi' 'bestll' 'gmdip' 'gmdipres' 'gmfit1d' 'gmfit' 'gmfits' 'aspectratio'...
        'next' 'pos' 'TemplateUsed' 'mumeanUsed' 'MeanSpike' 'eventrate' 'exptreadmethod'...
        'mydip' 'mydipsize' 'manual'};
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
    else
        cl.clst = ones(DATA.nevents,1);
    end
    if isempty(cl.clst)
        cl.clst = ones(DATA.nevents,1);
    end
        
    if isfield(E,'cluster')
        cluster.cluster = E.cluster;
    elseif DATA.hidecluster
        cluster.cluster = DATA.hidecluster+1;
    else
        cluster.cluster = 1;
    end
    cluster.cluster = clnum-1;
    if clnum == 2
        C = cluster;
    elseif length(cluster.next) > clnum-3
        C = cluster.next{clnum-2};
    else %new cut 
        C = [];
    end
    
    if DATA.verbose > 1
        fprintf('Classify: %s Space %s\n',IDStr(DATA),sprintf('%d ',E.space));
    end
    if isfield(E,'usegmcluster') && E.usegmcluster == 1
        xy = DATA.ndxy;
        id = find(E.bestcl == 1);
        nid = find(E.bestcl == 2);
        e(1) = mean(DATA.energy(1,id));
        e(2) = mean(DATA.energy(1,nid));
        if diff(e) > 0
            id = find(E.bestcl == 2);
            nid = find(E.bestcl == 1);
            cluster.sign = -1;
        else
            cluster.sign = 1;
        end
        r = xy(:,1);
        cluster.space = E.space;
        cluster.crit = (mean(xy(id,1))+mean(xy(nid,1)))/2;
        x = FitGaussMeans(xy(DATA.uid,:),2,'clusterid',id);
        cluster.shape = 3;
    elseif E.shape == 2 || (E.space(1) == 6 && E.shape == 1) %line in n-D space
        cluster.angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),cluster.angle);
        cluster.y = exy(:,2);
        cluster.crit = mean(exy(:,1));
        cluster.sign = E.sign;
        cluster.space = E.space;
        xy = xyrotate(DATA.ndxy(:,1),DATA.ndxy(:,2),cluster.angle);
        r = xy(:,1);
        if cluster.sign < 0
            id = find(xy(DATA.uid,1) < cluster.crit);
            nid = find(xy(DATA.uid,1) >= cluster.crit);
        else
            cluster.sign = 1;
            id = find(xy(DATA.uid,1) > cluster.crit);
            nid = find(xy(DATA.uid,1) <= cluster.crit);
        end
        aid = id;
        id = DATA.uid(id);
        nid = DATA.uid(nid);
        e(1) = mean(DATA.energy(1,id));
        e(2) = mean(DATA.energy(1,nid));
        if e(2) > e(1);
            cluster.sign = -cluster.sign;
            cid = nid;
            nid = id;
            id = cid;
        end
        cluster.bmc = BimodalCoeff(r,1.5);
        cluster.dprime = CalcDprime(r(id),r(nid));
        cluster.hdip = HartigansDipTest(sort(r));
        x = FitGaussMeans(xy(DATA.uid,:),2,'clusterid',aid);
%        [a,b] = FindDip(r,DATA.energy(1,:),'eval',cluster.crit);
        if ~isfield(E,'gmdprime')
            [a,b] = GMDip(xy(DATA.uid,1),DATA.energy(1,DATA.uid),'eval',cluster.crit,'label',DATA.idstr);
            cluster.gmdprime = b.gmdprime;
            cluster.mahal(4) = b.gmdprime;
            cluster.autodipsize = b.dipsize;
            cluster.dipsize = b.cdipsize;
            cluster.gmfit1d = b.G{b.best};
        end

    elseif E.shape == 1 %line
% for line cuts, rotate space so that x > crit defines Cluster
% ClsuterDetails.xy will then have rotated values. Need rotating back in
% Plotclsuters.
        angle = atan(diff(E.pos([1 3]))/diff(E.pos([2 4])));
        exy = xyrotate(E.pos([1 3]),E.pos([2 4]),angle);
        cluster.y = exy(:,2);
        crit = mean(exy(:,1));
        cluster.angle = angle;
        cluster.crit = crit;
        if isempty(p) %Var/E plot
            y = DATA.spkvar(ispk,DATA.uid)./DATA.energy(1,DATA.uid);
            xy = xyrotate(DATA.energy(1,:),y,angle);
           cluster.space = [5];
        elseif DATA.plottype == 2 %voltage pairs
            AllV = GetAllV(DATA);
            xy = xyrotate(AllV(p(1),p(2),:),AllV(p(3),p(4),:),angle);
           cluster.space = [2 p(1:4)];
        elseif ismember(DATA.plottype,[3 4 7])
            xy = xyrotate(DATA.TemplateScores(:,p(1)),DATA.TemplateScores(:,p(2)),angle);
            cluster.space = [3 p(1:2)];
        elseif E.space(1) == 10
            xy = xyrotate(DATA.icas(:,p(1)),DATA.icas(:,p(2)),angle);
            cluster.space = [E.space(1) p(1:2)];
        else
            xy = xyrotate(DATA.pcs(:,p(1)),DATA.pcs(:,p(2)),angle);
            cluster.space = [1 p(1:2)];
        end
        if forcesign == 1
            id = find(xy(DATA.uid,1) > crit);
            nid = find(xy(DATA.uid,1) <= crit);
            cluster.sign = 1;
        elseif forcesign == -1
            id = find(xy(DATA.uid,1) < crit);
            nid = find(xy(DATA.uid,1) >= crit);
            cluster.sign = -1;
        elseif crit < 0 && DATA.plottype == 1
            id = find(xy(DATA.uid,1) < crit);
            nid = find(xy(DATA.uid,1) >= crit);
            cluster.sign = -1;
        elseif crit <0 && mean(xy(:,1)) < 0
            id = find(xy(DATA.uid,1) < crit);
            nid = find(xy(DATA.uid,1) >= crit);
            cluster.sign = -1;
        else
        id = find(xy(DATA.uid,1) > mean(exy(:,1)));
        nid = find(xy(DATA.uid,1) <= mean(exy(:,1)));
            cluster.sign = 1;
        end
        aid = id;
        anid = nid;
        id = DATA.uid(id);
        nid = DATA.uid(nid);

        e(1) = mean(DATA.energy(1,id));
        e(2) = mean(DATA.energy(1,nid));
        if e(2) > e(1) && forcesign == 0
            cluster.sign = -cluster.sign;
            cid = nid;
            nid = id;
            id = cid;
            cid = anid;
            anid = aid;
            aid = cid;
        end
        if  length(id) < 2
            
        end
        if quickmode.quickest == 1
            cluster.dprime = NaN;
            cluster.hdip = NaN;
            cluster.dipsize = NaN;
            cluster.autodipsize = NaN;
            cluster.gmdprime = 0;
        else
%if quck ==2, do the 1d fit, so we can see fit quality             
           if quickmode.fit1cut == 1 || quickmode.quick == 0
            [a,b] = GMDip(xy(DATA.uid,:),DATA.energy(1,DATA.uid),'crit',crit,'label',DATA.idstr);
            cluster.autodipsize = b.dipsize;
            cluster.dipsize = b.cdipsize;
            cluster.gmdprime = b.gmdprime;
            cluster.gmfit1d = b.G{b.best};
           end
            if quickmode.quick == 0
                cluster.dprime = CalcDprime(xy(id,1),xy(nid,1));
                cluster.hdip = HartigansDipTest(sort(xy(:,1)));
                x = FitGaussMeans(xy(DATA.uid,:),2,'clusterid',aid);
            end
        end
        r = xy(:,1);
    else  %E.shape == 0   = ellpse
        cluster.shape = 0;
        xy = zeros(length(DATA.clst),2);
        if E.space(1) == 6
            xy(DATA.uid,:) = DATA.ndxy(DATA.uid,:);
            cluster.space = E.space;
        elseif isempty(p) %Var/E plot
            xy(DATA.uid,2) = DATA.spkvar(ispk,DATA.uid)./DATA.energy(1,DATA.uid);
            xy(DATA.uid,1) = DATA.energy(1,DATA.uid);
            cluster.space = [5];
        elseif ismember(E.space(1),[3 4])
            xy(DATA.uid,1) = DATA.TemplateScores(DATA.uid,p(1));
            xy(DATA.uid,2) = DATA.TemplateScores(DATA.uid,p(2));
            cluster.space = [3 p(1:2)];
        elseif E.space(1) == 2
            AllV = GetDataStruct(DATA, 'AllV');
            if DATA.verbose > 1
                fprintf('%s AllV %s(%d)\n',IDStr(DATA),sprintf('%d ',size(AllV)),DATA.toplevel);
            end
            xy(DATA.uid,1) = AllV(p(1),p(2),DATA.uid);
            xy(DATA.uid,2) = AllV(p(3),p(4),DATA.uid);
            cluster.space = [2 p(1:4)];
        elseif E.space(1) == 10
            xy(DATA.uid,1) = DATA.icas(DATA.uid,p(1));
            xy(DATA.uid,2) = DATA.icas(DATA.uid,p(2));
            cluster.space = [E.space(1) p(1:2)];

        else
            xy(DATA.uid,1) = DATA.pcs(DATA.uid,p(1));
            xy(DATA.uid,2) = DATA.pcs(DATA.uid,p(2));
            cluster.space = [1 p(1:2)];
        end
        if isfield(E,'angle')
            cluster.angle = E.angle;
        else
            cluster.angle = 0;
        end
        r = CalcClusterDistance(cluster, xy);
        id = find(r(DATA.uid) < 1);
        nid = find(r(DATA.uid) >1);
        if quickmode.quick == 0
        x = FitGaussMeans(xy(DATA.uid,:),2,'clusterid',id);
        end
        if quickmode.quick == 0 || quickmode.fit1cut
            if cluster.shape == 0
                [a,b] = GMDip(Rprime(r(DATA.uid)),DATA.energy(1,DATA.uid),'crit',1,'label',DATA.idstr);
            else
                [a,b] = GMDip(r(DATA.uid),DATA.energy(1,DATA.uid),'crit',1,'label',DATA.idstr);
            end
            cluster.autodipsize = b.dipsize;
            cluster.dipsize = b.cdipsize;
            cluster.gmdprime = b.gmdprime;
            cluster.gmfit1d = b.G{b.best};
        end
        id = DATA.uid(id);
       nid = DATA.uid(nid);
    end
    
    %first reset any previous events to cluster 0
    cl.clst(cl.clst == clnum) = 1;
    %nid and clid refer to posttions in whole array, not just in uid
    if length(DATA.uid) < length(cl.clst)
        DATA.nid = nid;
        DATA.clid = id;
        X = setdiff(1:length(cl.clst),DATA.uid);
        xid = find(cl.clst(X) > 0);
        cl.clst(X(xid)) = 0;
        cl.clst(id) = clnum;
        if cluster.cluster == 1 %why did we do this - unsets clsuter 2 if its tehre
            %            cl.clst(nid) = 1;
        end
    else
        cl.clst(id) = clnum;
        DATA.nid = nid;
        DATA.clid = id;
    end
    DATA.clst = cl.clst;
    cl.id = DATA.clid;
    cl.nid = DATA.nid;
    
    cluster.ncut = length(id);
    cluster.clst = cl.clst;
    if cluster.space(1) == 3
        cluster.spacescale = DATA.TemplateScaling(cluster.space(2:end));
    end
    
    if quickmode.quick  
        cluster.bmc = 0;
        cluster.mahal = [0 0 0 0];
        if quickmode.fit2gauss
            [dp, fits, dpdetails] = Fit2Gauss(cluster, r, DATA);
            cluster.fitdprime = [dp dpdetails.fitpos dpdetails.dx dpdetails.minxpt];
            cluster.fitdpparams = cat(1,fits{1}.params, fits{2}.params);
        end
    else
        [dp, fits, dpdetails] = Fit2Gauss(cluster, r, DATA);
        cluster.fitdprime = [dp dpdetails.fitpos dpdetails.dx dpdetails.minxpt];
        cluster.fitdpparams = cat(1,fits{1}.params, fits{2}.params);
        cluster.bmc = BimodalCoeff(r,1.5);
        dip = MyDip(r);
        cluster.mydipsize = dip.dipsize;
        cluster.mydip = dip.x(dip.dip);;
        a = FitGaussMeans(xy(DATA.uid,:),2);
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
        if isfield(cluster,'gmdprime')
            cluster.mahal(4) = cluster.gmdprime;
        elseif isfield(E,'gmdprime')
            cluster.mahal(4) = E.gmdprime;
        elseif isfield(E,'gmfit1d') && isfield(E.gmfit1d,'mu')
            cluster.mahal(4) = GMdprime(E.gmfit1d);
        else
            cluster.mahal(4) = 0;
        end
        cluster.times = DATA.t(id);
        cluster.quick = 0;
    end



%don't rebuild MeanSpike if quick. But if there is no meanspike, it needs
%buidling
    if quickmode.quick == 0 || ~isfield(C,'MeanSpike') || quickmode.mean
        cl.MeanSpike = PlotMeanSpike(DATA,'recalc');
        cluster.MeanSpike = cl.MeanSpike;
        if DATA.plot.xyseq
            C = cluster;
            C.clst = cl.clst;
            C.xy = xy;
            PlotXYSequence(DATA,C);
        end
    end
    if DATA.check.dropi(1) > 0
        quickmode.dropi = 1;
        quickmode.fit2gauss = 1;
    end
    cluster.quick = quickmode.quick; %if this is set, need to call again
    cluster.r = r;
    if quickmode.triggerhist
        cluster =  PlotTriggerHist(DATA,cluster,quickmode);
    end
    cluster.minspke = prctile(DATA.energy(1,cl.id),1) .* 0.95;
    cluster.minspkvar = prctile(DATA.spkvar(DATA.probe(1),cl.nid),1) .* 0.95;
    if recluster == 0  %If cluster  boundary has changed, will need new template
%if reclsuter > 0 means this is applying an existing cluster, so if template has been 
%calcualted, no need to set this.
        cluster.neednewtemplate = 1;
    end
    if isfield(DATA,'TemplateUsed') && ~isempty(DATA.TemplateUsed) 
        cluster.TemplateUsed = DATA.TemplateUsed;
        if isfield(DATA,'DprimeUsed')
            cluster.DprimeUsed = DATA.DprimeUsed;
        end
        if isfield(DATA,'mumeanUsed')
            cluster.mumeanUsed = DATA.mumeanUsed;
        end
    end
    if cluster.cluster > 1
        clnum = cluster.cluster;
        a = rmfields(cluster,'next'); %avoid recursion
        cluster = DATA.cluster;
        cluster.next{clnum-1} = a;
    end
 
    E.h = [];
    if DATA.watchplots  
        if replot 
        DATA = ReplotPCs(DATA,cluster);
        SetFigure(DATA.tag.spikes, DATA);
        hold off;
        DATA.spkst = DATA.uid;
        QuickSpks(DATA, 1000);
        if (cluster.space(1) == 3) && ...
                size(DATA.TemplateScores,2) > 12 
            SetFigure(DATA.tag.tmplscore, DATA);
            subplot(1,1,1);
            hold off;
            plot(DATA.TemplateScores(cl.nid,1),DATA.TemplateScores(cl.nid,13),'.');
            hold on;
            plot(DATA.TemplateScores(cl.id,1),DATA.TemplateScores(cl.id,13),'r.');
        elseif cluster.space(1) == 6
            drawnow;
            SetFigure(DATA.tag.tmplscore, DATA);
            subplot(1,1,1);
            PlotXY(DATA.ndxy,DATA.clst);
            set(gca,'UserData',[NaN cluster.space]);
        end
        end
        if DATA.plot.isi
            PlotISI(DATA,0,0);
        end
    end
    if DATA.auto.checkxcorr
        cls = unique(DATA.clst);
        if length(cls) > 2
            CalcXcorr(DATA,[],'clusters');
        end
    elseif DATA.plot.xcorrprobes
        CalcXcorr(DATA,[],'probes');
    end

function QuickSpks(DATA,nspk)
    step = max([1 round(length(DATA.uid)/nspk)]);
    id = 1:step:length(DATA.uid);
    PlotSpikes(DATA,DATA.uid(id),'setcontext');
    FinishSpikePlot(DATA);

    function cluster = PlotTriggerHist(DATA, cluster,varargin)
        j = 1;
        calcdropi = 1;
        plotallcl = 0;
        while j <= length(varargin)
            if isfield(varargin{j}, 'quickest')
                quickmode = varargin{j};
                if quickmode.quick && quickmode.dropi == 0
                    calcdropi = 0; 
                end
            elseif strncmpi(varargin{j},'quick',5)
                calcdropi = 0;
            elseif strncmpi(varargin{j},'showall',6)
                plotallcl = 1;
            end
            j = j+1;
        end        
         
 %If call this for differnt clusters, DATA.nid may not be right list
 if 0
     nid = DATA.nid;
     id = DATA.clid;
 else
     id = find(DATA.clst == cluster.cluster+1);
     nid = find(DATA.clst ~= cluster.cluster+1);
 end
 if DATA.interactive >=0 
     SetFigure(DATA.tag.vhist, DATA);
     subplot(1,1,1);
     hold off;
 end
 t = find(DATA.spts == 0);
 if length(id) < 10
     nbins(1) = 1;
     cluster.dropi = [0 0 NaN 0];
     cluster.trigsd = 0;
     cluster.minspke = min(DATA.energy(1,DATA.clid));
     cluster.minspkvar = min(DATA.spkvar(DATA.probe(1),DATA.clid));
 elseif length(id) < 100
     nbins(1) = 10;
 elseif length(id) < 5000
     nbins(1) = round(length(DATA.clid)./20);
 else
     nbins(1) = 250;
 end
 if length(nid) < 10
     nbins(2) = 1;
 elseif length(nid) < 100
     nbins(2) = 10;
 elseif length(nid) < 5000
     nbins(2) = round(length(nid)./20);
 else
     nbins(2) = 250;
 end
     
 if nbins(1) > 2 && ~isempty(id)
%     V = AllV(DATA.probe(1),t,id);
     V = DATA.rV(id);
     [a,b] = hist(V,nbins(1));
     cluster.vhist = a;
     cluster.vhistrange = minmax(b);
     [c,d] = hist(DATA.rV(nid),nbins(2));
     if DATA.interactive >= 0
         bar(b,a,1,'facecolor',DATA.colors{cluster.cluster+1},'edgecolor',DATA.colors{cluster.cluster+1});
         hold on;
         plot(d,c .*max(a)./max(c),'color',[0.5 0.5 0.5]);
     end
     cluster.muvhist = c;
     cluster.muvhistrange = minmax(d);
     if DATA.Trigger(1) < 0
         tid = find(b <= DATA.Trigger(1));
         ntid = find(d <= DATA.Trigger(1));
     else
         tid = find(b >=DATA.Trigger(1));
         ntid = find(d >= DATA.Trigger(1));
     end
     if calcdropi
         gfit = FitGauss(b(tid),a(tid),'maxiter',500);
         if nbins(2) > 1
             ngfit = FitGauss(d(ntid),c(ntid),'maxiter',500);
         else
             ngfit.sd = 0;
             ngfit.mean = d;
         end
         if length(tid) && DATA.interactive >= 0
             hold on;
             plot(b(tid),gfit.fitted,'k');
         end
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
     elseif ~isfield(cluster,'dropi')
         cluster.dropi = [0 0 0 0];
     end
 
     if plotallcl
         cls = unique(DATA.clst);
         for j = 3:length(cls)
             id = find(DATA.clst == cls(j));
              V = DATA.rV(id);
              [a,b] = hist(V,nbins(1));
              plot(b,a,'color',DATA.colors{j});
         end
     end

if DATA.interactive >= 0
    p = ProbeNumber(DATA);
     title(sprintf('P%d/%d %d/%d (%.0f/%.0fHx) spikes drop %.3f,%.2f',p,DATA.currentcluster,...
         length(id),length(id)+length(nid),length(id)/DATA.duration,...
         (length(id)+length(nid))/DATA.duration,cluster.dropi(1),cluster.dropi(3)));
     hold on;
     plot([DATA.Trigger(1) DATA.Trigger(1)],get(gca,'ylim'),'r');
end
 elseif nbins(2) > 1 && DATA.interactive >= 0
    V = DATA.rV(nid);
     [c,d] = hist(V,nbins(2));
     hold off;
     plot(d,c ,'color',[0.5 0.5 0.5]);

     title(sprintf('P%d/%d %d/%d (%.0f/%.0fHx) spikes drop %.3f,%.2f',ProbeNumber(DATA),DATA.currentcluster,...
         length(id),length(id)+length(nid),length(id)/DATA.duration,...
         (length(id)+length(nid))/DATA.duration,cluster.dropi(1),cluster.dropi(3)));
     hold on;
     plot([DATA.Trigger(1) DATA.Trigger(1)],get(gca,'ylim'),'r');
     
 end

function chspk = UseProbeList(DATA, nprobes)
    
    nx = floor(nprobes/2);
    chspk = [-nx:nx] + DATA.probe(1);
    chspk = chspk(chspk > 0 & chspk <= DATA.nprobes);
    
            

function MeanSpike = PlotMeanSpike(DATA, varargin)
    recalc = 0;
    clnum = DATA.currentcluster;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cluster',4)
            j = j+1;
            clnum = varargin{j};
        elseif strncmpi(varargin{j},'recalc',4)
            recalc = 1;
        end
        j = j+1;
    end
    id = find(DATA.clst == clnum+1);
    nid = find(DATA.clst == 1);
    MeanSpike.addmean = DATA.addmean;
    if DATA.interactive >= 0
        SetFigure(DATA.tag.meanspike, DATA);
        subplot(1,2,1);
        hold off;
        AllV = GetAllV(DATA);
        MeanV = getappdata(DATA.toplevel,'MeanV');
    else
        AllV = DATA.AllV;
        MeanV = DATA.meanV;
    end
    if DATA.keepsmooth && ~isempty(id)
        lfpts = -20000:10:20000;
        sampleid = DATA.trigtimes{1}(id);
        samplenid = DATA.trigtimes{1}(nid);
        Vall = getappdata(DATA.toplevel,'Vall');
        allid = RectAdd(lfpts, sampleid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.lfpmean = mean(Vall.Vsmooth(allid)+double(Vall.V(allid)),1);
        [st, tid] = ShuffleTimes(Vall.t(sampleid),DATA.Expt);
        sid  = FullVTimes2id(Vall,st);
        allid = RectAdd(lfpts, sid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.shufflemean = mean(Vall.Vsmooth(allid)+double(Vall.V(allid)),1);
        muid = RectAdd(lfpts, samplenid, 'range', [1 length(Vall.Vsmooth)]);
        MeanSpike.lfpumean = mean(Vall.Vsmooth(muid)+double(Vall.V(muid)),1);
    end
    if ~isfield(DATA,'MeanSpike') || recalc
        ms = mean(AllV(:,:,id),3);
        mu = mean(AllV(:,:,nid),3);
        if DATA.addmean == 0 %if 1, already added to AllV
            mum =mean(MeanV(:,nid)');
            msm = mean(MeanV(:,id)');
            for j = 1:DATA.nprobes
                ms(j,:) = ms(j,:) + msm;
                mu(j,:) = mu(j,:) + mum;
            end
        end
        MeanSpike.ms = ms;
        MeanSpike.mu = mu;
        xc = corrcoef(ms(:),mu(:));
        MeanSpike.muxc = xc(1,2);
    else
        if isfield(DATA.cluster,'MeanSpike')
            ms = DATA.cluster.MeanSpike.ms;
            mu = DATA.cluster.MeanSpike.mu;
        else
        ms = DATA.MeanSpike.ms;
        mu = DATA.MeanSpike.mu;
        end
    end

    chspk = DATA.probe(1)+ [-1:1]; 
    chspk = DATA.chspk;
    dj = 0;
%do csd first, then can do dvdt to CSD 
    if DATA.plotcsd
       ms = (diff(ms,2,1));
       mu = (diff(mu,2,1));
       csd = diff(AllV,DATA.plotcsd,1);
        chspk = chspk-1;
        dj = 1;
    end
    if DATA.plotdvdt
       ms = (diff(ms,1,2));
       mu = (diff(mu,1,2));
    end
    if DATA.interactive >= 0
    imagesc(ms);
    title(sprintf('P%d Cl %d',DATA.probe(1),clnum));
    if DATA.nprobes == 1
        if DATA.keepsmooth && isfield(MeanSpike,'lfpmean')
            subplot(2,1,2);
            hold off;
            plot(MeanSpike.lfpmean);
            hold on;
            plot(get(gca,'xlim'),[0 0],'k:');
            plot(MeanSpike.lfpumean,'g');
            plot(MeanSpike.shufflemean,'r');
            subplot(2,1,1);
        else
            subplot(1,1,1);
        end
    else
        subplot(1,2,2);
        if 0  %might add this one day
            hold off;
            imagesc(MeanSpike.mu);
            return;
        end
    end
    
    

    hold off;
    end
    chspk = chspk(chspk >0 & chspk <= size(ms,1));
    voff = DATA.voffset - DATA.voffset(DATA.probe(1));
    for j = chspk
        if DATA.plotdvdt && DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(csd(j,:,id),1,2),[],3) var(diff(csd(j,:,nid),1,2),[],3)]));
        elseif DATA.plotdvdt
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(diff(AllV(j,:,nid),1,2),[],3) var(diff(AllV(j,:,id),1,2),[],3)]));
        elseif DATA.plotcsd
                dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(csd(j,:,id),[],3) var(csd(j,:,nid),[],3)]));
        else
        dp = (ms(j,:)-mu(j,:))./sqrt(mean([var(AllV(j,:,nid),[],3) var(AllV(j,:,id),[],3)]));
        end
        dpscale = max(ms(:))./max(dp);
        if DATA.interactive >= 0
            plot(mu(j,:)+voff(j)/5,'color',[0.5 0.5 0.5]);
            hold on;
            plot(ms(j,:)+voff(j)/5,'r');
            if j == DATA.probe(1)
                plot(abs(dp).*dpscale,'m');
            else
                plot(abs(dp),'g');
            end
            text(size(ms,2),voff(j)/5,sprintf('%d',j+dj),'fontsize',DATA.gui.fontsize(1));
        end
        peaks = find(diff(sign(diff(abs(dp)))) > 0);
        [a,b] = sort(abs(dp(peaks)),'descend');
        MeanSpike.dpmax(j,1:length(b)) = peaks(b);
        MeanSpike.dp(j,:) = dp;
        if recalc
        MeanSpike.vdprime(j,:) = dp;
        end
    end
    if DATA.interactive >= 0
    set(gca,'xlim',[1 size(ms,2)+1]);
    end

function [Scores, T, details] = CalcScores(DATA, MeanSpike)
% this need to be modified to work like TemplatePlot now, with
%T being for all probes, but only calculating scores for chspk
 
details.tstart = now;

    AllV = GetAllV(DATA);

if isfield(MeanSpike,'ms')
    T = MeanSpike.ms;
    if isfield(MeanSpike,'mu')
        mT = MeanSpike.mu;
    else
        mT = T;
    end
else
    T = MeanSpike;
end

if isfield(MeanSpike,'othermeans')  && sum(size(MeanSpike.othermeans{1}) == size(MeanSpike.ms)) > 1
    oT = MeanSpike.othermeans;
else
    oT{1} =T;    
end

    if DATA.chspk(end) == DATA.nprobes
    else
        dspk = [DATA.chspk DATA.chspk(end) + 1];
    end
    if max(DATA.chspk) >= DATA.nprobes
        dspk = [DATA.nprobes-length(DATA.chspk):1:DATA.nprobes];
        csdspk = [DATA.nprobes-length(DATA.chspk)-1:1:DATA.nprobes];
    elseif min(DATA.chspk) <= 1
        dspk = [1:length(DATA.chspk)+1];
        csdspk = [1:length(DATA.chspk)+2];
    else
        dspk = [DATA.chspk DATA.chspk(end) + 1];
        csdspk = [DATA.chspk(1) - 1 DATA.chspk DATA.chspk(end) + 1];
    end
    
    
    if DATA.usestdtemplates
        j = DATA.probe(1);
        ispk = j;
        meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]);
        Scores(1,1,:) = MeanSpike(1,:) * squeeze(AllV(j,:,:) - meanV);
        Scores(1,2,:) = MeanSpike(2,:) * squeeze(AllV(j,:,:) - meanV);
        Scores(1,3,:) = squeeze(diff(DATA.StdTemplate(1,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        Scores(1,4,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        Scores(1,5,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        return;
    end

    if size(AllV,1) > 2
        if size(AllV,1) < max(csdspk) && DATA.loadfromspikes
            AllV(max(csdspk),end,end) = 0;
        end
    dv = diff(AllV(dspk,:,:),1,1);
    mdv = diff(T(dspk,:,:),1,1);
    csd = diff(AllV(csdspk,:,:),2,1);
    mcsd = diff(T(csdspk,:,:),2,1);
    end
%    meanV = repmat(mean(AllV(DATA.chspk,:,:),2),[1  size(AllV,2) 1]);
    vid = 1:min([size(T,2) size(AllV,2)]); %% in case old cut used different range
    %but using indices makes calculatins slower, so only do this if
    %necessary.
    if size(T,2) ~= size(AllV,2)
        sizediff = 1;
    else
        sizediff = 0;
    end
%does meanV need to know about vid?    
    for j = 1:length(DATA.chspk)
        c = DATA.chspk(j);
        meanV = repmat(squeeze(mean(AllV(c,:,:),2))',[size(AllV,2) 1]);
        AllVs{j} = squeeze(AllV(c,:,:)) - meanV;
    end
    
    Scores = zeros(length(DATA.chspk),7+length(oT),DATA.nevents);
    for j = 1:length(DATA.chspk)
        c = DATA.chspk(j);
        if size(AllV,1) > 2
        Scores(j,3,:) = squeeze(mdv(j,:)) * squeeze(dv(j,:,:));
        Scores(j,4,:) = squeeze(mcsd(j,:)) * squeeze(csd(j,:,:));
        end
        if sizediff
            sV = T(c,vid) - mean(T(c,vid));
            mV = mT(c,vid) - mean(mT(c,vid));
            Scores(j,2,:) = squeeze(diff(T(c,vid),1,2)) * squeeze(diff(AllV(c,vid,:),1,2));
            for k = 1:length(oT)
                oV = oT{k}(c,vid) - mean(oT{k}(c,vid));
                Scores(j,k+7,:) =  squeeze(oV) * AllVs{j};
            end
        else
            sV = T(c,:) - mean(T(c,:));
            mV = mT(c,:) - mean(mT(c,:));
            Scores(j,2,:) = squeeze(diff(T(c,:),1,2)) * squeeze(diff(AllV(c,:,:),1,2));
            for k = 1:length(oT)
                oV = squeeze(oT{k}(c,:) - mean(oT{k}(c,:)));
                Scores(j,k+7,:) =  oV * AllVs{j};
            end
        end
        Scores(j,1,:) = squeeze(sV) * AllVs{j};
        Scores(j,6,:) =  squeeze(mV) * AllVs{j};
        Scores(j,7,:) = sum(abs(repmat(sV',1,size(AllV,3)) - AllVs{j}));
        if size(MeanSpike.vdprime,1) >= c
        dp = MeanSpike.vdprime(c,:);
        Scores(j,5,:) = dp * squeeze(AllV(c,:,:) - repmat(mV,[1 1 DATA.nevents]));
        else
        Scores(j,5,:) = 0;
        end
    end
    if length(DATA.chspk) == 1
        meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]);
        Scores(1,3,:) = DATA.StdTemplate(1,:) * squeeze(AllV(DATA.chspk,:,:) - meanV);
        Scores(1,4,:) = DATA.StdTemplate(2,:) * squeeze(AllV(DATA.chspk,:,:) - meanV);
        Scores(1,6,:) = MeanSpike.mu * squeeze(AllV(DATA.chspk,:,:) - meanV);;
        Scores(1,7,:) = DATA.rV;
        Scores(1,8,:) = squeeze(diff(DATA.StdTemplate(1,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
    end
    if 0 %old way. Might want this back. Now store all channels (so can use CSD of template)
        T = T(DATA.chspk,:);
    end
    details.endtime = now;


function Labels = TemplateLabels(DATA, usestd)
    ispk = DATA.probe(1);
    if isfield(DATA,'chspk')
        chspk = DATA.chspk;
    else
        chspk = DATA.probe(1)+ [-1:1];
        chspk = chspk(chspk >0 & chspk <= DATA.allnprobes);
    end
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
    if DATA.nprobes == 1
        ispk = DATA.probelist(DATA.probe(1));
        Labels{1} = sprintf('%d:r',ispk);
        Labels{2} = sprintf('MU');
        Labels{8} = sprintf('%d:dt',ispk);
        Labels{9} = sprintf('Std1');
        Labels{4} = sprintf('std2');
        Labels{6} = sprintf('%d:dp',ispk);
        Labels{5} = sprintf('%dr-otherr',ispk);
        Labels{7} = sprintf('%d:dp',chspk(1));
        Labels{3} = sprintf('Std1',chspk(1));
        Labels{10} = sprintf('dp wieghted');
        Labels{11} = sprintf('trigger');
        Labels{12} = sprintf('dvdt Std1',ispk);
        Labels{14} = sprintf('MUsum',ispk);
        Labels{15} = sprintf('sum2',ispk);
        Labels{13} = sprintf('absdiff',ispk);
        Labels{16} = sprintf('sum2-sum1');
        Labels{17} = sprintf('sumP%d',DATA.probe(1)-1);
        Labels{18} = sprintf('sumP%d',DATA.probe(1)+1);
    else
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
            if max(chspk)  == DATA.nprobes
                xspk = min(chspk)-1;
            else
                xspk = max(chspk)+1;
            end
            Labels{4} = sprintf('%d:r',xspk);
        end
        Labels{10} = sprintf('sumdt');
        if DATA.trigdt == 4 %template triggering - look at trig values
            Labels{11} = sprintf('trigger');
        else
            Labels{11} = sprintf('sumdp',ispk);
        end
        Labels{12} = sprintf('sumdy',ispk);
        Labels{14} = sprintf('MUsum',ispk);
        Labels{15} = sprintf('sum2',ispk);
        Labels{13} = sprintf('absdiff',ispk);
        Labels{16} = sprintf('sum2-sum1');
        Labels{17} = sprintf('sumP%d',DATA.probe(1)-1);
        Labels{18} = sprintf('sumP%d',DATA.probe(1)+1);
    end
    
function [out, TemplateUsed, DprimeUsed] = TemplatePlot(DATA, varargin)

    projectout = 0;
    calcdips = 0;
    usemean = 0;
    usestd = 0;
    %normalized correlation problem is htat get lots close to 1 (or
    %std of template, so its compresed at that end. Not good for GM fits
    normalize = 0;
    plottype = 1;
    if DATA.profiling
        ts = now;
        profile on;
    end
    j = 1;
    while  j <= length(varargin)
        if strncmpi(varargin{j},'calcdip',7)
            calcdips = 1;
        elseif strncmpi(varargin{j},'nodip',4)
            calcdips = 0;
        elseif strncmpi(varargin{j},'noplot',6)
            plottype = 0;
        elseif strncmpi(varargin{j},'projectout',7)
            projectout = 1;
        elseif strncmpi(varargin{j},'recalc',6) %redo cal with 'templateused field
            usemean = 2;
        elseif strncmpi(varargin{j},'usemean',7)
            usemean = 1;
        elseif strncmpi(varargin{j},'stdtemplate',7)
            usestd = 1;
        end
        j = j+1;
    end

    AllV = GetAllV(DATA);
  
    DataClusters = getappdata(DATA.toplevel,'Clusters');
    %default in case there is no template to use
    othermeans{1} = repmat(DATA.StdTemplate(1,1:length(DATA.cluster.spts)),DATA.nprobes,1);
    othermeans{2} = repmat(DATA.StdTemplate(2,1:length(DATA.cluster.spts)),DATA.nprobes,1);
    othermeans{3} = repmat(DATA.StdTemplate(2,1:length(DATA.cluster.spts)),DATA.nprobes,1);
    if DATA.currentcluster > 1 && isfield(DATA.cluster,'next') && length(DATA.cluster.next) >= DATA.currentcluster-1
        C = DATA.cluster.next{DATA.currentcluster-1};
        C.cluster = DATA.currentcluster;
        if isfield(DATA.cluster,'MeanSpike')
            othermeans{1} = DATA.cluster.MeanSpike.ms;
        end
    else
        C = DATA.cluster;
        C.cluster = 1;
        if length(DATA.cluster.next) && isfield(DATA.cluster.next{1},'MeanSpike')
        othermeans{1} = DATA.cluster.next{1}.MeanSpike.ms;
        end
    end
    if C.probe(1) >1 && length(DataClusters) >= C.probe(1) && isfield(DataClusters{C.probe(1)-1},'MeanSpike')
        othermeans{2} = DataClusters{C.probe(1)-1}.MeanSpike.ms;
    end
    if ~isfield(C,'quick')
        C.quick = 0;
    end
    if isfield(C,'neednewtemplate') && C.neednewtemplate == 0
        SetFigureName(DATA.toplevel,sprintf('Calculating Template Scores for C%d...',DATA.currentcluster));
    else
        SetFigureName(DATA.toplevel,sprintf('Calculating New Templates for C%d...',DATA.currentcluster));
    end
    
    if C.quick || ~isfield(C,'MeanSpike')
        C.MeanSpike = PlotMeanSpike(DATA,'recalc');
        if C.cluster > 1
            DATA.cluster.next{C.cluster-1}= C;
        else
            DATA.cluster = C;
        end
    end
    DATA.usestdtemplates = usestd;
    id = find(DATA.clst == DATA.currentcluster+1);
    nid = find(DATA.clst == 1);
    nevents = DATA.nevents;
    %If a new cluster has stolen points from an old one, want to keep teh
    %old template scores for that old cluster (until its boudnary is
    %chaned)
    if isfield(C,'neednewtemplate') && C.neednewtemplate == 0
        SetFigureName(DATA.toplevel,sprintf('Calculating Template Scores for C%d...',DATA.currentcluster));
        ms = C.TemplateUsed;
        mu = C.MeanSpike.mu;
        if length(DATA.chspk) == size(C.mumeanUsed,1)
            mu(DATA.chspk,:) = C.mumeanUsed;
%        elseif max(DATA.chspk) == size(C.mumeanUsed,1)
%            mu(DATA.chspk,:) = C.mumeanUsed;
        end
    elseif usemean == 1 || isempty(id)
        SetFigureName(DATA.toplevel,sprintf('Calculating Templates from mean C%d...',DATA.currentcluster));
        ms = C.MeanSpike.ms;
        mu = C.MeanSpike.mu;
    else
        SetFigureName(DATA.toplevel,sprintf('Calculating New Templates for C%d...',DATA.currentcluster));
        for j = DATA.nprobes:-1:1
            ms(j,:) = mean(AllV(j,:,id),3);
            mu(j,:) = mean(AllV(j,:,nid),3);
        end
    end
    if size(mu,1) < size(ms,1)
        a = size(mu,1);
        b = size(ms,1);
        mu(a+1:b,:) = 0;
    end
    if C.cluster > 1
        DATA.cluster.next{C.cluster-1}.neednewtemplate = 0;
    else
        DATA.cluster.neednewtemplate = 0;
    end

   for j = DATA.nprobes:-1:1
            mu(j,:) = mu(j,:)-mean(mu(j,:));
            ms(j,:) = ms(j,:)-mean(ms(j,:));
            %        xc(j,1,:) = squeeze(TemplateScores(j,1,:))./(DATA.spkvar(j,:)' .* std(ms(j,:)));
            %        xc(j,2,:) = squeeze(TemplateScores(j,2,:))./(DATA.spkvar(j,:)' .* std(mu(j,:)));
   end

   if isnan(sum(ms(:)))
       fprintf('Cannot calculate template for empty cluster\n');
       DATA.MeanSpike.ms = zeros(size(ms));
       TemplateUsed = DATA.MeanSpike.ms;
       DprimeUsed = [];
       set(DATA.toplevel,'UserData',DATA);
       return;
   end
    ispk = find(DATA.chspk == DATA.probe(1));
    if usestd
        j = DATA.probe(1);
        meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]); %mean for each spike
        TemplateScores(ispk,1,:) = DATA.StdTemplate(1,:) * squeeze(AllV(j,:,:) - meanV);
        TemplateScores(ispk,2,:) = DATA.StdTemplate(2,:) * squeeze(AllV(j,:,:) - meanV);
        TemplateScores(ispk,3,:) = squeeze(diff(DATA.StdTemplate(1,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateScores(ispk,7,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateScores(ispk,8,:) = squeeze(diff(DATA.StdTemplate(2,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
        TemplateUsed = DATA.StdTemplate;
        DprimeUsed = [];
    else
        for k = length(DATA.chspk):-1:1
            j = DATA.chspk(k);
            meanV = repmat(mean(AllV(j,:,:),2),[1  size(AllV,2) 1]);
            TemplateScores(k,1,:) = squeeze(ms(j,:)) * squeeze(AllV(j,:,:) - meanV);
            TemplateScores(k,2,:) = squeeze(mu(j,:)) * squeeze(AllV(j,:,:) - meanV);
            for nt = 1:3
            if size(othermeans{nt},2) ==  size(AllV,2)
                TemplateScores(k,9+nt,:) = squeeze(othermeans{nt}(j,:)) * squeeze(AllV(j,:,:) - meanV);
            elseif size(othermeans{nt},2)>  size(AllV,2)
                TemplateScores(k,9+nt,:) = squeeze(othermeans{nt}(j,1:size(AllV,2))) * squeeze(AllV(j,:,:) - meanV);
            end
            end
            if length(id)
                dp(j,:) = (mean(AllV(j,:,id),3)-mean(AllV(j,:,nid),3))./sqrt(mean([var(AllV(j,:,nid),[],3); var(AllV(j,:,id),[],3)]));
            else
                dp(j,:) = DATA.cluster.MeanSpike.dp(j,:);
            end
            if normalize
                TemplateScores(k,1,:) = TemplateScores(k,1,:) ./ std(AllV(j,:,:));
                TemplateScores(k,2,:) = TemplateScores(k,1,:) ./ std(AllV(j,:,:));
            end
            TemplateScores(k,9,:) = sum(abs(repmat(ms(j,:)',1,nevents) -squeeze(AllV(j,:,:)-meanV)));
        end
        TemplateUsed = ms;
        DprimeUsed = dp(DATA.chspk,:);
        clear meanV;
        dpcrit = 1;
        for k = length(DATA.chspk):-1:1
            j = DATA.chspk(k);
            %        dp * squeeze(AllV(j,:,:) - repmat(mu(j,:),[1 1 DATA.nevents]));
            TemplateScores(k,3,:) = squeeze(dp(j,:)) * squeeze(AllV(j,:,:) - repmat(mu(j,:),[1 1 DATA.nevents]));
            id = find(abs(dp(j,:)) > dpcrit);
            if length(id)
                TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(AllV(j,id,:));
                TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(AllV(j,id,:));
            else
                [a,b] = sort(dp(j,:),'descend');
                id = b(1:2);
                TemplateScores(k,4,:) = squeeze(dp(j,id)) * squeeze(AllV(j,id,:));
                TemplateScores(k,5,:) = squeeze(sign(dp(j,id))) * squeeze(AllV(j,id,:));
            end
            TemplateScores(k,6,:) = squeeze(diff(ms(j,:),1,2)) * squeeze(diff(AllV(j,:,:),1,2));
            if normalize
                TemplateScores(k,6,:) = TemplateScores(k,6,:) ./ std(diff(AllV(j,:,:),1,2));
            end

        end
        if size(AllV,1) > 2 %look at spatial derivatives
            if max(DATA.chspk) < DATA.nprobes
                dyspk = [DATA.chspk DATA.chspk(end)+1];
            elseif min(DATA.chspk) <= 1
                dyspk = [1:length(DATA.chspk)+1];
            else
                dyspk = [DATA.nprobes-length(DATA.chspk):DATA.nprobes];
            end
            dyspk = dyspk(dyspk > 0);
            dv = diff(AllV(dyspk,:,:),1,1);  %diff along length
            mdv = diff(ms(dyspk,:),1,1);
            for k = 1:size(dv,1)
                TemplateScores(k,7,:) = squeeze(mdv(k,:)) * squeeze(dv(k,:,:));
                if normalize
                    TemplateScores(k,7,:) = TemplateScores(k,7,:) ./ std(squeeze(mdv(k,:)) * squeeze(dv(k,:,:)));
                end

            end
            clear dv;
            clear mdv;

            if min(dyspk) > 1
                csdspk = [dyspk(1)-1 dyspk];
            else
                csdspk = [1:length(dyspk)+1];
            end
            csdspk = csdspk(csdspk > 0 & csdspk <= DATA.nprobes);

            csd = diff(AllV(csdspk,:,:),2,1);
            mcsd = diff(ms(csdspk,:),2,1);
            for j = 1:size(csd,1)
                TemplateScores(j,8,:) = squeeze(mcsd(j,:)) * squeeze(csd(j,:,:));
            end
            clear csd;
            clear mcsd;
        else %Temporary - fill 7 adn 8 with somethins
            if size(AllV,2) >  size(DATA.StdTemplate,2)
                vpts = 1:size(DATA.StdTemplate,2);
            else
                vpts = 1:size(AllV,2);
            end
            meanV = repmat(mean(AllV(DATA.chspk,:,:),2),[1  size(AllV,2) 1]); %mean for each spike
            TemplateScores(1,7,:) = DATA.StdTemplate(1,vpts) * squeeze(AllV(DATA.chspk,vpts,:) - meanV);
            TemplateScores(1,8,:) = DATA.StdTemplate(2,vpts) * squeeze(AllV(DATA.chspk,vpts,:) - meanV);
            TemplateScores(1,12,:) = squeeze(diff(DATA.StdTemplate(1,vpts),1,2)) * squeeze(diff(AllV(j,vpts,:),1,2));
        end
    end
    DATA.DprimeUsed = DprimeUsed;
    DATA.TemplateUsed = TemplateUsed;
    DATA.mumeanUsed = mu(DATA.chspk,:);
    DATA.Template.othermeans = othermeans;
    %These shouls be copied into cluster ONLY if this space is used for
    %classification, so don't it here
    if 0 
    if DATA.currentcluster > 1 && isfield(DATA.cluster,'next')
        DATA.cluster.next{DATA.currentcluster-1}.TemplateUsed = TemplateUsed;
        DATA.cluster.next{DATA.currentcluster-1}.DprimeUsed = DprimeUsed;
        DATA.cluster.next{DATA.currentcluster-1}.mumeanUsed = mu(DATA.chspk,:);
    elseif DATA.currentcluster == 1
        DATA.cluster.TemplateUsed = TemplateUsed;
        DATA.cluster.DprimeUsed = DprimeUsed;
        DATA.cluster.mumeanUsed = mu(DATA.chspk,:);
    end
    end
    if length(DATA.chspk) > 2
    chspk = DATA.chspk;
    else
    chspk = DATA.probe(1)+ [-1:1];
    end
    chspk = chspk(chspk >0 & chspk <= DATA.nprobes);
    if min(chspk) > 1
    csdspk = chspk-1;
    else
        csdspk = chspk;
    end
    
    if max(chspk)  == DATA.nprobes
        xspk = min(chspk)-1;
    else
        xspk = max(chspk)+1;
    end
        
    if projectout  %doesn't seem much use...., and uses too much memory
%project out the templates, then redo the pca
    mg = sum(ms.^2,2);
    G = sum(AllV.*repmat(ms,[1 1 DATA.nevents]),2)./repmat(mg,[1 1 DATA.nevents]);
    nv = AllV - repmat(ms,[1 1 DATA.nevents]) .* repmat(G,[1 size(AllV,2) 1]);
    TV = nv(chspk(1),:,:);
    for j = 2:length(chspk)
        TV = cat(2,TV,nv(chspk(j),:,:));
    end
    TV = squeeze(TV)';
    [pc, E] = eig(cov(TV));
    pc = fliplr(pc); %put largest first;
    pcs = TV*pc;
    end
    ispk = find(DATA.chspk == DATA.probe(1));

    if projectout == 2
        TMPL.pcs(:,1) = sum(TemplateScores(:,1,:));
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
    elseif length(DATA.chspk) == 1
        TMPL.pcs(:,1) = TemplateScores(ispk,1,:);
        TMPL.pcs(:,2) = TemplateScores(1,2,:);
        TMPL.pcs(:,8) = TemplateScores(ispk,6,:); %dvdt
        TMPL.pcs(:,3) = TemplateScores(ispk,7,:); %std 1
        TMPL.pcs(:,4) = TemplateScores(ispk,8,:); %std 2
        TMPL.pcs(:,6) = TemplateScores(ispk,3,:); %dprime weighted
        TMPL.pcs(:,7) = TemplateScores(1,3,:); %dprime weighted
        TMPL.pcs(:,3) = TemplateScores(1,7,:);
        TMPL.pcs(:,4) = TemplateScores(1,8,:);
        TMPL.pcs(:,10) = TemplateScores(1,3,:); %dprime
        TMPL.pcs(:,11) =DATA.rV;
        TMPL.pcs(:,12) = (TemplateScores(1,12,:)); %std1 dvdt
        TMPL.pcs(:,13) = (TemplateScores(1,9,:)); %sum abs diffs
        TMPL.pcs(:,14) = (TemplateScores(1,2,:)); %sum mu score
        TMPL.pcs(:,15) = (TemplateScores(1,10,:)); %sum score for other template
        TMPL.pcs(:,5) = (TemplateScores(1,10,:)) - TemplateScores(1,1,:); %diff in sums
        TMPL.pcs(:,17) = (TemplateScores(1,11,:)); %TemplateScore for cluster above
        TMPL.pcs(:,18) = (TemplateScores(1,12,:)); %TemplateScore for cluster below
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
        if DATA.trigdt == 4
            TMPL.pcs(:,11) =DATA.rV;
        else
            TMPL.pcs(:,11) = sum(TemplateScores(:,3,:)); %dprime weighted
        end
        TMPL.pcs(:,12) = sum(TemplateScores(:,7,:)); %sum dy
        TMPL.pcs(:,13) = sum(TemplateScores(:,9,:)); %sum abs diffs
        TMPL.pcs(:,14) = sum(TemplateScores(:,2,:)); %sum mu score
        TMPL.pcs(:,15) = sum(TemplateScores(:,10,:)); %sum score for other template
        TMPL.pcs(:,16) = sum(TemplateScores(:,10,:)) - sum(TemplateScores(:,2,:)); %diff in sums
        TMPL.pcs(:,17) = sum(TemplateScores(:,11,:)); %TemplateScore for cluster above
        TMPL.pcs(:,18) = sum(TemplateScores(:,12,:)); %TemplateScore for cluster below
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
        sds(j) = std(DATA.TemplateScores(:,j));
        if sds(j) > 0
        DATA.TemplateScores(:,j) = DATA.TemplateScores(:,j)./sds(j);
        DATA.TemplateScaling(j) = sds(j);
        end
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

    if DATA.watchplots  && plottype == 1
    SetFigure(DATA.tag.tmplscore, DATA);
    subplot(1,1,1);
    PlotTemplateScores(TMPL,TemplateScores, DATA.chspk);
%   plot(DATA.TemplateScores(:,1),DATA.spkvar(:,1));

    DATA = ReplotPCs(DATA,[]);
    end
%currently sets DATA.TemplateScores, DATA.TemplateLabels, and DATA.tmpdips
%if output is requestsed, don't set the figure userdata
    if DATA.profiling
       fprintf('Templates tookd %.2f\n',mytoc(ts));
       profile viewer;
    end
    if nargout
        set(DATA.toplevel,'Name',get(DATA.toplevel,'Tag'));
        out = DATA.TemplateScores;
        return;
    end
    set(DATA.toplevel,'UserData',DATA);
    SetFigureName(DATA.toplevel,get(DATA.toplevel,'Tag'));

    
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
        
setids = [];
meanpos = [];
autospace = 0;
args = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'autospace',6)
        autospace = 1;
    elseif strncmpi(varargin{j},'setid',5)
        j = j+1;
        setids = varargin{j};
    elseif strncmpi(varargin{j},'showpos',7)
        j = j+1;
        meanpos = varargin{j};
        j = j+1;
        dimpos = varargin{j};
    end
    j = j+1;
end
if DATA.interactive < 0
    return;
end
if autospace
    DATA.plottype = WhichPlotType(DATA.cluster, DATA.currentcluster);
end

SetFigure(DATA.tag.top, DATA);
clusterplot = [];
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
elseif DATA.plottype == 8;
    vpts = DATA.dvpts;
    AllV=getappdata(DATA.toplevel,'AllV');
    DATA.dV= diff(AllV,1,2);
elseif DATA.plottype == 9
    vpts = [1 1 2 2; 1 1 3 3];
    vpts = repmat(vpts,4,1);
end
%need to get boundary anyway in case there are other clusters
if isempty(E) && isfield(DATA,'cluster') && isfield(DATA.cluster,'space')
    if DATA.cluster.space(1) == DATA.plottype %cluster is in this space
        E = BoundaryFromCluster(E,DATA.cluster, DATA.currentcluster);
    else
        E = BoundaryFromCluster(E,DATA.cluster,DATA.currentcluster);
    end
end

plots = DATA.pcplots(1:8,:);
if DATA.plottype == 3 || DATA.plottype ==4
    if DATA.usestdtemplates
    plots = DATA.tmplots(17:24,:);
    elseif DATA.plottype == 4
        plots = DATA.tmplots(9:16,:);
    else
        plots = DATA.tmplots(1:8,:);
    end
    if isfield(DATA,'tmpdips')
    dipvals = DATA.tmpdips;
    end
elseif DATA.plottype == 12 %Try new scores
    tpt = find(DATA.spts == 0);
    p = DATA.probe(1);
    AllV = GetAllV(DATA);
    tvals = squeeze(AllV(p,tpt,:));
    pcs(:,1) = squeeze(AllV(p,tpt+6,:))-tvals;
    pcs(:,2) = squeeze(AllV(p,tpt+19,:))-tvals;
    pcs(:,3) = squeeze(AllV(p,tpt-6,:))-tvals;
    if p < size(AllV,1)
        pcs(:,4) = squeeze(AllV(p+1,tpt,:));
    else
        pcs(:,4) = squeeze(AllV(p-1,tpt+8,:));
    end
    if p > 0
        pcs(:,5) = squeeze(AllV(p-1,tpt,:));
    else
        pcs(:,5) = squeeze(AllV(p+1,tpt+8,:));
    end
    plots = DATA.pcplots;
elseif DATA.plottype == 11 %show spaces needed for all cells
    pcs = [];
    for j = 1:length(DATA.xy)
        pcs = cat(2,pcs,DATA.xy{j}(:,1));
        pcs = cat(2,pcs,DATA.xy{j}(:,2));
    end
    if size(pcs,2) == 4
        plots = [1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 1; 1 1]; 
    elseif size(pcs,2) == 3
        plots = [1 2; 1 3; 1 4; 2 3 ; 2 4; 3 4; 1 1; 1 1]; 
    elseif size(pcs,2) > 8
        for j = 1:2:size(pcs,2)
            plots(:,ceil(j/2)) = [j j+1];
        end
        k = 3;
        for j = (ceil(j/2)+1):8
            plots(:,j) = [1 k];
            k = k+1;
        end
    else
        plots = DATA.pcplots;
    end
elseif DATA.plottype == 7
    plots = DATA.tmplots(17:24,:);
elseif DATA.plottype == 5
    plots = DATA.tmplots(9:16,:);
end

   
if DATA.usegmcid & length(DATA.gmcid)
    if isfield(DATA,'gmcid')
        clid = DATA.gmcid;
    elseif isfield(E,'bestcl')
        clid = E.bestcl;
    else
    clid = DATA.clid;
    DATA.usegmcid = 0; %don't use it if not defined
    end
elseif length(setids)
    clid = setids;
    args = {args{:} 'ptsz' [8 1]};
else
    args = {args{:} 'ptsz' DATA.ptsz};
    if isfield(DATA,'clst') 
        clid = DATA.clst;
    else
        clid = DATA.clid;
    end
end
p = DATA.probelist(DATA.probe(1));
if DATA.clplot ==1 && DATA.plot.scaledensity
    DATA.clplot = 2;
end
C = GetClusterDef(DATA.cluster, DATA.currentcluster);

for j = 1:size(plots,1)
    mysubplot(2,4,j);
    t = 0;
    if ismember(DATA.plottype, [3 4 7])
        PlotPCs(DATA.TemplateScores,plots(j,1),plots(j,2),DATA.clplot,clid,DATA.colors,C,'fixrange',args{:});
       t =  title(sprintf('%s vs %s',...
            DATA.TemplateLabels{plots(j,1)},...
            DATA.TemplateLabels{plots(j,2)}));
        if isfield(E,'gmfit') && E.space(1) == 6 && E.space(2) ==4 
            [a,xi] = ismember(plots(j,1),DATA.tmplspace(1,:));
            [b,yi] = ismember(plots(j,2),DATA.tmplspace(1,:));
            if a && b && size(E.gmfit.mu,2) > max([xi yi])
                hold on;
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
    elseif ismember(DATA.plottype,[2 6 8 9])
        PlotVals(DATA,vpts(j,[1 2]),vpts(j,[3 4]),DATA.clplot,clid,DATA.colors,'fixrange');
        set(gca,'UserData',vpts(j,:));
        if DATA.usegmcid && sum(ismember([vpts(j,2) vpts(j,4)],DATA.vspace)) == 2
            k = find(DATA.vspace == vpts(j,2));
            m = find(DATA.vspace == vpts(j,4));
            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)
        end
    else
        if ismember(DATA.plottype,[11 12]) %tests
            PlotPCs(pcs,plots(j,1),plots(j,2),DATA.clplot,clid,DATA.colors, C,'fixrange',args{:});
        elseif ismember(DATA.plottype,[10]) %ICA
            PlotPCs(DATA.icas,DATA.pcplots(j,1),DATA.pcplots(j,2),DATA.clplot,clid,DATA.colors, C,'fixrange',args{:});
        else
            PlotPCs(DATA.pcs,DATA.pcplots(j,1),DATA.pcplots(j,2),DATA.clplot,clid,DATA.colors, C,'fixrange',args{:});
        end
        set(gca,'UserData',DATA.pcplots(j,:));
        if DATA.usegmcid && sum(ismember(DATA.pcplots(j,:),[1:4])) == 2 && isfield(DATA.cluster,'gmfit')
            k = find(DATA.pcspace == DATA.pcplots(j,1));
            m = find(DATA.pcspace == DATA.pcplots(j,2));
            plot(DATA.cluster.gmfit.mu(:,k),DATA.cluster.gmfit.mu(:,m),'c+','linewidth',2)
        end
        if j == 1
            text(0,1.2,sprintf('E%dP%d',DATA.exptno,p),'units','normalized','fontsize',DATA.gui.fontsize(1));
        end
    end
    if ~ismember(j,[1 4])
        set(gca,'ytick',[]);
    end
    if ~ismember(j,[5:8])
        set(gca,'xtick',[]);
    end
    if isempty(setids)
    axis('tight');          
    end
    xl = get(gca,'Xlim');
    yl = get(gca,'Ylim');
    if t && ishandle(t)
        set(t,'VerticalAlignment','top','position',[mean(xl) yl(2)]);
    end
    
    if DATA.showdipvals
    text(mean(xl),yl(2),sprintf('%.2f',dipvals(j)),'VerticalAlignment','top','color','k');
    end
    clusterplot(:,j) = GetClusterPlots(DATA,E,plots,j);
    %would like to used cluster, not E in future, so it doesn't depend on
    %which is current cluster
%    clusterplot(:,j) = GetClusterPlots(DATA,DATA.cluster,plots,j);
end
   if isfield(DATA.Expt.Header,'expname')
       exname = DATA.Expt.Header.expname;
   else
       exname = [];
   end
if DATA.plottype == 1
    p = DATA.probelist(DATA.probe);
    if isfield(DATA,'alldips')
        if size(DATA.alldips,2) > 1
            b = max(DATA.alldips');
        else
            b = DATA.alldips;
        end
        tstr = sprintf('%s %.2f, dt%.2f, csd%.2f',num2str(p),b(1).*100,b(2).*100,b(3).*100);
    elseif DATA.dvdt
        tstr = sprintf('E%dP%d: dvdt %s',DATA.exptno,p(1),exname);
    elseif DATA.csd
        tstr = sprintf('E%dP%d: csd %s',DATA.exptno,p(1),exname);
    else
        tstr = sprintf('E%dP%d %s',DATA.exptno,p(1),exname);
    end
else

   tstr = sprintf('E%dP%d%s %s',DATA.exptno,p(1),DATA.probelabel,exname);
end
mysubplot(2,4,1);
th = text(0.5,0.0,tstr,'units','normalized','VerticalAlignment','Top','fontsize',DATA.gui.fontsize(1));
if isfield(DATA,'cluster') && isfield(DATA.cluster,'auto') && DATA.cluster.auto == 0
    set(th,'color','r');
end

DATA.maintitle = th;
DATA.clustericon = SetClusterIcon(DATA);
    
for j = 1:size(clusterplot,2)
for k = 1:size(clusterplot,1)
if clusterplot(k,j)
    mysubplot(2,4,clusterplot(k,j)); %need to find right graph
    hold on;
    if k > 1
        E.next{k-1}.color = DATA.colors{k+1};
        DATA.elmousept.h = DrawEllipse(E.next{k-1});
    else
    E.color = DATA.colors{2};
    DATA.elmousept.h = DrawEllipse(E);
    end
    DATA.elmousept.handles(k) = DATA.elmousept.h;
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
end
end
if length(DATA.elmousept.handles) >= DATA.currentcluster;
    DATA.elmousept.h = DATA.elmousept.handles(DATA.currentcluster);
end

[iscell, cellid] =  isacell(DATA, DATA.exptno, ProbeNumber(DATA));
if iscell
    mysubplot(2,4,4);
    for j = 1:length(cellid)
        if cellid(j) > 0
        text(1,0.1*j,sprintf('Cell %d',cellid(j)),'units','normalized','color',DATA.colors{j+1},...
            'HorizontalAlignment','Right','fontsize',DATA.gui.fontsize(1));
        end
    end
end
if isfield(DATA,'energy')
SetFigure(DATA.tag.vare, DATA);
subplot(1,1,1);
PlotVarE(DATA);
end

function clusterplot = GetClusterPlots(DATA,E, plots,pt)
clusterplot(1) = GetClusterPlot(DATA,E, plots,pt);
j = 1;
if isfield(E,'next')
for j = 1:length(E.next)
    clusterplot(j+1) = GetClusterPlot(DATA,E.next{j},plots,pt);
end
end


function h = SetClusterIcon(DATA)
h = NaN;    
    if DATA.elmousept.shape == 0
        str = 'O';
    elseif DATA.elmousept.shape < 0
        str = '';
    else
        str = '/';
    end
    if ishandle(DATA.clustericon)
        delete(DATA.clustericon);
    end
    if ishandle(DATA.maintitle)
        mysubplot(2,4,1);
    x = get(DATA.maintitle,'extent');
    h = text(0,0,str,'units','Normalized','VerticalAlignment','Top','fontweight','bold','fontsize',DATA.gui.fontsize(1));
    set(h,'color',DATA.colors{DATA.currentcluster+1});
    else
        mysubplot(2,4,1);
        h = text(0,0,str,'units','Normalized','VerticalAlignment','Top','fontweight','bold','fontsize',DATA.gui.fontsize(1));
        set(h,'color',DATA.colors{DATA.currentcluster+1});
    end
    DATA.clustericon = h;


function clusterplot = GetClusterPlot(DATA,E,plots, pt)
    clusterplot = 0; 
    spaces = [1 2 3 3 5 6 7 8 9 10];
    if ~isfield(E,'pcplot') || isempty(E.pcplot) || E.shape == 2 || length(E.pcplot) < 2
        clusterplot = 0;
    elseif DATA.plottype == 2 && length(E.pcplot) == 4 && sum(DATA.vpts(pt,:) == E.pcplot) == 4
        clusterplot = pt;
    elseif ismember(DATA.plottype,[1 3 4]) && sum(plots(pt,:) == E.pcplot(1:size(plots,2))) == 2 ...
        && spaces(DATA.plottype) == E.space(1)
        clusterplot = pt;
    end
    if isfield(E,'space') && length(E.space) > 2 && E.space(1) == DATA.plottype
        if E.space(1) == 2 && sum(DATA.vsmps(plots(pt,:)) == E.space([3 5])) == 2
            if DATA.vpts(pt,3) == E.space(4) %check probe too
                clusterplot= pt;
            end
        elseif E.space(1) ~= 2 && sum(plots(pt,:) == E.space(2:3)) == 2
            clusterplot= pt;
        end
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
    hold off;
    DensityPlot(DATA.energy(1,DATA.nid),DATA.spkvar(c,:)./DATA.energy(1,:));
else
hold off;
if DATA.usegmcid && length(DATA.gmcid)
    cls = unique(DATA.gmcid);
    for j = length(cls)
        id = find(DATA.gmcid ==cls(j));
        plot(DATA.energy(1,id),DATA.spkvar(c,id)./DATA.energy(1,id),'.','markersize',1,'color',colors(j));
    end
elseif isfield(DATA,'nid')
plot(DATA.energy(1,DATA.nid),DATA.spkvar(c,DATA.nid)./DATA.energy(1,DATA.nid),'.','markersize',1);
hold on;
plot(DATA.energy(1,DATA.clid),DATA.spkvar(c,DATA.clid)./DATA.energy(1,DATA.clid),'r.','markersize',1);
else
plot(DATA.energy(1,DATA.uid),DATA.spkvar(c,DATA.uid)./DATA.energy(1,DATA.uid),'.','markersize',1);
end
end
title(sprintf('%d/%d Spikes',length(DATA.clid),DATA.nevents));


function DATA = RestrictTimeRange(DATA, t)
   DATA.uid = find(DATA.t >= t(1) & DATA.t <= t(2));
    DATA.restricttimerange = t;
    if isfield(DATA,'clid')
       DATA.clid = DATA.clid(ismember(DATA.clid,DATA.uid));
       DATA.nid = DATA.nid(ismember(DATA.nid,DATA.uid));
    end
    id = find([DATA.Expt.Trials.TrialStart] < t(1).*10000 | [DATA.Expt.Trials.TrialStart] > t(end).*10000);
    DATA.excludetrialids = [DATA.Expt.Trials(id).id];
    set(DATA.toplevel,'UserData',DATA);
    SetTrialList(DATA);

function DATA = ExcludeTrials(DATA, trials, add)

    if add == 1
        DATA.excludetrialids = cat(2, DATA.excludetrialids,[DATA.Expt.Trials(trials).id]);
    else
        DATA.excludetrialids = trials;
    end
        
    id = find(ismember([DATA.Expt.Trials.id],DATA.excludetrialids));
    iid = [];
    for j = id
        oid = find(DATA.t > DATA.Expt.Trials(j).Start(1)./10000 - DATA.preperiod & ...
            DATA.t < DATA.Expt.Trials(j).End(end)/10000 + DATA.postperiod);
        iid = [iid oid];
    end
    uid = unique(iid);
    DATA.clst(uid) = -1;
    DATA.uid = setdiff(1:DATA.nevents,uid);
    set(DATA.toplevel,'UserData',DATA);
    SetTrialList(DATA);

function SetCellCompare(a,b, cellid)
DATA = GetDataFromFig(a);
D = get(gcf,'UserData');
if ~CellIsEmpty(DATA.CellDetails,'MeanSpike', cellid)
    DATA.comparecell(cellid) = ~DATA.comparecell(cellid);
    if DATA.comparecell(cellid)
        AddCellMean(DATA,cellid);
    elseif strcmp(D.plottype,'QuickSpikes')s %redraw without
        QuickSpks(DATA,1000);
    end
end
set(DATA.toplevel,'UserData',DATA);
    
function AddCellMean(DATA, cellid)
    
    if isempty(cellid) || sum(cellid) ==0
        return;
    end
for c = 1:length(cellid)    
if ~CellIsEmpty(DATA.CellDetails,'MeanSpike', cellid(c))
    Ms = DATA.CellDetails.MeanSpike{cellid(c)};
    voff = DATA.voffset(DATA.chspk)-DATA.voffset(ProbeNumber(DATA));
    hold on;
    for j = 1:size(Ms.ms,1);
        plot(Ms.ms(j,:)+voff(j),'k-','linewidth',2);
    end
end 
end


function FinishSpikePlot(DATA)
    AddCellMean(DATA,find(DATA.comparecell));

function AllV = GetAllV(DATA)
    if DATA.interactive < 0 && isfield(DATA,'AllV')
        AllV = DATA.AllV;
    else
        AllV = getappdata(DATA.toplevel,'AllV');
        if isempty(AllV) && isfield(DATA,'AllV')
            AllV = DATA.AllV;
        end
    end
        
        
function ph = PlotSpikes(DATA,spkid, varargin)

dvdt = 1;
fixy = 0;
colors = DATA.colors;
scale = DATA.plotspk.muscale;
yl = [];
j = 1;
quicktest = 1;
showall = 0;
setcontext = 0;
tpt = find(DATA.spts ==0);
ph = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'fixy',3)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            yl = varargin{j};
        else
            fixy = 1;
        end
    elseif strncmpi(varargin{j},'setcontext',7)
        setcontext = 1;
    elseif strncmpi(varargin{j},'showall',7)
        showall = 1;
    end
    j = j+1;
end

if DATA.interactive < 0
    return;
end

AllV = GetAllV(DATA);
SetFigure(DATA.tag.spikes, DATA);
if fixy
    yl = get(gca,'ylim');
end
subplot(1,1,1);
if ~isfield(DATA,'clid')
    DATA.clid = [];
    DATA.nid = 1:DATA.nevents;
end

if DATA.plotspk.oneprobe
    chspk = DATA.probe(1);
else
    chspk = DATA.plotspk.probes;
end

if max(spkid) > DATA.nevents
    spkid = spkid(spkid<=DATA.nevents);
end
spkid = spkid(spkid > 0);
        if DATA.plotspk.includeprepost
            pre = DATA.preperiod;
            post = DATA.postperiod;
        else
            pre = 0;
            post = 0;
        end
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
    if DATA.plotspk.bytrial

        nt = DATA.currenttrial;
        T = DATA.Expt.Trials(nt);
        title(sprintf('Trial %d %.2f-%.2f: No Spikes  ed%.2f',T.Trial,T.Start(1)-pre,T.End(end)+post,T.ed));
    end
    return;
end



if isfield(DATA,'clst') && length(DATA.clst) >= max(spkid)   
    clst = DATA.clst;
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
elseif twoclusters
    ids = {[] []}; %min necessary
    for j = 1:length(nc)
    ids{nc(j)} = spkid(find(clst(spkid) == nc(j)));
    end
    if length(ids) > DATA.currentcluster
        id = ids{1+DATA.currentcluster};
    else
        id = ids{2};
    end
   nid = ids{1};
else
    id = find(ismember(spkid,DATA.clid));
    nid = find(ismember(spkid,DATA.nid));
    ids{2} = spkid(id);
    ids{1} = spkid(nid);
end
ispk = DATA.probe;
for j = 1:length(ids)
    V{j} = AllV(:,:,ids{j});
end

voff = DATA.voffset - DATA.voffset(ispk(1));
if DATA.plotcsd 
    AllCSD = GetCSD(DATA, DATA.plotcsd);
    for j = 1:length(ids)
        V{j} = AllCSD(:,:,ids{j});
    end
    if min(ispk) > 1
    ispk = ispk-1;
    end
    voff = (DATA.voffset - DATA.voffset(DATA.probe(1)));
elseif DATA.plotdvdt
for j = 1:length(ids)
    V{j} = diff(AllV(:,:,ids{j}),1,2);
end
end

if DATA.plotspk.submax
    for j = 1:length(V)
        V{j} = V{j} - repmat(max(V{j},[],2),[1 size(V{j},2) 1]);
    end
elseif DATA.plotspk.subtrigger
    for j = 1:length(V)
        V{j} = V{j} - repmat(V{j}(:,tpt,:),[1 size(V{j},2) 1]);
    end
elseif DATA.plotspk.submin
    for j = 1:length(V)
        V{j} = V{j} - repmat(min(V{j},[],2),[1 size(V{j},2) 1]);
    end
elseif DATA.plotspk.submean
    for j = 1:length(V)
        V{j} = V{j} - repmat(mean(V{j},[],2),[1 size(V{j},2) 1]);
    end
end

V{j} = V{j} .* scale;

if quicktest
    l = size(V{1},2);
    hold off;
    x = [1:l NaN];
for c= chspk;
    if c > 0 & c <= DATA.nprobes
        for j = 1:length(V)
            h(j) = 0;
        nV = squeeze(V{j}(c,:,:)) + voff(c);
        if length(nV) == 0
        elseif size(nV,1) > 1 
            h(j) = plot([0 30],[0 0],'color',colors{j});
            nV(l+1,:) = NaN;
            set(h(j),'Ydata',reshape(nV,1,prod(size(nV))),'Xdata',repmat(x,1,size(nV,2)));
            hold on;
        else
            h(j) = plot(1:l,nV,'color',colors{j});
            hold on;
            end
        end
        if setcontext && sum(h) > 0 && j > 0
            DATA = AddAxisContextMenu(DATA, gca);
            AddLineContextMenu(DATA, h);
        end
        if c == DATA.probe(1)
            ph = h;
        end
        h = text(size(V{1},2)-2,voff(c)+0.5,sprintf('%d',c),'fontsize',DATA.gui.fontsize(1));
        set(h,'fontweight','bold');
        if DATA.plotspk.showmean
            if isfield(DATA.cluster,'MeasSpike')
                plot(1:l,DATA.cluster.MeanSpike.ms(c,:)+voff(c),'-','linewidth',2,'color',colors{2});
            end
            for j = 1:length(DATA.cluster.next)
                if isfield(DATA.cluster.next{j},'MeanSpike')
                    plot(1:l,DATA.cluster.next{j}.MeanSpike.ms(c,:)+voff(c),'-','linewidth',2,'color',colors{j+2});
                end
            end
        end
    end

end
if DATA.plotrv
    rV = getappdata(DATA.toplevel,'AllrV');
    plot(rV(:,spkid),'k-');
end
else
hold off;


for c= [ispk(1)-1:1:ispk(1)+1];
    if c > 0 & c <= DATA.nprobes
        for j = 1:size(V)
            plot(squeeze(voff(c)+V{j}(c,:,:)),'color',colors{j});
            hold on;
        end
        h = text(size(V{1},2),voff(c),sprintf('%d',c),'fontsize',DATA.gui.fontsize(1));
    end
end
end
if length(yl) == 2
    set(gca,'ylim',yl);
end
hold off;

if DATA.plotspk.bytrial && isfield(DATA.Expt,'Trials');
    nt = max([DATA.currenttrial 0]);
    T = DATA.Expt.Trials(nt);
title(sprintf('Trial %d %.2f-%.2f: %d/%d ed%.2f',T.Trial, (T.Start(1)./10000)-pre,(T.End(end)./10000)+post, length(id),length(spkid),DATA.Expt.Trials(nt).ed));
else
title(sprintf('Spikes %d-%d(%.3f-%.3f): %d/%d',...
    spkid(1),spkid(end),DATA.t(spkid(1)),DATA.t(spkid(end)),length(id),length(spkid)));
end

if DATA.plotspk.showfullv && DATA.plotspk.bytrial
    tt = [T.Start(1) T.End(end)]./10000 +[-pre post];
    SetFigure('FullV', DATA);
    hold off;
    Vall = PlotFullV(DATA,tt);

    st = DATA.cluster.spts .* Vall.samper;
    VT.ids = ids;
    VT.toplevel = DATA.toplevel;
    VT.probe = chspk;
    VT.spkid = 0;
    VT.t = DATA.t;
    chspk = UseProbeList(DATA,DATA.plotspk.nfullvprobes);
    for j = 1:length(V)
        for k = 1:size(V{j},3)
            for p = 1:length(chspk)
                plot(st+DATA.t(ids{j}(k)),V{j}(chspk(p),:,k)+voff(chspk(p)),'color',DATA.colors{j});
            end
        end
    end

    id = find(Vall.t > tt(1)  & Vall.t < tt(2));
    if isfield(Vall,'intscale') && isinteger(Vall.V)
        vscale = Vall.intscale(1)./Vall.intscale(2);
    else
        vscale = 1;
    end
    for p = 1:length(chspk)
        plot(Vall.t(id),double(Vall.V(chspk(p),id)).*vscale+voff(chspk(p)),'k');
        hold on;
    end
    st = DATA.cluster.spts .* Vall.samper;
    for j = 1:length(V)
        for k = 1:size(V{j},3)
            for p = 1:length(chspk)
                plot(st+DATA.t(ids{j}(k)),V{j}(chspk(p),:,k)+voff(chspk(p)),'color',DATA.colors{j});
            end
        end
    end
    axis('tight');
    if isfield(DATA,'onespiketime')
        for j = 1:length(DATA.onespiketime)
            plot([DATA.onespiketime./10000 DATA.onespiketime./10000],get(gca,'ylim'),'b:');
        end
    end
    set(gcf,'UserData',VT,'KeyPressFcn',{@FullVKeyPressed});

end

function ShowFullV(src,b, fcn)
   DATA = GetDataFromFig(src);
   pos = get(gca,'currentpoint');
   x = get(src,'Xdata');
   [a,b] = min(abs(x-pos(1,1)));
   t(1) = x(b)-0.001;
   t(2) = x(b)+0.001;

   SetFigure(DATA.tag.fullv,DATA);
   [V,id] = PlotFullV(DATA,t);
   [V,id] = PlotFullV(DATA,t,'addmean','hold','color','g');
   T = getappdata(DATA.toplevel,'TriggerV');
   if ~isempty(T)
       plot(V.t(id),T(id));
   end
   c = find(V.t(id) < x(b));
   c = c(end) + DATA.cluster.spts(1)-1;
   if DATA.cluster.trigdt == 0
       plot(V.t(id([1 end])),[DATA.cluster.Trigger DATA.cluster.Trigger]);
   end
   if isfield(DATA,'oldSpikes')
       [a,b] = min(abs(DATA.oldSpikes.times-x(b)));
       hold on;
       plot(V.t(id(c+[1:size(DATA.oldSpikes.values,2)])),DATA.oldSpikes.values(b,:),'r');
       if DATA.cluster.tsmooth
           smv = smooth(DATA.oldSpikes.values(b,:),DATA.cluster.tsmooth,'gauss');
           plot(V.t(id(c+[1:size(DATA.oldSpikes.values,2)])),smv,'m');
       end
       axis('tight');
   end
   
function [Vall, id] = PlotFullV(DATA, t, varargin)

    addmean = 0;
    dohold = 0;
    color = 'k';
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'addmean',7)
            addmean = 1;
        elseif strncmpi(varargin{j},'color',5)
            j = j+1;
            color = varargin{j};
        elseif strncmpi(varargin{j},'hold',4)
            dohold = 1;
        end
        j = j+1;
    end
    Vall = getappdata(DATA.toplevel,'Vall');
    if size(Vall.V,1) > 1
        voff = DATA.voffset- DATA.voffset(ProbeNumber(DATA));
    else
        voff = 0;
    end
    id = find(Vall.t > t(1)  & Vall.t < t(2));
    chspk = UseProbeList(DATA,DATA.plotspk.nfullvprobes);
    if isfield(Vall,'intscale') && isinteger(Vall.V)
        vscale = Vall.intscale(1)./Vall.intscale(2);
    else
        vscale = 1;
    end
    if dohold
        hold on;
    else
        hold off;
    end
    for p = 1:length(chspk)
        V = double(Vall.V(chspk(p),id)).*vscale+voff(chspk(p));
        if addmean && isfield(Vall,'meanV')
            V = V + Vall.meanV(id);
        end
        plot(Vall.t(id),V,'color',color);
        hold on;
    end
    axis('tight');
        
function OldSetMenuCheck(F, tag, value)

    
onoff = {'off' 'on'};
if ischar(F) %tag for a figure  
    F = findobj('type','figure','tag',F);
    if isempty(F);
        return;
    end
end

it = findobj(F,'Tag',tag);
if length(it) == 1
    set(it,'Checked',onoff{1+value});
end

function SetGUI(DATA)

    if DATA.interactive < 0
        return;
    end
    onoff = {'off' 'on'};
    c = get(DATA.toplevel,'Children');
    SetMenuCheck(DATA.tag.spikes,'PlotByTrial',DATA.plotspk.bytrial);
    SetMenuCheck(DATA.toplevel,'PlotISI',DATA.plot.isi);
    SetMenuCheck(DATA.toplevel,'PlotExpt',DATA.plot.expt);
    SetMenuCheck(DATA.toplevel,'advanceprobe',DATA.auto.advanceprobe);
    SetMenuCheck(DATA.toplevel,'checkxcorr',DATA.auto.checkxcorr);
    SetMenuCheck(DATA.toplevel,'csd',DATA.csd > 0);
    SetMenuCheck(DATA.toplevel,'dvdy',DATA.csd == 1);
    SetMenuCheck(DATA.toplevel,'dvdt',DATA.dvdt);
    SetMenuCheck(DATA.toplevel,'ProbeModeMenu',DATA.probeswitchmode,'exclusive');

function PlotTemplateScores(DATA, TemplateScores, probes)

p = 3;
if length(probes) > 1
p = 2;
slope = mean(squeeze(TemplateScores(p,1,:)))./mean(squeeze(TemplateScores(p,2,:)));
diffs = squeeze(TemplateScores(p,1,:) - slope.* TemplateScores(p,2,:));
hold off;
p = 1;
plot(squeeze(TemplateScores(p,1,DATA.nid)),diffs(DATA.nid),'.','markersize',1)
hold on;
plot(squeeze(TemplateScores(p,1,DATA.clid)),diffs(DATA.clid),'r.','markersize',1)
hold off;
end



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


    function  sgn = CheckSign(C, x, energy)
        e(1) = mean(energy(find(x > C.crit)));
        e(2) = mean(energy(find(x < C.crit)));
        sgn = 0;
    if e(2) > 2 * e(1)
        sgn = -1;
    elseif e(1) > 2 * e(2)
        sgn = 1;
    end
    if sgn == 0
        if mean(x) < 0 && skewness(x) < 0
            sgn = -1;
        else
            sgn = 1;
        end
    end
        
    
function [dip, details] = oldFindDip(values, energy, varargin)
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
        j = j+1;
        tag = varargin{j};
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
    SetFigure(tag);
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
    SetFigure(f);
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
subdir = [];
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
        elseif strcmp(varargin{j},'subdir')
            j = j+1;
            subdir = varargin{j};;
        elseif strcmp(varargin{j},'auto')
            filetype = 1;
        elseif strcmp(varargin{j},'DistanceMatrix')
            filetype = 4;
        elseif strcmp(varargin{j},'log')
            filetype = 2;
        elseif strcmp(varargin{j},'clusterlog')
            filetype = 3;
        end
        j = j+1;
    end
    if isdir(name)
        a = name;
    else
        [a,b] = fileparts(name);
    end
    if length(subdir)
        a = [a '/' subdir];
    end
    if filetype == 2
        cname = [a '/' prefix 'AutoLog.txt'];
    elseif filetype == 3
        cname = [a '/ClusterLog.txt'];
    elseif filetype == 4
        cname = [a '/' prefix 'DistanceMatrix.mat'];
    elseif filetype == 1
        cname = [a '/' prefix 'AutoClusterTimes.mat'];
    else
        cname = [a '/' prefix 'ClusterTimes.mat'];
    end

        
    
function HistButtonPressed(src, data)
DATA = get(src,'UserData'); %not the main fig data. hust for this one

start = get(gca,'CurrentPoint');
bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'})

if bt == 2
    mode = 1;  %allow angled line
else
    mode = 3; %just translate criterion
end
%mode = 3;  %only allow line translation at the moment
DATA.elmousept.mode = mode;
DATA.elmousept.shape = 1;
DATA.elmousept.pos(1) = start(1,1);
yl = get(gca,'ylim');
if mode == 3
DATA.elmousept.pos(3) = start(1,1);
DATA.elmousept.pos(4) = yl(2);
DATA.elmousept.pos(2) = yl(1);
else
DATA.elmousept.pos(2) = start(1,2);
DATA.elmousept.pos(3) = start(1,1);
DATA.elmousept.pos(4) = start(1,2);
end
if start(1,2) < yl(2) && start(1,2) > yl(1)
    DATA.elmousept.down = 1;
else
    return;
end
DATA.elmousept.axis = gca;
DATA.elmousept.plotargs = {};
DATA.elmousept.h = [];
DATA.elmousept.h = DrawEllipse(DATA.elmousept);
subplot(2,1,1);
hold on;
DATA.elmousept.h2 = plot(DATA.elmousept.pos([1 3]),get(gca,'ylim'),'r-');
subplot(2,1,2);
set(src,'UserData',DATA);

function HistButtonReleased(src, data)
HDATA = get(src,'UserData');  % Just for this figure
start = get(gca,'CurrentPoint');
F = gcf;
if isfield(HDATA,'elmousept') && HDATA.elmousept.down
HDATA.elmousept.down = 0;
HDATA.elmousept.done = 1;
p = HDATA.elmousept.pos;
if HDATA.elmousept.mode == 3
HDATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) diff(p([1 3]))/2 diff(p([2 4]))/2]; 
else
HDATA.elmousept.xyr = [mean(p([1 3])) mean(p([2 4])) diff(p([1 3]))/2 diff(p([2 4]))/2]; 
end
set(src,'UserData',HDATA);

%now get main data 
DATA = GetDataFromFig(src);
mode = HDATA.elmousept.mode;
if HDATA.elmousept.mode == 3
DATA.cluster.crit = p(1);
E = BoundaryFromCluster(HDATA.elmousept, DATA.cluster, DATA.currentcluster);
else
E = BoundaryFromCluster(HDATA.elmousept, DATA.cluster, DATA.currentcluster);
E.pos = p;
end
[cl, DATA.cluster, DATA.xy{DATA.currentcluster}] = ClassifySpikes(DATA,E,DATA.quickcutmode);
DATA.cluster.auto = 0;
DATA.clid = cl.id;
DATA.nid = cl.nid;
DATA.clst = cl.clst;


if isfield(cl,'MeanSpike')
    DATA.MeanSpike = cl.MeanSpike;
end
set(DATA.toplevel,'UserData',DATA);
E.h = [];
if DATA.watchplots == 0 %otherwise plotted in Replot
    DATA = ReplotPCs(DATA,E);
end
if mode == 1
    PlotHistogram(DATA,[]);
else
    figure(F);
    title(sprintf(' G2 %.2f',DATA.cluster.fitdprime(1)));
end

end

function C = ClusterInfo(DATA)

if DATA.currentcluster > 1
    if length(DATA.cluster.next) >= DATA.currentcluster-1
        C = DATA.cluster.next{DATA.currentcluster-1};
        if isempty(C) %can happen after deleting
            C.space = [0 0];
        end
    else
        C.space = [0 0];
    end
else 
    C = DATA.cluster;
end


function HistButtonDragged(src, data)
DATA = get(src,'UserData');

if isfield(DATA,'elmousept') && DATA.elmousept.down > 0
    if DATA.elmousept.mode == 3
        start = get(gca,'CurrentPoint');
        DATA.elmousept.pos(3) = start(1,1);
        DATA.elmousept.pos(1) = start(1,1);
        DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
        if isfield(DATA.elmousept,'h2')
            set(DATA.elmousept.h2, 'Xdata', DATA.elmousept.pos([1 3]));
        end
    else
        start = get(gca,'CurrentPoint');
        DATA.elmousept.pos(3) = start(1,1);
        DATA.elmousept.pos(4) = start(1,2);
        DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});

    end
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
if isfield(E,'pos')
a = (E.pos(3)-E.pos(1))/2; %radius 
b = (E.pos(4)-E.pos(2))/2;
cntr(1) = mean(E.pos([1 3]));
cntr(2) = mean(E.pos([2 4]));
elseif isfield(E,'xyr')
    a = E.xyr(3);
    b = E.xyr(4);
    cntr(1) = E.xyr(1);
    cntr(2) = E.xyr(2);
end
    
sn = 0;
cn = 1;

x = linspace(0,a);
sn = sin(E.angle);
cn = cos(E.angle);
if isfield(E,'aspectratio') && E.aspectratio > 0
    b  = b ./E. aspectratio;
    y =  sqrt(b.^2 - (x.*b/a).^2);
else
    y =  sqrt(b.^2 - (x.*b/a).^2);
end

x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+cntr(1);
if isfield(E,'aspectratio') && E.aspectratio > 0
    y = yr.*E.aspectratio+cntr(2);
else
    y = yr+cntr(2);
end

if isfield(E,'h') & ishandle(E.h)  & get(E.h,'parent') == gca
%    fprintf('Using handle %f\n',get(E.h,'parent'));
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
%fprintf('%.3f,%.3f -> %.3f,%.3f\n',x(1),y(1),x(2),y(2));
if isfield(E,'h') & ishandle(E.h) & get(E.h,'parent') == E.axis; 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end

function [distance, cluster] = FindNearestCluster(DATA, pos)
    distance = NaN;
    cluster = 0;
    
    if DATA.cluster.space(1) == DATA.plottype && sum(DATA.elmousept.pcplot == DATA.cluster.space(2:end)) == length(DATA.elmousept.pcplot)
        d(1) = DistanceToEllipse(DATA.cluster, pos);
    else
        d(1) = NaN;
    end
    for j = 1:length(DATA.cluster.next)
        if ~isempty(DATA.cluster.next{j}) && DATA.cluster.next{j}.space(1) == DATA.plottype && sum(DATA.elmousept.pcplot == DATA.cluster.next{j}.space(2:end)) == length(DATA.elmousept.pcplot)
            d(j+1) = DistanceToEllipse(DATA.cluster.next{j}, pos);
        else
            d(j+1) = NaN;
        end
    end
    [distance, cluster] = min(d);
    


function ButtonPressed(src, data)
DATA = GetDataFromFig(src);

DATA.ts = now;
start = get(gca,'CurrentPoint');
if InGraph(start,gca)
    mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
%    DATA.elmousept.mode = mode;
%only check for in ellipse if in correct space
    C = ClusterInfo(DATA);
    DATA.elmousept.pcplot = get(gca,'UserData');
    distance = DistanceToEllipse(DATA.elmousept,start(1,1:2));
%Not so simple. When an ellipse is drawn in a new space, C.space is not
%updated until clicking inside. 
%    if length(DATA.elmousept.pcplot) == 2 && length(C.space) > 1
%        if isfield(C,'space') &&  sum(C.space(end-1:end) == DATA.elmousept.pcplot) < 2
%            distance = 100
%        end
    %end
    if ishandle(DATA.elmousept.h)
        if get(DATA.elmousept.h,'Parent') == gca  %only inside if its teh right subplot
            axisok = 1;
        else
            axisok = 0;
        end
    else
        axisok = 1;
    end
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    set(gca,'xlimmode','manual','ylimmode','manual'); %dont rescale for ellispes
    DATA.elmousept.aspectratio = diff(yl)./diff(xl);
    if  mode  == 2 % R button press. Inside and existing ellispe sets cutting to that space
        DATA.elmousept.down = 3;
        DATA.elmousept.done = 0;
        DATA.elmousept.start = start(1,1:2);
        [x, c] =  FindNearestCluster(DATA, start(1,1:2));
        if x < 1
            DATA = SetEllipseDrawing(DATA,0,'cluster',c);
        end
    elseif distance < 1.05 && axisok%test fails for NaN
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
        DATA.elmousept.pos =[start(1,1) start(1,2) 0 0 ];
        DATA.elmousept.axis = gca;
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        DATA.elmousept.aspectratio = diff(yl)./diff(xl);
  %if drawing in a new space, set cluster angle to zero.    
        DATA.elmousept.pcplot = get(gca,'UserData');
        if length(DATA.elmousept.pcplot) == 2 && length(C.space) > 1
            if isfield(C,'space') &&  sum(C.space(end-1:end) == DATA.elmousept.pcplot) < 2
                DATA.elmousept.angle = 0;
            end
        end
        
        if 0
        if ~isfield(DATA.elmousept,'angle')
            DATA.elmouept.angle = 0 ;
        end

        if ~isfield(DATA.elmousept,'plotargs')
            DATA.elmouept.plotargs = {} ;
        end
        end
    end
set(DATA.toplevel,'UserData',DATA);
end

 function len = LineLength(l)
     
     if length(l) == 4
         len = abs(diff(l([1 3])) + i * diff(l([2 4])));
     end

function distance = DistanceToEllipse(E, pos);
   
%When called from a new graph, ther is nothing to check, but E may have
%leftover bits in
if isempty(E) | ~isfield(E,'pos') | ~isfield(E,'angle');
    distance = NaN;
    return;
end

if E.shape == 1
r(1) = LineLength(E.pos)/3;
r(2) = r(1)./10;
a(1) = (E.pos(3)+E.pos(1))/2; %x radius
a(2) = (E.pos(4)+E.pos(2))/2;

xy = pos - a;
xy = xy ./r;
cn = cos(-E.angle);
sn = sin(-E.angle);
p(1) = xy(1) * cn + xy(2) * sn;
p(2) = xy(2) * cn - xy(1) * sn;

distance = sum(p.^2);
else
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
end

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
if length(DATA.elmousept.pcplot) > 2 && isnan(DATA.elmousept.pcplot(1))
    DATA.elmousept.space = DATA.elmousept.pcplot(2:end);
else
    DATA.elmousept.space = [DATA.plottype DATA.elmousept.pcplot];
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
if mode == 2  || DATA.elmousept.shape  ==1 
PCCluster(DATA,DATA.elmousept,1);
end

function ScrollWheel(src, evnt)
DATA = GetDataFromFig(src);
RotateCluster(DATA, 0.02*evnt.VerticalScrollCount);


function DATA = RotateCluster(DATA, angle)
%just in case its called before cluster is set
if isfield(DATA,'elmousept') && isfield(DATA.elmousept,'angle')
DATA.elmousept.angle = DATA.elmousept.angle+angle;
fprintf('Angle %.2f AR %.3f\n',DATA.elmousept.angle,DATA.elmousept.aspectratio);
DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
set(DATA.toplevel,'UserData',DATA);
end

 function PCButtonDragged(src, data)
DATA = GetDataFromFig(src);
if isfield(DATA,'elmousept')
if  DATA.elmousept.down == 1
    start = get(gca,'CurrentPoint');
    DATA.elmousept.pos(3) = start(1,1);
    DATA.elmousept.pos(4) = start(1,2);
    DATA.elmousept.h = DrawEllipse(DATA.elmousept, DATA.elmousept.plotargs{:});
%    drawnow;
%    mytoc(DATA.ts);
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
set(DATA.toplevel,'UserData',DATA);
end


function pos = xyr2pos(xyr)

    pos(1) = xyr(1)-xyr(3);
    pos(2) = xyr(2)-xyr(4);
    pos(3) = xyr(1)+xyr(3);
    pos(4) = xyr(2)+xyr(4);

function C = GetClusterDef(cluster, n)
%return relevant member of teh clsuter struc    
C = [];
    if n == 1
        C = cluster;
    elseif n-1 <= length(cluster.next)
        C = cluster.next{n-1};
    else 
        C = cluster;
    end
        
function DATA = SetEllipseDrawing(DATA, shape,varargin)

cluster = 1;
boundarytype=1;
plotargs = {};
F = DATA.toplevel;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cluster',4)
            j = j+1;
            cluster = varargin{j};
        elseif strncmpi(varargin{j},'boundarytype',10)
            j = j+1;
            boundarytype = varargin{j};
        elseif strncmpi(varargin{j},'figure',6)
            j = j+1;
            F = findobj('type','figure','tag',varargin{j});
        else
            plotargs = {plotargs{:} varargin{j}};
        end
        j = j+1;
    end

    tic;
DATA.elmousept.boundarytype = boundarytype;
DATA.currentcluster = cluster;
DATA.elmousept.h= -1;
if length(DATA.elmousept.handles) >= cluster && ishandle(DATA.elmousept.handles(cluster))
    DATA.elmoustpt.h = DATA.elmousept.handles(cluster);
end
DATA.elmousept.shape= shape;
DATA.elmousept.down = 0;
DATA.elmousept.done = 0;
DATA.elmousept.angle = 0;
DATA.elmousept.cluster = cluster;
DATA.elmousept.plotargs = plotargs;
DATA.elmousept.color = DATA.colors{cluster+1};
C = GetClusterDef(DATA.cluster, cluster);
if isfield(C,'shape') && C.shape == 0
    DATA.elmousept.pos = xyr2pos(C.xyr);
    DATA.elmousept.angle = C.angle;
    DATA.elmousept.xyr = C.xyr;
end
DATA.elmousept.dragfcn = get(F,'WindowButtonMotionFcn');
DATA.clustericon = SetClusterIcon(DATA);

%should get old Fcns, then reset them after button release
set(F, 'WindowButtonDownFcn',@ButtonPressed);
set(F, 'WindowButtonMotionFcn',@PCButtonDragged);
set(F, 'WindowButtonUpFcn',@ButtonReleased);
set(F, 'WindowScrollWheelFcn',@ScrollWheel);
set(F,'UserData',DATA);
toc

sv = findobj('type','figure','Tag',DATA.tag.vare);
X.parentfig = DATA.toplevel;
X.elmousept = DATA.elmousept;
set(sv, 'WindowButtonDownFcn',@ButtonPressed);
set(sv, 'WindowButtonMotionFcn',@PCButtonDragged);
set(sv, 'WindowButtonUpFcn',@ButtonReleased);
set(sv, 'WindowScrollWheelFcn',@ScrollWheel);
set(sv,'UserData',X);

sv = findobj('type','figure','Tag',DATA.tag.tmplscore);
X.parentfig = DATA.toplevel;
X.elmousept = DATA.elmousept;
set(sv, 'WindowButtonDownFcn',@ButtonPressed);
set(sv, 'WindowButtonMotionFcn',@PCButtonDragged);
set(sv, 'WindowButtonUpFcn',@ButtonReleased);
set(sv, 'WindowScrollWheelFcn',@ScrollWheel);
set(sv,'UserData',X);

function MiscMenu(a, b, type)
    DATA = GetDataFromFig(a);
    onoff = {'off' 'on'};
    if sum(strcmp(type,{'savelayout' 'savedefaultlayout'}))
        f = fields(DATA.tag);
        for j = 1:length(f);
            it = findobj('type','figure','Tag', DATA.tag.(f{j}));
            if length(it) == 1
                Figpos.(DATA.tag.(f{j})) = get(it,'Position');
            end
        end
        if isempty(DATA.layoutfile) || strcmp(type,'savedefaultlayout') 
            DATA.layoutfile = DATA.defaultlayout;
        end
        [outname, pathname] = uiputfile(DATA.layoutfile);
        if outname
            DATA.layoutfile = [pathname outname];
            save(DATA.layoutfile,'Figpos');
            fprintf('Layout saved to %s\n',DATA.layoutfile);
        end
    elseif strcmp(type,'loadlayout')
        [afile, pathname] = uigetfile(DATA.layoutfile);
        if ischar(afile)
        DATA.layoutfile = [pathname afile];
        ApplyLayout(DATA);
        end
    elseif sum(strcmp(type,{'saveconfig' 'savedefaultconfig'}))
        if strcmp(type,'savedefaultconfig') || ~isfield(DATA,'configfile')
            DATA.configfile = DATA.defaultconfig;
        end
        [outname, pathname] = uiputfile(DATA.configfile);
        if outname
            DATA.configfile = [pathname outname];
            savefields = {'quickcutmode' 'probeswitchmode' 'auto' 'ptsz' 'plot' 'plotspk' 'gui'};
            DATA.plotspk.probes = [];
            SaveConfig(DATA, DATA.configfile,savefields,'verbose');
        end
    elseif strcmp(type,'scaledensity')
        DATA.plot.scaledensity = ~DATA.plot.scaledensity;
        set(a,'checked',onoff{1+DATA.plot.scaledensity});
        ReplotPCs(DATA,[]);
        set(DATA.toplevel,'UserData',DATA);
    elseif strcmp(type,'loadconfig')
        [afile, pathname] = uigetfile(DATA.configfile);
        if ischar(afile)
        DATA.configfile = [pathname afile];
        DATA = ReadConfig(DATA, DATA.configfile);
        end
    elseif strcmp(type,'tofront')
            FiguresToFront(DATA.tag);
    elseif strcmp(type,'plotcelllist')
        PlotCellList(DATA, 'showfig','reload');
    end
figure(DATA.toplevel);

    
 function SetOption(a, b, type)
     onoff = {'off' 'on'};
     DATA = GetDataFromFig(a);
     if nargin < 3
         type = get(a,'tag');
     end
     DATA.(type) = ~DATA.(type);
     set(a,'Checked',onoff{DATA.(type)+1});
     if DATA.dvdy
         DATA.csd = 1;
     elseif DATA.csd
         DATA.csd = 2;
     end
     if strcmp(type,'usegmcid')
         PlotHistogram(DATA,[]);
     elseif sum(strcmp(type,{'csd' 'dvdy' 'dvdt'}))
         DATA = SetPCs(DATA, 1);
     end
     set(DATA.toplevel,'UserData',DATA);
     
function PlotResult(a, b, type)
     onoff = {'off' 'on'};
     DATA = GetDataFromFig(a);
     if strcmp(type,'plotexpt')
         DATA.Expt = PlotExptCounts(DATA);
         plotexpt = 1;
     else
         DATA.plot.expt = ~DATA.plot.expt;
         set(a,'Checked',onoff{DATA.plot.expt+1});
         plotexpt = DATA.plot.expt;
     end
     if plotexpt
         DATA.Expt = PlotExptCounts(DATA);
     end
     set(DATA.toplevel,'UserData',DATA);
     
     
function HitTrial(data,b, cell)
    DATA = GetDataFromFig(data);
    Expt = DATA.Expt;
    f = gcf;
    D = get(f,'UserData');
    pos = get(gca,'currentpoint');
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    x = get(data,'xdata');
    y = get(data,'ydata');
    [a,id] = min(abs(pos(1,1)-x));
    
    PlotOneTrial(DATA, id);

     
function [Expt, plotres] = PlotExptCounts(DATA)
    Expt = DATA.Expt;
    fitres = [];
     SetFigure(DATA.tag.Expt, DATA);
     spkt = DATA.t(DATA.uid(DATA.clst == DATA.currentcluster+1)) .*10000;
     preperiod = DATA.preperiod .* 10000;
     postperiod = DATA.postperiod .* 10000;
     for j = 1:length(DATA.Expt.Trials)
         id = find(spkt > Expt.Trials(j).Start(1)-preperiod & spkt < Expt.Trials(j).End(end)+postperiod);
         spks = round(spkt(id)'-Expt.Trials(j).Start(1));
         Expt.Trials(j).Spikes = spks;
         Expt.Trials(j).count = sum(spks > 500 & spks < Expt.Trials(j).dur);
     end
     if strncmpi(DATA.plot.expttype,'trialcounts',10)
         PlotRateSequence(Expt,'callback',@HitTrial);
        %plotres =  PlotExpt(Expt,'seqt');
     else
         plotres = PlotExpt(Expt,'shown','fbox','rcnmin',10);
     end
     if DATA.plot.exptfit
        fitres = FitExpt(plotres,'plotfit');
     end
     t = get(gca,'title');
     if isfield(fitres,'fitstr')
         s = fitres(1).fitstr;
     else
         s = get(t,'string');
     end
     
     
     p = DATA.probelist(DATA.probe(1));
     set(t,'string',sprintf('P%dCl%d %s',p,DATA.currentcluster, s));
     
 function PlotMenu(a, b, type)
 
     DATA = GetDataFromFig(a);
     onoff = {'off' 'on'};
     if sum(strcmp(type,{'xyseq' 'xcorrprobes'}))
         DATA.plot.(type) = ~DATA.plot.(type);
         set(a,'Checked',onoff{DATA.plot.(type)+1});
         if DATA.plot.(type)
             if strcmp(type,'xyseq')
                 PlotXYSequence(DATA,DATA.cluster);
             elseif strcmp(type,'xcorrprobes')
                 CalcXcorr(DATA,[],'probes');
             end
         end
     elseif strcmp(type,'rateseq')
         if strcmp(DATA.plot.expttype,'trialcounts')
             DATA.plot.expttype = 'means';
             DATA.plot.expt = 0;
             set(a,'checked', 'off');
         else
             DATA.plot.expttype = 'trialcounts'
             DATA.plot.expt = 1;
             set(a,'checked', 'on');
             PlotExptCounts(DATA);
         end
         
     elseif strcmp(type,'xcorrprobe')
         DATA.plot.xcorrprobe = ~DATA.plot.xcorrprobe;
         set(a,'Checked',onoff{DATA.plot.xyseq+1});
         if DATA.plot.xyseq
             PlotXYSequence(DATA,DATA.cluster);
         end
     end
     set(DATA.toplevel,'UserData',DATA);
     
 function sdx = PlotXYSequence(DATA, probe, varargin)
    plottype = 1;
    clst = [];
    if length(varargin) && isnumeric(varargin{1})
        clst = varargin{1};
        j = 1;
    else
    j = 1;
    end
    while j <= length(varargin)
        if strncmpi(varargin{j},'noplot',6)
            plottype = 0;
        end
        j = j+1;
    end
    
    Clusters = getappdata(DATA.toplevel,'Clusters');
    if isstruct(probe)
        C = probe;
    elseif probe == 0 
        C = DATA.cluster;
    else
        DataClusters = getappdata(DATA.toplevel,'Clusters');
        C = DataClusters{probe};
    end
    if length(clst)
        C.clst = clst;
    end
    if ~isfield(C,'clst')
        C.clst = DATA.clst;
    end
    if ~isfield(C,'xy')
        C.xy = DATA.xy{1};
    end
    nc = unique(C.clst(:))';
    nt = length(DATA.Expt.Trials);
    T = DATA.Expt.Trials;
    if plottype > 0
    F = SetFigure(DATA.tag.oldxy, DATA);
    subplot(1,1,1);
    set(F,'name','XY sequence');
    hold off;
    end
%    need to use time or trial bins for this test, as rate changes 

    for j = nc
        if isfield(C,'excludetrialids')
         xcl = C.excludetrialids;
        else
            xcl = [];
        end
        id = find(C.clst == j);
        if plottype == 1
            plot(DATA.t(id),C.xy(id,1),'.','color',DATA.colors{j});
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
                id  = find(C.clst == j & DATA.t.*10000 > T(k).Start(1) & DATA.t.*10000 < T(k+smw).End(end));
            else
                id  = find(C.clst' == j & DATA.t.*10000 > T(k).Start(1) & DATA.t.*10000 < T(k+smw).End(end));
            end
            sds(k) = std(C.xy(id,1));
            sms(k) = mean(C.xy(id,1));
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
    if yl(2) > 0
        plot(ts(sid),sds .*yl(2)./max(sds),'ms');
    else
        plot(ts(sid),sds .*yl(1)./max(sds),'ms');
    end
    plot(ts(sid),sms,'r-');
    if C.shape ~= 0
        plot([DATA.t(1) DATA.t(end)],[C.crit C.crit],'r-');
    end
%    MarkTrialStarts(DATA.Expt,0,xcl);
    title(sprintf('SDindex %.2f (%.2f/%.2f) CV %.2f',sdx(end),ssds(end),msds(end),mdx(end)));
end


 function PlotISI(a, b, type)
 
     DATA = GetDataFromFig(a);
     if isfield(DATA,'clst')
         tid = find(DATA.clst == 2);
         t = DATA.t(tid);
     else
         DataClusters = getappdata(DATA.toplevel,'Clusters');
         t = DataClusters{DATA.probe(1)}.times;
     end
     isi = diff(t);
     SetFigure('ISI', DATA);
     set(gcf,'UserData',DATA.toplevel);
     if type == 1
         if DATA.plot.isi == 1
             DATA.plot.isi = 0;
             set(a,'Checked','off');
             set(DATA.toplevel,'UserData',DATA);
             return;
         else
             DATA.plot.isi = 1;
             set(a,'Checked','on');
             set(DATA.toplevel,'UserData',DATA);
         end
     elseif type == 0
         type = 1;
     end
     if type == 1
     [a,b] = hist(isi,[0:0.0005:0.1]);
     hold off;
     bar(b(1:end-1),a(1:end-1));
     id = find(isi < 0.002);
     if size(tid,1) == 1
         sid = cat(2,tid(id), tid(id+1));
     else
         sid = cat(1,tid(id), tid(id+1));
     end
%     sid = tid(id+1);
     ReplotPCs(DATA,[],'setid',sid);
     elseif type == 2
         hold off;
         for j = 1:length(isi)
         plot(t(j+1),isi(j),'o','buttondownfcn',{@HitISI, t(j+1)});
         hold on;
         end
         set(gca,'ylim',[0 0.005]);
     end
         
function HitISI(a,b, t)
DATA = GetDATAFromFig(a);
sid(2) = find(DATA.t == t);
fprintf('t=%.2f event %d\n',t,sid(2));
if size(DATA.t,2) == size(DATA.clst,2)
id = find(DATA.clst == DATA.clst(sid(2)) & DATA.t < t);
else
id = find(DATA.clst == DATA.clst(sid(2)) & DATA.t' < t);
end
sid(1) = id(end);

if length(sid)
PlotSpikes(DATA,sid);
ReplotPCs(DATA,[],'setid',sid);
SetFigure(DATA.tag.fullv, DATA);
PlotFullV(DATA,[DATA.t(sid(1))-0.0005 DATA.t(sid(2))+0.0005]);

end


function [isis, trials, spkids] = CalcISI(Trials, varargin)
%
%[isis, trials, spkids] = CalcISI(Trials, varargin)
%spkids are actually times
latency = 500;

isis = [];
trials = [];
spkids = [];
for j = 1:length(Trials)
    duration = Trials(j).End(end) - Trials(j).Start(1);
 spks = find(Trials(j).Spikes > latency & ...
     Trials(j).Spikes < duration+latency);
 isis = [isis diff(Trials(j).Spikes(spks)')];
 trials = [trials ones(size(spks(2:end)')) * j];
 spkids = [spkids Trials(j).Spikes(spks(2:end))'+Trials(j).Start(1)];
end


 function [C, fits] = OptimizeEllipse(DATA)
     c = DATA.currentcluster;
     state.cluster = c;
     if DATA.currentcluster > 1
         guess(1:4) = DATA.cluster.next{c-1}.xyr;
         guess(5) = DATA.cluster.next{c-1}.angle;
%         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
         if isfield(DATA.cluster.next{c-1},'aspectratio')
             state.aspectratio = DATA.cluster.next{c-1}.aspectratio;
         else
             state.aspectratio = 1;
         end
         C = DATA.cluster.next{c-1};
     else
         C = DATA.cluster;
         guess(1:4) = DATA.cluster.xyr;
         guess(5) = DATA.cluster.angle;
%         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
         if isfield(DATA.cluster,'aspectratio')
             state.aspectratio = DATA.cluster.aspectratio;
         else
             state.aspectratio = 1;
         end
     end
     C.y = minmax(DATA.xy{c});
setappdata(DATA.toplevel,'fitparams',[]);
setappdata(DATA.toplevel,'fithists',[]);

maxiter = 100;
state.mintype = 1; %eucliean sum to weight smaller
state.mintype = 2; %use whichever is smallest
state.mintype = 3;
state.initialar = guess(4)./guess(3);
state.rx = guess(3);
state.ry = guess(4);
state.angle = guess(5);
guess = guess([1:3 5]);
%guess = guess([1 2]);
options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
if strcmp(DATA.refinemode,'full')
    dipguess = MinimiseEllipseDip(C, DATA);
    [dpb, bfits] = MinimiseEllipseb(dipguess, DATA,state);
    [dpa, afits] = MinimiseEllipseb(guess, DATA,state);
    if dpb < dpa
        fprintf('Changing starting X,Y\n');
        guess = dipguess;
    end
    nloop = 0;
    while min(afits.r) > 1 && nloop < 10 %no events in cluster
        cc = mean(DATA.xy{state.cluster});
        dvec = guess(1:2)-cc;
        guess(1:2) = guess(1:2)-dvec * min(afits.r)/4;
        [dpa, afits] = MinimiseEllipseb(guess, DATA,state);
        nloop = nloop + 1;
    end
else
    [dp, fits] = MinimiseEllipseb(guess, DATA,state);
end
[fittedparams,fval,exitflag, output] = fminsearch(@MinimiseEllipseb,guess,options,DATA,state);
state.finished = 1;
[dp, fits] = MinimiseEllipseb(fittedparams, DATA,state);
if length(guess) == 2
    C.xyr(1:2) = fittedparams(1:2);
elseif length(guess) == 4
    C.xyr = fittedparams(1:4);
    C.xyr(4) = C.xyr(3).*state.initialar;
    C.angle = fittedparams(4);
else
    C.xyr = fittedparams(1:4);
    C.angle = fittedparams(5);
end

function [SSD,dipr, details ] = EllipseDip(params, DATA, state)
%[SSD, details ] = MinimiseEllipseb(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d
%does not allow change to ellipse aspectratio
cx = params(1);
cy = params(2);
if length(params) == 2
    rx = state.rx;
    ry = state.ry;
    a = state.angle;
elseif length(params) == 4
rx = params(3);
ry = params(3) * state.initialar;
a = params(4);
end
xy = DATA.xy{state.cluster};
xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
a = MyDip(sqrt(r));
SSD = a.dipsize(2);
dipr = a.x(a.dip(2));
details.r = sqrt(r);


function guess = MinimiseEllipseDip(C, DATA)
state.initialar = C.xyr(4)./C.xyr(3);
state.aspectratio = C.aspectratio;
state.rx = C.xyr(3);
state.ry = C.xyr(4);
state.angle = C.angle;
state.cluster = DATA.currentcluster;
guess(1:3) = C.xyr(1:3);
guess(4) = C.angle;
xs = -guess(3)*2:guess(3)/2:guess(3)*2;
ys = -C.xyr(4)*2:C.xyr(4)/2:C.xyr(4)*2;
rs = [C.xyr(3)/2 C.xyr(3) C.xyr(3) .* 2];
for l = 1:length(rs)
    guess(3) = rs(l);
for j = 1:length(xs)
    for k = 1:length(ys)
        guess(1) = C.xyr(1)+xs(j);
        guess(2) = C.xyr(2)+ys(k);
        [dips(j,k,l) dipr(j,k,l)] = EllipseDip(guess, DATA, state);    
    end
end
end
[a,b] = min(dips(:));
[j,k,l] = ind2sub(size(dips),b);
        guess(1) = C.xyr(1)+xs(j);
        guess(2) = C.xyr(2)+ys(k);
        guess(3) = rs(l);
   [a,b,c] =     EllipseDip(guess, DATA, state);
        guess(3) = rs(l).*dipr(j,k);
   [a,b,c] =     EllipseDip(guess, DATA, state);
        
function [SSD, details ] = MinimiseEllipse(params, DATA, state)
%[SSD, details ] = MinimiseEllipse(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d


cx = params(1);
cy = params(2);
rx = params(3);
ry = params(4);
a = params(5);
xy = DATA.xy{state.cluster};
xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
%r = Rprime(r);
details.r = r;

C = DATA.cluster;
C.clst(r < 1) = 2;
C.clst(r>=1) = 1;
[dp, fits] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster

if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (1 - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-1)./fits{2}.sd;
details.fits = fits;
%don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swamp the smaller one
if state.mintype == 1
    dp = -(sqrt(dpa)+sqrt(dpb)).^2;
else
    dp = -(dpa+dpb);
end
if dpa < 0 || dpb < 0
    SSD = 0;
else
    SSD = dp; 
end
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);
rhist = histc(r,linspace(0, max(r)));
fithists = getappdata(DATA.toplevel,'fithists');
fithists = cat(2,fithists,rhist);
setappdata(DATA.toplevel,'fithists',fithists);

function [SSD, details ] = MinimiseEllipseb(params, DATA, state)
%[SSD, details ] = MinimiseEllipseb(params, DATA, state)
%find ellipze that maximizes separation of 2 Gaussian fits in 1d
%does not allow change to ellipse aspectratio
cx = params(1);
cy = params(2);
if length(params) == 2
    rx = state.rx;
    ry = state.ry;
    a = state.angle;
elseif length(params) == 4
rx = params(3);
ry = params(3) * state.initialar;
a = params(4);
end
xy = DATA.xy{state.cluster};
xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);
r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;
%r = Rprime(r);
if isfield(state,'finished')
    details.r = r;
end
details.r = r;

C = DATA.cluster;
C.shape(1) = 0;
[dp, fits, fdetails] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster


if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (1 - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-1)./fits{2}.sd;
details.fits = fits;
details.x = fdetails.x;
if isfield(state,'finished')
    details.xrange = fdetails.xrange;
    details.nbins = fdetails.nbins;
end
%
%if mintype == 1 don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swamp the smaller one
%if mintype == 2, just use whichever dp is the worst. Otherwise can do
%funny things to make dp for inside ellispe very good, and allow dp for
%outsiders to be much owrks
%mintype = 3. Use fitted values at bounary, relative to amplitude. So its
%minimizing the "dip"
if state.mintype == 1
    dp = -(sqrt(dpa)+sqrt(dpb)).^2;
elseif state.mintype == 2
    dp = -min([dpa./sqrt(fits{1}.amp) dpb./sqrt(fits{2}.amp)]);
elseif state.mintype == 3
      if ~isfield(fits{1},'fitted') || ~isfield(fits{2},'fitted')
          dp = 1e6;
      else
          dp = fits{1}.fitted(end)./sqrt(fits{1}.amp)+fits{2}.fitted(1)./sqrt(fits{2}.amp);
          dp = (fits{1}.fitted(end)+fits{2}.fitted(1))./min([fits{1}.amp fits{2}.amp]);
%          dp = -dp;
      end
 else
    dp = -(dpa+dpb);
end
if dpa < 0 || dpb < 0
    if state.mintype == 3
        SSD = 1e6;
    else
        SSD = 0;
    end
else
    SSD = dp; 
end
end
if abs(imag(SSD)) > abs(SSD)/1000 %shouldn't happen
    SSD
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);
rhist = histc(r,linspace(0, max(r)));
fithists = getappdata(DATA.toplevel,'fithists');
fithists = cat(2,fithists,rhist);
setappdata(DATA.toplevel,'fithists',fithists);

function [C, fits] = OptimizeLine(DATA)
    c = DATA.currentcluster;
    state.cluster = c;
    if DATA.currentcluster > 1
        guess(1:4) = DATA.cluster.next{c-1}.xyr;
        guess(5) = DATA.cluster.next{c-1}.angle;
        %         guess(5) = 0;  %DATA.xy{} is rotated by ellipse angle already
        if isfield(DATA.cluster.next{c-1},'aspectratio')
            state.aspectratio = DATA.cluster.next{c-1}.aspectratio;
        else
            state.aspectratio = 1;
        end
        C = DATA.cluster.next{c-1};
    else
        C = DATA.cluster;
        guess(1) = DATA.cluster.crit;
        guess(2) = DATA.cluster.angle;
        state.aspectratio = 1;
    end
    setappdata(DATA.toplevel,'fitparams',[]);
    setappdata(DATA.toplevel,'fithists',[]);
    %First find best starting point - fitting routine likes local minima
    state.mintype = 1; %eucliean sum to weight smaller
    [a,b] = GMDip(DATA.xy{state.cluster},[]);
    crits = [a(1) a(2) DATA.cluster.crit];
    setappdata(DATA.toplevel,'fitparams',[]);
    if strcmp(DATA.refinemode,'test')
        B = MyDip(DATA.xy{state.cluster}(:,1),'idlist',DATA.cluster.clst);
        if B.dip(2)./min(B.dip(3:4)) < 0.97 && B.d > 0.1 & B.d < 0.9
            C.crit =  B.x(B.dip(1));
        end
        fittedparams(1) = C.crit;
        fittedparams(2) = C.angle;
        [dp, fits] = MinimiseLine(fittedparams, DATA,state);
        return;
    end
    if strcmp(DATA.refinemode,'full')
        for j = 1:length(crits)
            guess(1) = crits(j);
            [dps(j), afits] = MinimiseLine(guess, DATA,state);
        end
    else
        j = 0;
    end
    dip = MyDip(DATA.xy{state.cluster}(:,1));
    j = j+1;
    guess(1) = dip.dip(2);
    [dps(j), afits] = MinimiseLine(guess, DATA,state);
    crit(j) = guess(1);
    
    [a,besti] = min(dps);
    guess(1) = crits(besti);
    if min(afits.counts) ==0
        [d, details] = GMDip(afits.r,[]);
        guess(1) = d(1);
    end
    maxiter = 100;
    options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
    nloop = 0;
    %first optimize line position in1D
    %setappdata(DATA.toplevel,'fitparams',[]);
    %[fittedparams,fval,exitflag, output] = fminsearch(@MinimiseLine,guess(1),options,DATA,state);
    %guess(1) = fittedparams(1);
    setappdata(DATA.toplevel,'fitparams',[]);
    [fittedparams,fval,exitflag, output] = fminsearch(@MinimiseLine,guess,options,DATA,state);
    [dp, fits] = MinimiseLine(fittedparams, DATA,state);
    dx = mean(diff(fits.x));
    if fits.counts(2) < C.ncut/5  %Check counts(2) is alwasy the classified group
        PrintMsg(DATA.logfid,sprintf('Optimization Failed nspkds %d -> %d',C.ncut,fits.counts(2)));
        return;
    end

        
        C.crit = fittedparams(1);
        x = fits.fits{1}.mean:dx:fits.fits{2}.mean;
        if length(x) > 2
            a = FitGauss(x,fits.fits{1}.params,fits.fits{1},'eval');
            b = FitGauss(x,fits.fits{2}.params,fits.fits{2},'eval');
            % Best value for fit can be with silly criterion that is far enough
            %into one group that the fit is odd. So if there is a miniumum, put
            %ciriterion there
            [c,d] = min(a+b);
            b = a+b;
            if d > 1 && d < length(x) && c < 0.95 * min(b([1 end]))
                C.crit = x(d);
                fittedparams(1) = x(d);
                [dp, fits] = MinimiseLine(fittedparams, DATA,state);
            end
        end
        C.angle = fittedparams(2);
        C.fitdprime = fits.dpa+fits.dpb;
        C.y = minmax(fits.y);
        

function Plot2GaussFit(params, DATA, state)

    if length(params) <= 2 %line
        C = DATA.cluster;
        C.crit = params(1);
        xy = DATA.xy{state.cluster};
        if length(params) == 2
            a = params(2);
            xys = xyrotate(xy(:,1),xy(:,2),a-C.angle);
        else
            xys = xy;
        end
        
        r = xys(:,1);
    end
    [dp, fits] = Fit2Gauss(C, r, DATA,'plothist')

function [SSD, details ] = MinimiseLine(params, DATA, state)

C = DATA.cluster;
xy = DATA.xy{state.cluster};
if length(params) == 2
a = params(2);
xys = xyrotate(xy(:,1),xy(:,2),a-C.angle);
else
    xys = xy;
end
C.crit = params(1);
crit = C.crit;
r = xys(:,1);
details.r = r;
C.clst(r < C.crit) = 2;
C.clst(r>=C.crit) = 1;
details.counts = [sum(r < C.crit) sum(r >= C.crit)];

[dp, fits, fdetails] = Fit2Gauss(C, r, DATA);
%don't just fit on traditional dprime.  If noise is far away, can
%pay to find samll sd for cluster, even by putting boundary in the
%middle of the cluster
details.x = fdetails.x;
if isempty(fits{1}) || isempty(fits{2})
    SSD = 1e10;
else
dpa = (crit - fits{1}.mean)./fits{1}.sd;
dpb = (fits{2}.mean-crit)./fits{2}.sd;
%dpa and dpb are both positive if fit means are the
%correct side of boundary. If this is not true, the value is
%uninterpretable
details.fits = fits;
details.y = xys(:,2);
details.dpa = dpa;
details.dpb = dpb;
%don't just sum. Favor the smaller number, If one dp is large, don' t let improvements
%in that swampt th esmaller one. Weight accrding to Gaussian, so v large
%sdss become irrelevant
if state.mintype == 1
    if dpa > 0 &&  dpb > 0
        w =  4-[dpa dpb];
        w(w < 0.1) = 0.1;
        dp = -([dpa dpb] * w')./sum(w);
        if dp < fdetails.dprime
           dp = fdetails.dprime;
        end
    elseif dpa > 0
        dp = -dpa/10;
    else
        dp = -dpb/10;
    end
else
    if dpa > 0 &&  dpb > 0
        dp = -(dpa+dpb);
    elseif dpa > 0
        dp = -dpa/10;
    else
        dp = -dpb/10;
    end
end
SSD = dp;
end
fitparams = getappdata(DATA.toplevel,'fitparams');
fitparams = [fitparams; [params dpa dpb dp]];
setappdata(DATA.toplevel,'fitparams',fitparams);

function PlotGauss2Fit(fits);
    
    hist(Rprime(fits.r),fits.nbins);
    hold on;
    n = length(fits.fits{1}.fitted);
    plot(fits.x(fits.xrange(1):fits.xrange(2)),fits.fits{1}.fitted,'r');
    plot(fits.x(fits.xrange(2):fits.xrange(3)),fits.fits{2}.fitted,'r');
        
function [dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)
%[dp, fits, details] = Fit2Gauss(C, r, DATA, varargin)
%negative values of dp are good separation
    plottype = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'plothist',6)
            plottype = varargin{j};
        end
        j = j+1;
    end
    fx = linspace(min(r),max(r),200);
    if isfield(C,'clst')
        id = find(C.clst == DATA.currentcluster+1);
        nid = find(C.clst ~= DATA.currentcluster+1 & C.clst > 0);
    else
    id = find(DATA.clst == DATA.currentcluster+1);
    nid = find(DATA.clst ==1);
    end
    if C.shape(1) == 0 || DATA.currentcluster > 1
        r = Rprime(r);
        fx = linspace(min(r),max(r),200);
        [y,x] = hist(r,fx);
        [a,b] = min(abs(fx-1));
    elseif isfield(C,'crit')
        [a,b] = min(abs(fx-C.crit));
        [y,x] = hist(r,fx);
    else
        [a,b] = min(abs(fx-mean(fx))); %temporary
        [y,x] = hist(r,fx);
    end
    details.dx = mean(diff(fx));
   details.x = x;
    details.y = y;
    details.nbins = 200;
    crit = x(b);
    dp = abs((mean(y(1:b)) - x(b))./std(y(1:b)));
    guess = [mean(y(1:b)) std(y(1:b)) prctile(y(1:b),95)];
    if length(id) <= 1 || length(nid) <= 1
        dp = NaN;
        fits{1}.params = [NaN NaN NaN];
        fits{1}.mean = NaN;
        fits{1}.sd = NaN;
        fits{2}.params = [NaN NaN NaN];
        fits{2}.mean = NaN;
        fits{2}.sd = NaN;
        details.fitpos = [0 0];
        details.dprime = 0;
        details.xrange = [1 1 details.nbins];
        return;
    end
    details.dprime = -abs(mean(r(id))-mean(r(nid)))./sqrt(var(r(id))+var(r(nid)));
% If RHS of distribution is very bimodeal, just fit to the first peak
    lasti = length(y);
    firsti = 1;
    if isfield(C,'sign') && C.sign > 0
        [peak, peaki] = max(y(1:b));
        if peaki > 1
            hid = find(y(1:b) > peak/6);
            lid = find(y(1:hid(end)) < peak/20);
            if ~isempty(lid) && lid(end) > 5
                firsti = lid(end);
            end
        end
    else
        [peak, peaki] = max(y(b:end));
        if b+peaki < lasti
            hid = find(y(b+peaki:end) > peak/6); %don't really need this.
            lid = find(y(b+peaki:end) < peak/20);
            if ~isempty(lid) && lid(1) > 5
                lasti = b+peaki+lid(1);
            end
        end
    end
    details.xrange = [firsti b lasti];
    fits{1} = FitGauss(x(firsti:b),y(firsti:b),'meanlimit',[x(firsti) x(b)]);
    fits{2} = FitGauss(x(b:lasti),y(b:lasti),'meanlimit',[x(b) x(lasti)]);
     if isfield(fits{1},'params') && isfield(fits{2},'params')
         details.fitpos = [fits{1}.params(1) < fx(b) fits{2}.params(1) > fx(b)];
%dont need abs fit 1 ix for x < b, fit 2 is x > b. If mean 1 > mean 2, dprime is bad         
         dp = (fits{1}.params(1)-fits{2}.params(1))./sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
         if sum(details.fitpos) < 2  %one mean is wrong side of criterion. Must be bad
             dp = abs(dp);
         end
         details.diff = abs((fits{1}.params(1)-fits{2}.params(1)));
         details.sd = sqrt((fits{1}.params(2).^2+fits{2}.params(2).^2)/2);
         fya = FitGauss(fx, fits{1}.params, 'eval');
         fyb = FitGauss(fx, fits{2}.params, 'eval');
         fy = fya+fyb;
         id = find(fx > fits{1}.mean & fx < fits{2}.mean);
         if ~isempty(id)
             [a,b] = min(fy(id));
             details.minxpt = fx(id(b));
         else
             details.minxpt = NaN;
         end
         if strcmp(plottype,'plothist')
             GetFigure('FitGauss');
             hold off;
             bar(x,y);
             hold on;
             bar(x(1:b),y(1:b),'r')
             plot(fx,fya,'r-','linewidth',2);
             plot(fx,fyb,'b-','linewidth',2);
             set(gca,'ylim',[0 max(y)]);
             title(sprintf('Dprime from fits %.2f',dp));
         end
     else
         if isempty(fits{1})
             fits{1}.params = [NaN NaN NaN];
             fits{1}.mean = NaN;
             fits{1}.sd = NaN;
         end
         if isempty(fits{2})
             fits{2}.params = [NaN NaN NaN];
             fits{2}.mean = NaN;
             fits{2}.sd = NaN;
         end
         dp = NaN;
         details.fitpos = [0 0];
         details.diff = NaN;
         details.sd = NaN;
     end
     details.fits = fits;

     
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
    
function    AddComment(a,b,str)
    
    DATA = GetDataFromFig(a);
    n = length(DATA.Comments) + 1;
    DATA.Comments(n).ex = DATA.exptno;
    DATA.Comments(n).p = ProbeNumber(DATA);
    DATA.Comments(n).exptno = DATA.exptno;
    
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

    
    
function C= CondenseClusters(C, go, varargin)

allfields = {};
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
    allfields = unique({allfields{:} f{:}});
    end
end
for j = 1:length(C)
    if ~isempty(C{j})
    f = fields(C{j});
        for k = 1:length(f)
            CC(j).(f{k}) = C{j}.(f{k});
        end
    end
end
if go == 0
    return;
end
C = CC;
    

function p = GetProbeFromName(name)
p = 0;
id = regexp(name,'\.p[0-9]*');
if ~isempty(id)
    p = sscanf(name(id(1)+2:end),'%d');
end


    
function ShowTaggedProbes(DATA)

    if ~isfield(DATA,'TaggedProbes') || sum(DATA.TaggedProbes) == 0 || DATA.interactive <1
        return;
    end
    tmenu = findobj(DATA.toplevel,'Tag','ProbeSwitchMenu');
    c = get(tmenu,'Children');
    id = find(DATA.TaggedProbes > 0);
    for j = 1:length(id)
        set(c(id(j)),'foregroundcolor',DATA.colors{DATA.TaggedProbes(id(j))});
    end
    
function DATA = QuickAutoCut(a,b)
    DATA = GetDataFromFig(a);
    AllV = GetAllV(DATA);
    tpt(1) = find(DATA.spts == 0);
    p = DATA.probe(1);
    tpt(2) = tpt(1)+5;
    sz = abs(AllV(p,tpt(1),:) - AllV(p,tpt(2),:));
    nspk = size(AllV,3);
    if nspk > 10000
        prc =99.5;
    else
        prc = 98;
    end
    id = find(sz > prctile(sz,prc));
    DATA.cluster.clst = ones(1,size(AllV,3));
    DATA.cluster.clst(id) = 2;
    DATA.clst = DATA.cluster.clst;
    DATA.cluster.neednewtemplate = 1;
    TemplatePlot(DATA);
    DATA = get(DATA.toplevel,'UserData');
    %1 is r, 8 is dt
    x = (DATA.TemplateScores(:,1) + DATA.TemplateScores(:,8))./sqrt(2);
    a = FindDip(x);
    DATA.cluster.crit = a(4);
    DATA.cluster.angle = -pi/4;
    DATA.cluster.shape = 1;
    DATA.cluster.space = [3 1 8];
    [cl, DATA.cluster, DATA.xy{1}] = ClassifySpikes(DATA, DATA.cluster);
    SetFigure(DATA.tag.hist,DATA);
    hold off;
    hist(x,500);
    hold on;
    plot([a(4) a(4)],get(gca,'ylim'),'r');
    set(DATA.toplevel,'UserData',DATA);
    
    