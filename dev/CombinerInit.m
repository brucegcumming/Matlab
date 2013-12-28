function DATA = CombinerInit(DATA, name, TOPTAG, layout)


DATA.layout.top = [];
if ~isfield(DATA,'nevdir')
    DATA.nevdir = [];
end
    DATA.layout.spkv = [];
    DATA.layout.spkxy = [];
    DATA.layout.lfp = [];
    DATA.layout.combineplot = [];
    DATA.layout.options = [];
    DATA.gui.prefsdir = [GetFilePath('preferences') '/Combine/'];
    DATA.defaultconfig = [DATA.gui.prefsdir '/' gethostname '.' GetUserName '.config'];
    DATA.defaultlayout = [DATA.gui.prefsdir '/' gethostname '.' GetUserName '.layout.mat'];
    if ~isfield(DATA,'configfile') || isempty(DATA.configfile)
        DATA.configfile = DATA.defaultconfig;
    end
    if ~isfield(DATA,'layoutfile') || isempty(DATA.layoutfile)
        DATA.layoutfile = DATA.defaultlayout;
    end
    DATA.bysuffix = 0;
    if ~isempty(layout)
        f = fields(layout);
        for j = 1:length(fields)
            DATA.layout.(f{j}) = layout.(f{j});
        end
    end
DATA.expstrs = {'dx' 'xo' 'dp' 'dO'  'dxP' 'dxXceC' 'dxXce' 'sf' 'dxXId' 'dxXIdP' 'or' 'Op' 'Pp' 'sz' 'jv'...
    'orXob' 'orXobP' 'orXme' 'szXob' 'DcXorP' 'DcXor' 'DcXorXUsPRC' 'DcXorXUsP' 'OpXme' 'PpXme' 'co' 'CtfXip' 'tfXip' 'IB' 'backMov' 'backMovXannTyp' 'IBXannTyp' 'dirXme' ...
    'OpXip' 'PpXip' 'me' 'OlXOr' 'ce' 'dxXUs' 'dxXUsP' 'dxXpi' 'dxXor' 'orXdx' 'TwoCylDispXhxPRB' 'TwoCylDispXhxP' 'TwoCylDispXhx' 'dxXIdPD' 'dxXIdD'  ...
    'orXobXUs' 'orXobXUsP' 'dxXceXFr' 'yo' 'orXobXfyP' 'orXobXfy'...
    'orXobXfxP' 'orXobXfx' 'orXId' 'dirXId' 'orXIdD' 'dirXIdD' ...
    'szXobXme' 'orXip' 'dxXceXrC' 'meXFr' 'meXFrRC' 'dxXces' 'dxXcesRC' 'dxXmixac' 'dxXcesXFr' ...
    'dxPRB' 'dxRB' 'bdXIdP' 'dxXIdPDB' 'sMXId' 'sMXdx' ...
    'sMXPd' 'TwoCylDispXhxXbhPTWO' 'TwoCylDispXhxXbhPDID'  'jvXsz' 'sOXor' 'sOXorXar' 'fxXfy' ...
    'jvXpi' 'dxXdy' 'dxXst' 'OpXdx' 'PpXdx' 'szXdx' 'wiXdx' 'OpXdxnoback' 'PpXdxnoback'};
DATA.expnames = {'DT', 'XO', 'DP' 'ODX' 'DTP' 'AC' 'AC' 'SF' 'DTID' 'ABD' 'OT' 'OP' 'PP' 'SZ' 'VE' ...
    'OTOB' 'ORBW' 'OXM' 'SZOB' 'DCOR' 'DCOR' 'UDCOR' 'UDCOR' 'OPM' 'PPM' 'CO' 'CTFIP' 'TFIP' 'IB' 'Movie' 'MovieANN' 'IB' ....
    'OD' 'OPIP' 'PPIP' 'OD' 'BOP' 'CE' 'UADT' 'UADT' 'DPI' 'DXO' 'DXO' 'TWO' 'TWO' 'TWO' 'DID' 'DID' 'UORBW' 'UORBW' 'FAC' 'YO' ...
    'ORBW' 'ORBW' 'ORBW' 'ORBW' 'SRID' 'SRID' 'DRID' 'DRID' ...
    'SFMOB' 'ORP' 'AC' 'FOD' 'FODRC' 'XAC' 'XACRC' 'MIXAC' 'FXAC' 'DXRB' 'DXRB' 'BDID' 'DIDB' 'FLID'...
    'DTPL' 'PDPL' 'TWO' 'DID' 'VESZ' 'OPPP' 'OPPP' 'EyeCal' 'JPI' 'DXDY' 'DXST' 'OPDXB' 'PPDXB' 'SZDX' 'WIDX' 'OPDX' 'PPDX'};
DATA.expstrs = {DATA.expstrs{:} 'sfXob' 'sfXme' 'szXme' 'dxXdx' 'dOXdO' 'dOXceXFr' 'dOXce' 'orRC','pRP' 'szRC'};
DATA.expnames = {DATA.expnames{:} 'SFOB' 'SFM' 'SZM' 'DT' 'ODX' 'FOXAC' 'OXAC' 'OTRC' 'CR' 'SZRC'};
DATA.expstrs = {DATA.expstrs{:} 'dxXrb' 'dxXrbP'};
DATA.expnames = {DATA.expnames{:} 'DTRW' 'DTRW'};
DATA.monkeys = {'lem' 'ica' 'jbe' 'dae' 'jbe' 'bgy'};
DATA.monkey = 'none';
DATA.stimnames = {'none' 'grating' 'bar' 'cylinder' 'rds' 'rls' 'image'};
DATA.stimnamecodes = [0 3 4 11 2 15 21];
    DATA.subprobe = 0;
    DATA.subprobes = 0;
    DATA.xprobes = [];
DATA.listbycell = 0;
DATA.usenev = 0;


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
    DATA.tag.rcfigc= [TOPTAG 'RCc'];
    DATA.tag.showvals= [TOPTAG 'ShowVal'];
    DATA.tag.celllist= [TOPTAG 'CellList'];
    DATA.tag.psych = [TOPTAG 'Psych'];
    DATA.tag.cptag = [TOPTAG 'CombineCP'];
    DATA.tag.allexpts = [TOPTAG 'AllExpts'];
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.tag.dataplot = [TOPTAG 'Combiner'];
    DATA.appending = 0;
    DATA.appendexptid = NaN;
    
    DATA.wsc= [1 1];
    DATA.defaults.fz = 96;
    DATA.state.fixrange = 0;
    DATA.spooling = 0;
    DATA.densityplot = 0;
    DATA.alldensityplot = 0;
    DATA.plot.showspikeshapes = 0;
    DATA.plot.showallxy = 0;
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
    DATA.show.ntrials = 0;
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
    DATA.plot.nminrc = 10; %= include all by default
    DATA.plot.sdfw = 166;
    DATA.plot.showISI = 0;
    DATA.plot.quickspks = 0;
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
    DATA.probelist = [1];
    DATA.probe = 1;
    DATA.oldprobe = 1;
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
    DATA.profiling = 0;
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
    DATA.state.usensx = 0;
    DATA.state.somespikes = 0; %set when loading only spikes for some expts/probes - eg griddata
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
    DATA.progname = 'Combine';