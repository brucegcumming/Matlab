function varargout = PlotClusters(name, varargin)
%PlotClusters(dir, ......)
%reads in Cluster files made by AllVPcs and plots up clusters for each
%expt/cell
if isunix
 prefsdir = '/b/group/matlab/preferences/PlotClusters';
else
 if isdir('Z:/group')
  prefsdir = 'Z:/group/matlab/preferences/PlotClusters';
 else
  prefsdir = '/b/group/matlab/preferences/PlotClusters';
 end
end
defaultconfig = [prefsdir '/' gethostname '.' GetUserName '.config'];
defaultlayout = [prefsdir '/' gethostname '.' GetUserName '.layout.mat'];
TAGTOP = 'PlotClusters';
loadargs = {};
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
warnmode = 'nowarn';
expts = [];

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
 it = FindFig(TAGTOP);
 if isempty(it)
  initcall = 1;
  DATA.tag.top = TAGTOP;
  DATA = PC.SetDefaults(DATA);
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
 DATA = PC.InitInterface(DATA); 
 DATA.defaultlayoutfile = defaultlayout;
 DATA.defaultconfigfile = defaultconfig;
 PC.ApplyLayout(defaultlayout);
 PC.ApplyConfig(defaultconfig);
 [DATA, Clusters] = PC.ReadClusterResults(DATA, Clusters);
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
 DATA = PC.SetDefaults(DATA);
 DATA = PC.InitInterface(DATA);
elseif strncmpi(name,'close',4) %must be in front of exist test because of close.m
 f = fields(DATA.tag);
 for j = 1:length(f)
  PC.CloseTag(DATA.tag.(f{j}));
 end
 return;
elseif strncmpi(name,'profile',4) %must be in front of exist test because of close.m
 DATA.profiling = ~DATA.profiling;
 fprintf('Profiling %d\n',DATA.profiling);
 
 set(DATA.toplevel,'UserData',DATA);
 return;
elseif exist(name,'file')
 DATA.tag.top = TAGTOP;
 DATA = PC.SetDefaults(DATA);
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
 DATA.strings = strings;
 DATA.dates = zeros(size(DATA.strings));

 DATA = PC.InitInterface(DATA);
end

if strcmp(name,'getstate')
 varargout{1} = DATA;
 return;
end
j = 1;
while j <= length(varargin)
 if strncmpi(varargin{j},'LoadAuto',8)
  DATA = PC.LoadAll(DATA,'auto','force');
  set(DATA.lstui,'string',DATA.strings);
 elseif strncmpi(varargin{j},'config',4)
  j = j+1;
  if regexp(varargin{j},'^[A-Z]:') | strfind(varargin{j},'/') %real path
DATA.configfile = varargin{j};
  else
DATA.configfile = [prefsdir '/' varargin{j} '.config'];
  end
  if exist(DATA.configfile)
DATA = PC.ApplyConfig(DATA);
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
  PC.ApplyLayout(DATA);
  end
 elseif strncmpi(varargin{j},'expts',5)
  j = j+1;
  expts = varargin{j};
 elseif strncmpi(varargin{j},'exptload',7)
  DATA = LoadExpts(DATA);
 elseif strncmpi(varargin{j},'Load',4)
  
  if ~isfield(DATA,'name') && ischar(name)
fprintf('%s is not a directory\n',name);
return;
  elseif isdir(DATA.name)
 if strncmpi(varargin{j},'LoadSpikes',8)
DATA = PC.LoadAll(DATA,[],'loadspikes',warnmode);
 else
DATA = PC.LoadAll(DATA,[],'force', warnmode,'expts',expts,loadargs{:});
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
%  DATA.templatesrc = 1; 
 elseif strncmpi(varargin{j},'templateid',10)
  j = j+1;
  DATA.templateid = varargin{j};
 elseif strncmpi(varargin{j},'useauto',7)
  DATA.useautoclusters = 1;
 elseif strncmpi(varargin{j},'usealltrials',7)
  loadargs = {loadargs{:} varargin{j}};
 end
 j = j+1;
end

if strncmpi(name,'close',4)
 f = fields(DATA.tag);
 for j = 1:length(f)
  PC.CloseTag(DATA.tag.(f{j}));
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
 PC.CheckExclusion(DATA);
elseif strncmpi(name,'convertlist',7)
 DATA = PC.ConvertExclusion(DATA);
 if strcmp(name,'convertlistsave')
  SaveCellList(DATA);
 end
elseif strncmpi(name,'correlogram',7)
 PlotCorrelogram(Clusters{DATA.currentpoint(1)}{DATA.currentpoint(2)});
elseif strncmpi(name,'findcells',7)
 CellFinderTest(DATA,1);
elseif strncmpi(name,'quickspikes',7)
% QuickSpikes(DATA, DATA.currentpoint);
 PC.SpoolSpikes(DATA, DATA.currentpoint);
elseif strncmpi(name,'reloadtrials',7)
 DATA = LoadTrialLists(DATA);
elseif strncmpi(name,'loadextra',7)
 DATA = PC.LoadExtra(DATA,0);
elseif strncmpi(name,'reloadcellist',9)
 DATA = PC.LoadCellFile(DATA);
elseif strncmpi(name,'readres',7)
 [DATA, Clusters] = PC.ReadClusterResults(DATA, Clusters);
 setappdata(DATA.toplevel,'Clusters',Clusters);
elseif strncmpi(name,'PC.CallAllVPcs',7)


elseif strncmpi(name,'loadspikes',7)
 PC.PC.LoadAllSpikes(DATA);

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
 PC.OptionMenu(DATA, [], name);
 return; %Don't set UserData below
elseif strncmpi(name,'setentry',7)
 id = get(gui,'value');
 DATA.id = id;
 Clusters = getappdata(DATA.toplevel,'Clusters');
 if id > length(Clusters) || isempty(Clusters{id})
 Clusters{id} = LoadCluster([DATA.name '/' DATA.strings{id}],'getxy');
 end
 figure(DATA.fig.clusters);
 DATA = doentry(DATA, Clusters, id);
elseif strncmpi(name,'followcorr',7)
 DATA.plot.alltype = 'followcorr';
 DATA.templateid = fliplr(DATA.currentpoint);
 DATA = PC.PlotAllClusters(DATA,[]);
elseif strncmpi(name,'MakeCells',6)
  cellps = getappdata(DATA.toplevel,'CellPs');
  MakeCellId(cellps);
elseif strncmpi(name,'plotexpts',6)
  PlotExpts(DATA);
elseif strncmpi(name,'plottest',6) %for testing new plots
 PC.OptionMenu(DATA,[],'xcorrallprobes');
 return;
 
 PC.OptionMenu(DATA,[],'findconnect');
 eid = DATA.currentpoint(1);
 probes = find(DATA.selectprobe(eid,:));
 PC.SetFigure(DATA,DATA.tag.spikes,'front');
 PC.PlotSyncSpikes(DATA, eid, fliplr(probes),[1 1], varargin{:});
 return;
 DATA.plotcells.showmahal = 2;
 DATA.plotcells.showfitdp = 2;
 PC.PlotCellList(DATA, 'showfig');
elseif strncmpi(name,'plotxcorrprobes',12)
 PC.SetFigure(DATA, DATA.tag.xcorr);
 PlotXcorrs(DATA, DATA.xcorrs, 1:length(DATA.exptid), 0);
elseif strncmpi(name,'plotxcorr',6)
 PC.SetFigure(DATA, DATA.tag.xcorr);
 PlotXcorrs(DATA, DATA.xcorrs, 1:length(DATA.exptid),1);
elseif strncmpi(name,'popplot',7)
 if length(varargin) && ischar(varargin{1})
  DATA.plot.alltype = varargin{1};
 end
 DATA = PC.PlotAllClusters(DATA,[]);
elseif strncmpi(name,'plotmenu',7)
 PlotMenu(DATA,[], varargin{1}, varargin{2});
elseif strncmpi(name,'setplot',7)
 if gui
  strs = get(gui,'string');
 DATA.plot.onetype = deblank(strs(get(gui,'value'),:));
 end
 DATA = doentry(DATA, Clusters, DATA.id);
elseif strncmpi(name,'showdata',7)
 PC.ShowData(DATA, DATA.currentpoint(1),DATA.currentpoint(2));
end
if initcall 
 if isfield(DATA,'exptid') %Do pop plot if data loaded
 DATA = PC.PlotAllClusters(DATA,[]);
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
 DATA = PC.ApplyLayout(DATA, defaultlayout);
 DATA= PC.ApplyConfig(DATA, defaultconfig);
 DATA = PC.LoadCellFile(DATA);
 DATA = PC.LoadExtra(DATA,0);
 DATA.ArrayConfig = GetArrayConfig(DATA.name);
 if isfield(DATA.ArrayConfig,'X')
  DATA.nprobes = length(DATA.ArrayConfig.X);
 elseif size(DATA.CellList,2) > 0
  DATA.nprobes = size(DATA.CellList,2);
 end
 DATA = PC.PlotCellList(DATA,'showfig');
 if ~isfield(DATA,'exptid')
  DATA.exptid = DATA.CellDetails.exptids;
%  DATA = LoadExpts(DATA);
  for j = 1:length(DATA.exptid)
Clusters{j} = {};
  end
  setappdata(DATA.toplevel,'Clusters',Clusters);
 end
 DATA.nexpts = length(DATA.exptid);
end

PC.SetGUI(DATA);
set(DATA.toplevel,'UserData',DATA);


