function DATA = ReadFile(DATA, name, varargin)

args = {};
j = 1;
while j <= length(varargin)
args = {args{:} varargin{j}};
j =j+1;
end
DATA.datafilename = name;
DATA.state.online = 0;
[a,b] = fileparts(fileparts(name));
DATA.prefix = b;
DATA.monkey = GetMonkeyName(name);

DATA.logfid = fopen(strrep(name,'.mat', '.log'),'a');
ts = now;
if isfield(DATA,'toplevel')
set(findobj(DATA.toplevel, 'Tag','FileName'),'string',name);
end
if strfind(name,'.txt') %% an online text file
Expts = ReadOnlineTxt(name, args{:});
Expts = cmb.CountTxtSpikes(Expts,DATA.probe,DATA.spikelist);
DATA.state.online = 2;
DATA.lastread = 0;
DATA.lastsize = 0;
if ~isfield(DATA,'timerobj')
DATA.timerobj = timer('TimerFcn',@cmb.timerfn, 'Period', 1.0, ...
'Tag',DATA.tag.top, 'ExecutionMode','FixedSpacing');
DATA.lastread = 0;
end
else
ts = now;
[Trials, Expts, All] = APlaySpkFile(name, args{:});
fprintf('APlaySpkFile Took %.2f\n',mytoc(ts));
if isfield(Trials,'DataType')
DATA.filetype = Trials.DataType;
end
DATA.state.nospikes = Trials.state.nospikes;
if isfield(All,'Events')
sonid = find(All.Events.codes(:,1) ==48); %'0' = storage on
soffid = find(All.Events.codes(:,1) ==49); %'1' = storage off
else
sonid = [];
soffid = [];
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
if DATA.state.usensx
DATA.state.somespikes = 2;
else
DATA.state.somespikes = 1;
end
DATA.state.nospikes = 2;
DATA.probe = DATA.probelist(1);
DATA.AllClusters = {};
for j = 1:np
DATA.probes(j).probe = j;
end
if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
DATA.probe = DATA.probelist(1);
end
DATA.ArrayConfig = GetArrayConfig(name);

if np > 1
DATA = cmb.AddMultiProbeGUI(DATA);
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
DATA = cmb.AddMultiProbeGUI(DATA);
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
DATA.probesource = cmb.FindProbeSources(DATA);
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
elseif DATA.bysuffix
[a,b] = fileparts(name);
DATA.cellfile = [a '/CellList.mat']
elseif strncmp(DATA.filetype,'GridData',8)
[a,b] = fileparts(name);
DATA.cellfile = [a '/CellList.mat'];
else
DATA.cellfile = strrep(name,'.mat','.cells.mat');
end

end
DATA = cmb.cListExpts(DATA,Expts);

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
% when it all comes from one file, should be no need for trialCaoffset or suffix.        
Expt{nexp}.Header.trialoffset = 0;
Expt{nexp}.Header.suffix= 0;

if length(sonid)
id = find(All.Events.times(sonid) > Expts{nexp}.Header.trange(1) & ...
All.Events.times(sonid) < Expts{nexp}.Header.trange(2));
Expts{nexp}.gridstoreon = All.Events.times(sonid(id));
id = find(All.Events.times(soffid) > Expts{nexp}.Header.trange(1) & ...
All.Events.times(soffid) < Expts{nexp}.Header.trange(2));
Expts{nexp}.gridstoreoff = All.Events.times(soffid(id));
end
if isfield(All,'DigMark')
id = find(All.DigMark.times > Expts{nexp}.Header.trange(1) & ...
All.DigMark.times < Expts{nexp}.Header.trange(2));
Expts{nexp}.DigMark.times = All.DigMark.times(id);
Expts{nexp}.DigMark.codes = All.DigMark.codes(id,:);
end
end
if strncmp(DATA.filetype,'Grid',4)
DATA.grididx = BuildGridIndex(name, Expts);
DATA.ArrayConfig = GetArrayConfig(name);
for j = 1:length(DATA.ArrayConfig.id)
DATA.probenames{j} = sprintf('%d(E%d %d,%d)',j,DATA.ArrayConfig.id(j),DATA.ArrayConfig.X(j),DATA.ArrayConfig.Y(j));
end
DATA.state.includeprobename = 1;
end
if isempty(Expts)
cmb.NotBusy(DATA);
return;
end
if DATA.appending 
DATA.Expts = [DATA.Expts Expts];
else
DATA.Expts = Expts;
end
if isempty(DATA.exabsid)
DATA.exabsid = 1:length(DATA.Expts);
end
%need to load cell file AFTER listing expts so that can check exptids
%match
DATA = cmb.LoadCellFile(DATA);
if DATA.state.online
DATA.AllData = [];
else
if isfield(All,'AllClusters')
DATA.AllClusters = All.AllClusters;
All = rmfield(All,'AllClusters');
DATA.state.nospikes = 2;
end
if isfield(All,'AllSpikes')
DATA.AllSpikes = All.AllSpikes;
All.Spikes = [];
All = rmfield(All,'AllSpikes');
end
DATA.AllData = All;
if isfield(Trials,'Spkid')
DATA.AllData.SpikeIdx = Trials.Spkid;
end
DATA.AllData.Trialids = Trialids;
DATA.AllData.TrialTimes = TrialTimes;
DATA.AllData.pcs = [];
if isfield(Trials,'Comments')
DATA.Comments = Trials.Comments;
end

%if AllClusteres, means loding xy, times for many probes. Don't load
%the spikes files at startup for default probe

if ~isfield(DATA,'AllClusters')
DATA = cmb.SetSpkLists(DATA);
end
% load saved cluster params

if DATA.state.useonlineclusters
DATA = cmb.LoadOnlineClusters(DATA, cmb.ClusterFile(DATA,'getonline'));
end
DATA = cmb.LoadClusters(DATA, cmb.ClusterFile(DATA));
DATA.name = name;
if isdir(DATA.name)
DATA.datadir = DATA.name;
else
DATA.datadir = fileparts(DATA.name);
end

if isfield(DATA,'toplevel')
DATA = cmb.ListSubExpts(DATA,1,'relist');
end
if length(DATA.probelist) > 1 && DATA.state.NeedAllProbeFile
DATA = cmb.LoadClusters(DATA, cmb.ClusterFile(DATA,'allprobes'),'allprobes');
end
cfile = cmb.CombinerLst(DATA);
if exist(cfile,'file')
load(cfile);
DATA.combines = combines;
end
if isfield(Trials,'setprobe') & length(Trials.setprobe) > 1
DATA.plot.useprobe = zeros(size(DATA.probelist));
DATA.plot.useprobe(Trials.setprobe) = 1;
DATA.currentexpt = 1;
DATA = cmb.LoadAllProbeSpikes(DATA, 0, 'select');
else
%        DATA.plot.useprobe = zeros(size(DATA.probelist));
%        DATA.plot.useprobe(DATA.probe) = 1;
end


end
ShowExptErrs(DATA.Expts)


