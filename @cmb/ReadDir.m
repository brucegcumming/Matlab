function DATA = ReadDir(DATA, name, varargin)  
%ReadDir(DATA, name, varargin) read expts from dir (online, or split files)

d = dir(name);
[a,b] =sort([d.datenum]);
nx = 1;
nf = 1;
for j = 1:length(d)
    if strfind(d(j).name,'ClusterTimes')
        nx = nx+1;
    elseif regexp(d(j).name,'A.[0-9]*.mat');
        nf = nf+1;
    end
end

if nx > 5 %temporary - need to cherk expected number
    DATA.bysuffix = 1;
    DATA.state.online = 0;
end

functs = now;
AllExpts = {};
d = d(b);
reindex = 0;
args = {};
j = 1;
nskip = 0;
if name(end) == '/' || name(end) == '\'
    name = name(1:end-1);
end
[dirname, nameroot] = fileparts(name);
if DATA.bysuffix
    DATA.cellfile = [name, '/CellList.mat'];
    DATA.monkey = GetMonkeyName(name);
    DATA.prefix = [DATA.monkey nameroot];
else
    DATA.cellfile = [name, '/Cells.mat'];
end
DATA = cmb.LoadCellFile(DATA);
DATA.datadir = name;
if isempty(DATA.nevdir)
    DATA.nevdir = name;
end

quicksuffix = 0;
dirtype = expt.Type(name);
if strcmp(dirtype,'online')
    DATA.state.online = 1;
else
    DATA.state.online = 0; %will be overridden by commandline below
end

j = 1;
while j <= nargin-2
    if strncmpi(varargin{j},'quicksuffixtest',13)
        quicksuffix = 1;
    elseif strncmpi(varargin{j},'quick',5)
        quicksuffix = 1;
        DATA.bysuffix = 1;
        DATA.AllClusters = {};
        DATA.currentexpt = 1;
    elseif strncmpi(varargin{j},'relist',3)
        reindex =1;
        args = {args{:} 'relist'};
    elseif strncmpi(varargin{j},'online',6)
        DATA.bysuffix = 0;
        DATA.state.nospikes = 0;
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
    DATA.logfid = fopen([name,'/Combinelog.txt'],'a');
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
suffixlist = [];
suffixid = [];
DATA.ArrayConfig = GetArrayConfig(DATA.nevdir);

if DATA.state.online
    method = 1;
else
    method = 2;
end

if quicksuffix || method == 2
    DATA = cmb.ReadFinishedDir(DATA, name, varargin{:});
    fprintf('Quick Read took %.2f\n',mytoc(functs));
    if DATA.state.online
        DATA.bysuffix = 0;
    end
    if quicksuffix
    return;
    end
end
ts = now;
if method == 2
    oldstate = DATA.state;
    DATA.state.showspikes = 0;
    DATA.AllClusters = {};
    sumload = 0;
    for e = 1:length(DATA.Expts)
       DATA.currentexpt = e;
       fprintf('Loading Cluster Times for Expt %d (#%d)\n',GetExptNumber(DATA.Expts{e}),e);
       [DATA, details] =  cmb.CheckClusterLoaded(DATA,e);
       sumload = sumload + details.loaddur(1);
    end
    nexp = length(DATA.AllClusters);
    DATA.state.showspikes = oldstate.showspikes;
    Expts = DATA.Expts;
else
    for j = 1:length(d)
        if (length(regexp(d(j).name,'Expt[0-9]*.mat')) | length(regexp(d(j).name,'\.[0-9]*\.mat'))) & ...
                d(j).bytes > 128 & .....
                isempty(strfind(d(j).name,'idx.mat')) & ...    %exclude the .idx file
                d(j).datenum < now-(1/(24 * 60 * 60)) & ... %at least  1 sec old
                isempty(strmatch(d(j).name,{expnames{:}}))%dont read if we already have
            nexpts = nexpts+1;
            id = regexp(d(j).name,'\.[0-9]*\.mat');
            if length(id)
                suffix = sscanf(d(j).name(id(1)+1:end),'%d');
                suffixlist = [suffixlist suffix];
                suffixid = [suffixid j];
            end
        elseif regexp(d(j).name,[nameroot '.*\.[0-9]*.mat'])
            nexpts = nexpts+1;
            %            fprintf('%s\n',d(j).name);
        end
    end
    if DATA.bysuffix
        [suffixlist,sid] = sort(suffixlist);
        sid = suffixid(sid);
        suffixlist = unique(suffixlist);
        d = d(sid);
        DATA.suffixlist  = suffixlist;
    end
fprintf('%d Expts in %s (check took %.1f)\n',nexpts,name,mytoc(ts));
sumload = 0;
ts = now;
%Expts = {}; bad for online files
for j = 1:length(d)
    if (length(regexp(d(j).name,'Expt[0-9]*.mat')) | length(regexp(d(j).name,[nameroot '*\.[0-9]*.mat']))) & ...
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
            fprintf('Reading %s\n',d(j).name);
            if DATA.bysuffix
                args = {args{:}, 'bysuffix'};
            else
                args = {'online' args{:}};
            end
            ts = now;
            if isfield(DATA,'toplevel') && isfigure(DATA.toplevel)
                DATA.defaults.toplevel = DATA.toplevel;
            end
            [trls, exps, All] = APlaySpkFile([name '/' d(j).name],'Defaults',DATA.defaults,args{:});
            %            fprintf('APlaySpkFile took %.2f\n',mytoc(ts));
            if DATA.state.autocutatstart
                if isfield(All,'Spike2Ch')
                    setappdata(DATA.toplevel,'SpikeCh',All.Spike2Ch);
                end
            end
            if isfield(All,'Spikes') && ~isfield(All.Spikes,'exptid')
                All.Spikes.exptid = GetExptNumber(d(j).name);
            end
            if isfield(trls,'loaddur')
                sumload = sumload + trls.loaddur(1);
            end
            soffid = [];
            sonid = [];
            if isfield(trls,'DataType')
                DATA.filetype = trls.DataType;
                if strmatch('GridData',trls.DataType);
                    sonid = find(All.Events.codes(:,1) ==48); %'0' = storage on
                    soffid = find(All.Events.codes(:,1) ==49); %'1' = storage off
                    
                end
            end
            
            if isempty(exps)
                fprintf('%s No expts\n',d(j).name);
                DATA.badnames{nbad} = d(j).name;
                DATA.badtimes(nbad) = d(j).datenum;
                nbad = nbad+1;
                badidx = 1:nbad-1;
                suffix = GetExptNumber(d(j).name);
                if isfield(DATA,'suffixlist')
                    DATA.suffixlist(DATA.suffixlist == suffix) = []
                end
            elseif iscell(exps)
                if length(exps) > 1
                    nt = [];
                    for k = 1:length(exps)
                        nt(k) = length(exps{k}.Trials);
                    end
                    [a,b] = max(nt);
                    if b > 1
                        fprintf('Using Expt %d\n',b);
                    end
                    Expts{nexp} = exps{b};
                else
                    Expts{nexp} = exps{1};
                end
                Expts{nexp}.gui.classified = 0;
                Expts{nexp}.gui.counted = 0;
                Expts{nexp}.gui.clustertype = 0;
                Expts{nexp}.gui.firsttrial = 1+length(Trialids);
                Expts{nexp}.errs = trls.errs;
                if length(sonid)
                    Expts{nexp}.gridstoreon = All.Events.times(sonid);
                    if length(soffid)
                        Expts{nexp}.gridstoreoff = All.Events.times(soffid);
                    else
                        Expts{nexp}.gridstoreoff = 0;
                    end
                end
                if isfield(All,'DigMark')
                    Expts{nexp}.DigMark = All.DigMark;
                end
                if isfield(trls,'bstimes')
                    id = find(trls.bstimes > Expts{nexp}.Header.Start & trls.bstimes < Expts{nexp}.Header.End);
                    Expts{nexp}.bstimes = trls.bstimes(id);
                end
                if isfield(trls,'ClusterLoadTimes')
                    Expts{nexp}.ClusterLoadTimes = trls.ClusterLoadTimes;
                    sumload = sumload + sum(trls.ClusterLoadTimes);
                    sumloads(nexp,1:length(trls.ClusterLoadTimes)) = trls.ClusterLoadTimes;
                end
                SpkId = [SpkId; trls.Spkid];
                for k = 1:length(Expts{nexp}.Trials)
                    news(k) = Expts{nexp}.Trials(k).Start(1);
                end
                
                TrialTimes = [TrialTimes news];
                Expts{nexp}.gui.ntrials = length(news);
                if DATA.state.online == 1 && isfield(All,'Spikes');
                    Spikes = All.Spikes;
                elseif nexp > 1 && isfield(Spikes,'values') && isfield(All.Spikes,'values') && size(Spikes.values,2) == size(All.Spikes.values,2)
                    Spikes.values = cat(1,Spikes.values, All.Spikes.values);
                    if isfield(All.Spikes,'dVdt')
                        Spikes.dVdt = cat(1,Spikes.dVdt, All.Spikes.dVdt);
                    end
                    Spikes.codes = [Spikes.codes; All.Spikes.codes];
                    Spikes.times = [Spikes.times; All.Spikes.times];
                elseif isfield(All,'AllClusters')
                    DATA.AllClusters{nexp} = All.AllClusters;
                    DATA.Clusterfile{nexp}.quick = All.quickload;
                    DATA.Clusterfile{nexp}.loadtime = now;
                    DATA.Clusterfile{nexp} = CopyFields(DATA.Clusterfile{nexp},All, {'datenum' 'loadname'});
                    
                    if isfield(All,'FullVData')
                        DATA.Clusterfile{nexp}.FullVData = All.FullVData;
                    else
                        DATA.Clusterfile{nexp}.FullVData = [];
                    end
                    nprobes(nexp) = length(All.AllClusters);
                    All = rmfield(All,'AllClusters');
                    DATA.state.nospikes = 2;
                    DATA.state.somespikes = 1;
                elseif isfield(All,'Spikes')  && length(All.Spikes.times) > 1
                    Spikes = All.Spikes;
                else
                    Spikes = Spikes;
                    if isfield(trls.Probes,'probes')
                        nprobes(nexp) = max([trls.Probes.probe]);
                    end
                    if isfield(DATA,'AllClusters') %Will be loading from clustertimes
                        DATA.state.nospikes = 2;
                        DATA.state.somespikes = 1;
                    end
                end
                
                DATA.defaults.starttrial = 1+ Expts{nexp}.Trials(end).Trial;
                if isfield(trls,'Probes');
                    for j = 1:length(trls.Probes)
                        if isfield(trls.Probes,'var')
                            probenames{trls.Probes(j).probe} = trls.Probes(j).var;
                        end
                        Probes(trls.Probes(j).probe) = trls.Probes(j);
                    end
                    if DATA.state.autocutatstart && isfield(All,'Spike2Ch')
                        All.toplevel = DATA.toplevel;
                        DATA.AllData.toplevel = DATA.toplevel;
                        for j = 1:length(trls.Probes)
                            if trls.Probes(j).source == All.Spike2Ch.source
                                [DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, All.Spike2Ch.filename , trls.Probes(j).var,trls.Probes(j).probe);
                                DATA.AllData.Spikes.exptid = nexp;
                                DATA.spklist = 1:length(DATA.AllData.Spikes.times);
                                DATA.probe = trls.Probes(j).probe;
                                DATA = CalcClusterVars(DATA,DATA.spklist,'expt',nexp);
                                if DATA.state.autocutatstart == 2
                                    DATA = cmb.AutoCut(DATA,nexp);
                                end
                            end
                        end
                    end

                end
                nexp = nexp+1;
            else
                msgbox(sprintf('%s Incomplete Expt',d(j).name));
            end
        end
    end
end
Expts = SetTrialOffsets(Expts);
end %end if method 
fprintf('Loading took %.1f (%.1f)\n',mytoc(ts),sumload);

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
        DATA = cmb.AddMultiProbeGUI(DATA);
        DATA.probesource = cmb.FindProbeSources(DATA);
    else
        DATA.probesource = 1;
    end
elseif DATA.bysuffix  %Get here for Spike2 Array DATA. with or without quicksuffix
    if exist('nprobes','var')
        DATA.probelist = 1:max(nprobes);
    elseif isfield(DATA.ArrayConfig, 'id')
        DATA.probelist = 1:length(DATA.ArrayConfig.id);
    else
        DATA.probelist = 1;
    end
    if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
        DATA.probe = DATA.probelist(1);
    end
    for j = DATA.probelist
        DATA.probes(j).probe = j;
    end
    %            DATA.probes = Probes;
    if length(DATA.probelist) > 1
        DATA = cmb.AddMultiProbeGUI(DATA);
        DATA.probesource = cmb.FindProbeSources(DATA);
    else
        DATA.probesource = 1;
    end
    DATA.probenames = cellstr(int2str(DATA.probelist'));
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
    if DATA.state.usensx
        DATA.state.somespikes = 2;
    else
        DATA.state.somespikes = 1;
    end
    DATA.state.nospikes = 2;
    
    for j = 1:np
        DATA.probes(j).probe = j;
    end
    if sum(DATA.probelist == DATA.probe) == 0 % current probe not in list
        DATA.probe = DATA.probelist(1);
    end
    
    if isfield(DATA.ArrayConfig,'id')
        for j = 1:length(DATA.ArrayConfig.id)
            DATA.probenames{j} = sprintf('%d(E%d %d,%d)',j,DATA.ArrayConfig.id(j),DATA.ArrayConfig.X(j),DATA.ArrayConfig.Y(j));
        end
    end
    if np > 1
        DATA = cmb.AddMultiProbeGUI(DATA);
    end
end
if sum(DATA.probelist < 100) > 1 % count # of recording chans
    DATA.state.includeprobename = 1;
end

cmb.SetProbeList(DATA);

if nexp == 1
    acknowledge(sprintf('No expts in %s',name),'combine error');
else
    DATA = cmb.cListExpts(DATA,Expts);
    DATA.Expts = Expts;
    if method == 2
        DATA.AllData.Spikes = [];
    else
        DATA.AllData = All;
        DATA.AllData.Spikes = Spikes;
        DATA.AllData.SpikeIdx = SpkId;
        DATA.AllData.Trialids = Trialids;
        DATA.AllData.TrialTimes = TrialTimes;
        DATA.AllData.pcs = [];
    end
    if DATA.state.autosetlist && DATA.state.online == 1 && length(DATA.probelist) <= 1
        for j = lastn:nexp-1
            DATA = SetExptSpikes(DATA, j, 0);
        end
    end
end
DATA.name = name;
if DATA.bysuffix && isfield(DATA,'AllClusters') && DATA.state.online == 0
    DATA.state.online = 0;
    if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'excludetrials')
        for j = 1:length(DATA.CellDetails.exptids)
            eid = floor(DATA.CellDetails.exptids(j));
            cid = find(DATA.exptnos == eid);
            if cid <= length(DATA.AllClusters)
                for k = 1:size(DATA.CellDetails.excludetrials,2)
                    for c = 1:size(DATA.CellDetails.excludetrials,3)
                        if ~isempty(DATA.CellDetails.excludetrials{j,k,c})
                            if k <= length(DATA.AllClusters{cid})
                                if c > length(DATA.AllClusters{cid}(k).excludetrialids)
                                    DATA.AllClusters{cid}(k).excludetrialids{c} = DATA.CellDetails.excludetrials{j,k,c}(:)';
                                else
                                    DATA.AllClusters{cid}(k).excludetrialids{c} = union(DATA.AllClusters{cid}(k).excludetrialids{c}(:)',DATA.CellDetails.excludetrials{j,k,c}(:)');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    for j = 1:length(DATA.Expts)
        if isfield(DATA,'Clusterfile') && length(DATA.Clusterfile) >= j
            if isfield(DATA.Clusterfile{j},'excludetrialids')
                [DATA.Clusterfile{j}.excludetrialids{1}] = deal([]);
            end
            if isfield(DATA.Clusterfile{j},'FullVData')
                xcl = cmb.FindMissingTrials(DATA, DATA.Clusterfile{j}.FullVData, j);
                if length(xcl)
                    for k = 1:length(DATA.AllClusters{j})
                        for c = 1:length(DATA.AllClusters{j}(k).excludetrialids)
                            DATA.AllClusters{j}(k).excludetrialids{c} =  [DATA.AllClusters{j}(k).excludetrialids{c}(:)' xcl(:)'];
                        end
                    end
                end
            end
        end
        AllExpts{j} = DATA.Expts{j};
        AllExpts{j} = rmfields(AllExpts{j},'Trials');
        
        
        AllExpts{j}.ntrials = length(DATA.Expts{j}.Trials);
        AllExpts{j}.ids = [DATA.Expts{j}.Trials.id];
    end
    explistfile = [DATA.name '/AllExpts.mat'];
    if ~exist(explistfile,'file') && ~isempty(AllExpts)
        save(explistfile,'AllExpts');
    end
    allxcl = [];
    try
        for j = 1:length(DATA.AllClusters)
            xcl = [];
            for k = 1:length(DATA.AllClusters{j})
                xid = cat(2,DATA.AllClusters{j}(k).excludetrialids{:});
                xcl = cat(1,xcl,xid(:));
            end
            xcl = unique(xcl);
            allxcl = [allxcl xcl(:)'];
            if ~isempty(xcl)
                fprintf('Ex%d %d exclusions%s\n',DATA.exptnos(j),length(xcl),sprintf(' %d',xcl));
            end
        end
    catch ME
        cprintf('red', 'Error checking exclusion list\n');
    end
    if isfield(DATA,'CellDetails') && ~isfield(DATA.CellDetails,'checkedtimes')
        try
            CheckExptClusters(DATA.AllClusters,DATA.Expts, 'fix');
        catch
            cprintf('red','Error Checking Cluster Strutures\n');
        end
    end
else
    DATA.state.online = 1;
end
fprintf('Directory Read took %.2f\n',mytoc(functs));
DATA.exabsid = 1:length(DATA.Expts);
DATA = cmb.CheckCellExptList(DATA);


