function DATA = LoadSpikes(DATA, eid, varargin)

nocheck = 0;
usecache = DATA.state.usexycache;
if ~isfield(DATA,'probes') || DATA.state.nospikes == 1
    return;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nocheck',5) % to avoid recursion
        nocheck = 1;
    elseif strncmpi(varargin{j},'force',5) % don't use cache
        usecache =0;
        nocheck = 1;
    end
    j = j+1;
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
loaded = 1;
onexpt = 0; %temp kludge to allow quicker probe switching set 1 1 = only load spikes for 1 expt
if DATA.state.online == 0
    if DATA.bysuffix
        eid = DATA.currentexpt(1); %makes a mess called from CountSpikes
        [d, a,b]  = fileparts(DATA.Expts{eid}.Header.loadname);
        name = regexprep(a,'\.([0-9,a]*)$',['.p' num2str(DATA.probe) 't$1.mat']);
        filename = [d '/Spikes/' name];
        [DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, filename, [],[DATA.probe DATA.subprobe]);
        if isfield(DATA.AllData.Spikes,'times')
            DATA.Expts{DATA.currentexpt(1)}.gui.spks = 1:length(DATA.AllData.Spikes.times);
        else
            DATA.Expts{DATA.currentexpt(1)}.gui.spks = [];
        end
    elseif strncmp(DATA.filetype,'GridData',8)
        if ~isfield(DATA,'AllClusters')
            DATA.AllClusters = {};
        end
        %                DATA.AllData.Spikes = LoadGridSpikes(DATA, DATA.probe, DATA.currentexpt(1));
        if nocheck == 0
            DATA = cmb.CheckClusterLoaded(DATA, DATA.currentexpt(1), DATA.probe);
        end
    elseif length(id) > 1 || length(DATA.probes) > 2 %KLUDGE Need to detect new files with only 1 expt
        if onexpt
            DATA.AllData.Spikes = cmb.GetProbeFiles(DATA,DATA.probe,DATA.subprobe,'trange',DATA.Expts{eid}.Header.trange);
        else
            DATA.AllData.Spikes = cmb.GetProbeFiles(DATA,DATA.probe,DATA.subprobe);
        end
    else
        if DATA.probelist(id) > 16
            filename = strrep(DATA.datafilename,'.mat',sprintf('A.p%s.mat',DATA.probevars{id}(3:end)));
        else
            filename = strrep(DATA.datafilename,'.mat',sprintf('.p%s.mat',DATA.probevars{id}(3:end)));
        end
        [d, a,b]  = fileparts(filename);
        filename = [d '/Spikes/' a b];
        [DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, filename, DATA.probevars{id},[DATA.probe DATA.subprobe]);
    end
    DATA = cmb.SetSpkLists(DATA);
    DATA.spklist = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
    DATA = SetExptSpikes(DATA, eid, 0);
elseif strncmp(DATA.filetype,'Grid',4) && DATA.state.usensx > 0 %GetNS% deeals with NEV and NS5
    DATA =  cmb.GetNS5Spikes(DATA, DATA.currentexpt(1),  DATA.probe);
    DATA = SetExptSpikes(DATA, eid, 0);
elseif strncmp(DATA.filetype,'Grid',4)
else
    if isfield(DATA.Expts{eid}.Header,'loadname')
        filename = DATA.Expts{eid}.Header.loadname;
    else
        filename = ['C:' DATA.Expts{eid}.Header.Name];
    end
    if DATA.probelist(id) > 16  || DATA.probesource(id) == 2
        filename = strrep(filename,'/Expt','A/Expt');
        filename = regexprep(filename,'(\.[0-9]*.mat)','A$1');
    end
    if ~exist(filename,'file')
        [a,b] = fileparts(DATA.Expts{eid}.Header.Name);
        filename = [DATA.datafilename '/'  b '.mat'];
    end
    if strncmp(DATA.filetype,'Grid',4) && ~isfield(DATA,'probevars')
        if DATA.state.usensx == 1
            DATA =  cmb.GetNS5Spikes(DATA, DATA.currentexpt(1),  DATA.probe);
        else
            loaded = 0;
        end
    else
        ok = 0;
        if usecache
            [DATA, ok] = cmb.SpkCache(DATA,eid,DATA.probe','add');
        end
        if ~ok %either didnt' try cache, or its bad
            ts = now;
            fprintf('Loading Spikes P%d from %s',DATA.probe,filename);
            DATA.AllData.toplevel = DATA.toplevel;
            [DATA.AllData.Spikes]= cmb.GetProbeSpikes(DATA.AllData, filename , DATA.probevars{id},[DATA.probe DATA.subprobe]);
            DATA.AllData.Spikes.exptid = eid;
            if isfield(DATA.AllData.Spikes,'loaddur')
                loaddur = DATA.AllData.Spikes.loaddur;
            else
                loaddur = NaN;
            end
            DATA.spklist = [];
            fprintf('took %.2f (%.2f loading)\n',mytoc(ts),loaddur);
            DATA.AllData.Spikes.loaddur = [mytoc(ts), loaddur];
        end
%typically, this is where previous cluster is applied        
        if nocheck
%but if call is from SetSpkCodes ('force'), then avoid recursion            
            DATA = SetExptSpikes(DATA, eid, 'setrange');
        else
            DATA = SetExptSpikes(DATA, eid, 0);
        end
    end
    if loaded %otherwise recursively calls
        DATA = cmb.CountSpikes(DATA,eid);
    end
    if ~isfield(DATA.AllData.Spikes,'exptid')
        DATA.AllData.Spikes.exptid = eid;
    end
    %  DATA = SetSpkCodes(DATA, [1:length(DATA.AllData.Spikes.times)], 0)
end





