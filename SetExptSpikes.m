function [DATA, ispk, dprime, details] = SetExptSpikes(DATA, expid, show, varargin)
%[DATA, ispk, dprime, details] = SetExptSpikes(DATA, expid, show, varargin)
%Find Which Spikes fall within a given Expt for combine, then calculate
%cluster paramters
%SetExptSpikes(DATA, expid, 'setrange')  just finds spk list

    details = [];
    dprime = 0;
    ispk = [];
    if  DATA.state.online
        quickmode = 1;
    else
        quickmode = 0;
    end
    
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
        elseif strncmpi(varargin{j},'quick',5)
            quickmode = 1;
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
    probe = GetProbe(DATA, expid, DATA.probe);
    if probe == 0 %Cell not defined here
        counts = [];
        return;
    end
    if isempty(ispk) & ~DATA.state.psychonly
        if isfield(DATA,'AllSpikes')
            espk = [];
        elseif isfield(DATA,'AllClusters')
            espk = [];
        elseif isfield(DATA.AllData.Spikes,'times')
            espk = ExptSpikeListAll(DATA, expid, DATA.AllData.Spikes.times);
        else
            espk = [];
        end
        
        if isempty(espk)
            ispk = [];
        elseif quickmode
                times(1) = DATA.Expts{expid}.Trials(1).Start(1)-DATA.state.preperiod-10000;
                times(2) = DATA.Expts{expid}.Trials(end).End(end)+DATA.state.postperiod+10000;
                if DATA.bysuffix || strncmp(DATA.filetype,'Grid',4)
                    espk = expid;
                end
                ispk = FindSpikes(DATA, times, probe, espk);
        else
            for j = 1:length(DATA.Expts{expid}.Trials)
                times(1) = DATA.Expts{expid}.Trials(j).Start(1)-DATA.state.preperiod;
                times(2) = DATA.Expts{expid}.Trials(j).End(end)+DATA.state.postperiod;
                if DATA.bysuffix || strncmp(DATA.filetype,'Grid',4)
                    espk = expid;
                end
                tspk = FindSpikes(DATA, times, probe, espk);
                ispk = [ispk; tspk];
            end
        end
    end
    
DATA.Expts{expid}.gui.spkrange = [min(ispk) max(ispk)];
DATA.Expts{expid}.gui.spks = ispk;


%Don't reset codes if using Clusters.clst
if (ischar(show) & strcmp(show,'setrange')) || DATA.state.nospikes  %if 
    details.nc = 0;
    return;
end

if ~isfield(DATA,'spklist') | isempty(DATA.spklist)
    DATA.spklist = ispk;
end

DATA = CalcClusterVars(DATA,  ispk, 'expt',expid);
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
