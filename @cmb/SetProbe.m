function [DATA, ok] = SetProbe(DATA, probe)


DATA.probe = probe;
if DATA.listbycell == 2
    str = sprintf('Loading probe(Mu) %d...',DATA.probe);
elseif DATA.listbycell == 1
    str = sprintf('Loading Cell %d...',DATA.probe);
else
    str = sprintf('Loading Probe %d...',DATA.probe);
end
if isfield(DATA,'toplevel')
    set(DATA.toplevel,'Name',str);
end
drawnow;
fprintf('%s',str);
tstart = now;
if ~isfield(DATA,'currentexpt')
    DATA.currentexpt = 1;
end
if isfield(DATA,'AllClusters')
    ts = now;
    probe = GetProbe(DATA, DATA.currentexpt(1), DATA.probe);
    if probe > 0
        [DATA, D] = cmb.CheckClusterLoaded(DATA, DATA.currentexpt(1), probe);
        if sum(D.loaddur) > 0
            fprintf('(Loading%s)',sprintf(' %.2f',D.loaddur));
            DATA.newload = 1;
        end
        if iscell(DATA.AllClusters)
            C  = DATA.AllClusters{DATA.currentexpt(1)}(probe);
            DATA.Spikes.cx = C.cx;
            DATA.Spikes.cy = C.cy;
        else
            DATA.Spikes.cx = DATA.AllClusters(probe).cx;
            DATA.Spikes.cy = DATA.AllClusters(probe).cy;
        end
    end
    %    DATA.AllData.Spikes.times = DATA.AllClusters(probe).times;
    %    DATA.AllData.Spikes.codes = zeros(length(DATA.AllClusters(probe).times),4);
    DATA.currentprobe = DATA.probe;
    if probe ~= DATA.oldprobe
        for j = 1:length(DATA.Expts)
            DATA.Expts{j}.gui.classified = 0;
            DATA.Expts{j}.gui.counted = 0;
        end
    end
    DATA.oldprobe = DATA.probe;
    cmb.NotBusy(DATA);
    mytoc(ts);
    
    return;
elseif isfield(DATA,'AllSpikes')
    for j = 1:length(DATA.expid)
        ts(j,:) = DATA.Expts{DATA.expid(j)}.Header.trange;
    end
    times = [min(ts(:,1)) max(ts(:,2))];
    DATA.AllSpikes{probe} = cmb.GetProbeFiles(DATA, probe,DATA.subprobe,'trange',times/10000,'nodv');
elseif sum([DATA.probes.probe] < 100) == 1 %ustim expt with just 1 channel
    DATA.probe = probe;
else
    DATA = cmb.LoadSpikes(DATA, DATA.currentexpt(1));
    DATA.loaddur = DATA.AllData.Spikes.loaddur;
    if DATA.state.nospikes && isfield(DATA,'AllClusters')
        DATA = cmb.CheckClusterLoaded(DATA, DATA.currentexpt(1), DATA.probe);
    else
        DATA = cmb.LoadClusters(DATA,cmb.ClusterFile(DATA),'noclear');
    end
end
DATA.currentcluster = 1; %else can get blank clusters ->error
for j = 1:length(DATA.Expts)
    [DATA, DATA.Expts{j}.gui.spks] = SetExptSpikes(DATA,j,'setrange');
    if DATA.state.online
        DATA.Expts{j}.gui.counted = 0;
    end
end
DATA.spklist = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
mytoc(tstart);
if isfield(DATA,'toplevel') & DATA.state.autoplotnewprobe & length(DATA.probelist) > 2 & DATA.state.online == 0 && DATA.state.nospikes == 0
    cmb.PlotAllExpts(DATA);
end
cmb.NotBusy(DATA);

