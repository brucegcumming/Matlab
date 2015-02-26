function DATA = ReClassifyAll(DATA, varargin)
% apply defind clusters to all expts/probes
% useful if saved list of classification codes is missing/ wrong
% need to run this when loading all spikes (because, err....)
mkmean = 0;
j = 1;
probes = DATA.probelist;
while j <= length(varargin)
    if strncmpi(varargin{j},'mkmean',5)
        mkmean = 1;
    elseif strncmpi(varargin{j},'probes',5)
        j = j+1;
        probes = varargin{j};
    end
    j = j+1;
end
if 0
    [DATA, ispk] = SetExptSpikes(DATA, DATA.currentexpt(1),'setrange');
    DATA = SetSpkCodes(DATA, ispk, DATA.probe, 0);
end
oldprobe = DATA.probe;
ts = now;
if length(probes) == 1 && DATA.probe == probes(1)
    DATA = CalcClusterVars(DATA, 1:length(DATA.AllData.Spikes.times));
end
for k = 1:length(probes)
    if ~isfield(DATA,'AllSpikes')
        if DATA.probe ~= probes(k) || length(probes) > 1
            DATA = cmb.SetProbe(DATA, probes(k));
            DATA = CalcClusterVars(DATA, 1:length(DATA.AllData.Spikes.times));
        end
    elseif k <= length(DATA.AllSpikes) && ~isempty(DATA.AllSpikes{k})
        DATA.probe = k;
        DATA = CalcClusterVars(DATA, 1:length(DATA.AllSpikes{k}.times));
    end
    if DATA.probe == probes(k) % have spikes for this
        ts = [ts now];
        for j = 1:length(DATA.Expts)
            if isfield(DATA.Expts{j},'Cluster')
                DATA.cluster = DATA.Expts{j}.Cluster;
                if cmb.iscluster(DATA.cluster, 1, k) == 1
                    [DATA, ispk] = SetExptSpikes(DATA, j,'setrange');
                    tic;
                    DATA = SetSpkCodes(DATA, ispk, DATA.probe, 0);
                end
            end
            ts = [ts now];
            ctimes(k,j) = ts(end)-ts(end-1);
        end
        if mkmean
            DATA = cmb.CalcMeanSpike(DATA,1:length(DATA.Expts));
        end
    end
end
DATA.probe = oldprobe;


