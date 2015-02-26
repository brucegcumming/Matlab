function [ispk, spktimes, codes] = FindSpikes(DATA, times, probe, range)
% [ispk, spktimes, codes] = FindSpikes(DATA, times, probe, range)
% for Combine.  Finds Spikes that fall between times(1) and times(2)

if DATA.plot.showartifacts
    maxcl = 9;
else
    maxcl = 8;  %>8 is artifacts. Set this to 9 to include artifacts
end
if isfield(DATA,'AllSpikes')
    if ~isfield(DATA.AllSpikes{probe},'codes') | isempty(DATA.AllSpikes{probe}) | ~isfield(DATA.AllSpikes{probe},'times')
        ispk = [];
        return;
    end
    ispk = find(DATA.AllSpikes{probe}.times > times(1) &...
    DATA.AllSpikes{probe}.times < times(2) & ...
    DATA.AllSpikes{probe}.codes(:,2) < maxcl);
    if nargout > 1
    spktimes = DATA.AllSpikes{probe}.times(ispk);
    codes = DATA.AllSpikes{probe}.codes(ispk,2);
    end
elseif isfield(DATA,'AllClusters') && iscell(DATA.AllClusters)
    if length(range)
        nexp = range(1);
    else
        nexp = DATA.currentexpt;
    end
    if length(DATA.AllClusters) < min(nexp) || length(DATA.AllClusters{nexp(1)}) < probe
        ispk = [];
    else
        ispk = [];
        for j = 1:length(nexp)
            ispks{j} = find(DATA.AllClusters{nexp(j)}(probe).times > times(1) &...
                DATA.AllClusters{nexp(j)}(probe).times < times(2));
        end
        ispk = cat(1,ispks{:});
    end
elseif isfield(DATA,'AllClusters')
    ispk = find(DATA.AllClusters(probe).times > times(1) &...
    DATA.AllClusters(probe).times < times(2));
elseif isempty(DATA.AllData.Spikes) || ~isfield(DATA.AllData.Spikes,'times')...
        || isempty(DATA.AllData.Spikes.times)
    ispk = [];
else
    if isempty(range)
        ispk = find(DATA.AllData.Spikes.times > times(1) &...
            DATA.AllData.Spikes.times < times(2) ...
            & DATA.AllData.Spikes.codes(:,2) < maxcl);
    else
        ispk = find(DATA.AllData.Spikes.times(range) > times(1) &...
            DATA.AllData.Spikes.times(range) < times(2)...
            & DATA.AllData.Spikes.codes(range,2) < maxcl);
        ispk = range(ispk);
    end
end
if size(ispk,1) == 1
    ipsk = ispk';
end