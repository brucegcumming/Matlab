function out = SpkCache(DATA,e,p,action, varargin)
%Store times, classification, and Cluster xy for online data
%so that don't need to reload spikes

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',4)
        out = CheckCache(DATA,e,p);
        return;
    end
    j = j+1;
end

if strcmp(action,'check')
    out = CheckCache(DATA,e,p);
    return;
end

if strcmp(action,'set')
    DATA.spkcache{e,p}.codes = DATA.AllData.Spikes.codes(:,2);
    DATA.spkcache{e,p}.times = DATA.AllData.Spikes.times;
    DATA.spkcache{e,p}.ncut = sum(DATA.spkcache{e,p}.codes > 0);
    DATA.spkcache{e,p}.cx = DATA.Spikes.cx;
    DATA.spkcache{e,p}.cx = DATA.Spikes.cy;
    out = DATA;
elseif strcmp(action,'get')
    ok = CheckCache(DATA,e,p);
    if ok
        S = DATA.spkcache{e,p};
        DATA.AllData.Spikes.times = S.times;
        DATA.AllData.Spikes.codes(:,2) = codes;
        DATA.Spikes.cx = S.cx;
        DATA.Spikes.cy = S.cy;
    end
end

function ok = CheckCache(DATA,e,p)

ok = 0;
sz = size(DATA.spkcache);
if sz(1) >=e && sz(2) >= p 
    if isfield(DATA.spkcace{e,p},'ncut')
        if isfield(DATA.spkcace{e,p}.ncut > 1
            ok = 1;
        end
    elseif isfield(DATA.spkcace{e,p},'codes')
        if sum(DATA.spkcache{e,p}.codes > 0) > 1
            ok = 1;
        end
    end
end
