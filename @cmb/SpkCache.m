function [out, varargout] = SpkCache(DATA,e,p,action, varargin)
%SpkCache(DATA,e,p,'set') Store times, classification, and Cluster xy for online data
%                   'add') copies cached data into DATA.AllData.Spikes
%                   'addxy' also checks that saved cluster space matches current
%                   'get') returns cache data for this e,p
%                   'check') returns 0 if cache not available,1 if it is.
%                            varargout{1} is the cache data
%                   'init') makes sure size is at least e x p
%                
%so that don't need to reload spikes
%

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',4)
        out = CheckCache(DATA,e,p);
        return;
    end
    j = j+1;
end


if strcmp(action,'check') %0 if no spikes classified, 1 if data and is classified 
    %2 means data but no events so may consider reloading
    out = CheckCache(DATA,e,p);
    if out
        varargout{1} = DATA.spkcache{e,p};
    else
        varargout{1} = [];
    end
    return;
end
if strcmp(action,'init') %make sure array  is initialized
    if(CheckCache(DATA,e,p) ==0)
        DATA.spkcache{e,p} = [];
        out = DATA;
        return;
    end
end
       
    

if strcmp(action,'set')
    nspk = length(DATA.AllData.Spikes.times);
    if size(DATA.AllData.Spikes.codes,1) > 1
        DATA.spkcache{e,p}.codes = DATA.AllData.Spikes.codes(:,2);
    else
        DATA.spkcache{e,p}.codes = DATA.AllData.Spikes.codes;
    end
    DATA.spkcache{e,p}.times = DATA.AllData.Spikes.times;
    DATA.spkcache{e,p}.ncut = sum(DATA.spkcache{e,p}.codes > 0);
    
%cx and cy can be longer than the current expt. For SU files
%this keep track of all expts 
%Can also be shorter if online and expt has spk events after expt end
    if nspk > length(DATA.Spikes.cx) 
        fprintf('Only %d Cluster Vals. (%d spikes) for Cache\n',length(DATA.Spikes.cx),nspk);
        nspk = length(DATA.Spikes.cx);
    end
    DATA.spkcache{e,p}.cx = DATA.Spikes.cx(1:nspk);
    DATA.spkcache{e,p}.cy = DATA.Spikes.cy(1:nspk);
    DATA.spkcache{e,p}.espk = minmax(DATA.spklist);
    DATA.spkcache{e,p}.space = DATA.Spikes.space;
    if isempty(DATA.spklist) %no spikes  to cut
        DATA.spkcache{e,p}.ncut = NaN;
    end
    out = DATA;
elseif strncmp(action,'add',3)
    expts = e;
    for e = expts(:)';
        ok(e) = CheckCache(DATA,e,p);
        if ok(e) && strcmp(action,'addxy')
            S = DATA.spkcache{e,p};
            if ~strcmp(S.space{1},DATA.Spikes.space{1}) || ~strcmp(S.space{1},DATA.Spikes.space{1}) 
                fprintf('E%dP%d Cluster Space Changed since XY cached\n',e,p);
                ok(e) = 0;
            end
        end
        if ok(e)
            S = DATA.spkcache{e,p};
            DATA.AllData.Spikes.times = S.times;
            DATA.AllData.Spikes.codes = [];
            if ~isempty(S.codes)
                DATA.AllData.Spikes.codes(:,2) = S.codes;
            else
                DATA.AllData.Spikes.codes = [];
            end
            DATA.AllData.Spikes.cache = 1;
            DATA.Spikes.cx = S.cx;
            DATA.Spikes.cy = S.cy;
        end
        varargout{1} = sum(ok) > 0;
    end
elseif strcmp(action,'get')
    ok = CheckCache(DATA,e,p);
    if ok
        out = DATA.spkcache{e,p};
    else
        out = [];
    end
    varargout{1} = ok;
    return;
end
out = DATA;

function ok = CheckCache(DATA,e,p)


ok = 0;
if DATA.state.online == 0
    return;
end
if ~isfield(DATA,'spkcache')
    return;
end
sz = size(DATA.spkcache);
if sz(1) >=e && sz(2) >= p 
    if isfield(DATA.spkcache{e,p},'ncut')
        if DATA.spkcache{e,p}.ncut > 1
            ok = 1;
        elseif isnan(DATA.spkcache{e,p}.ncut)
            ok = 2;
        end
    elseif isfield(DATA.spkcache{e,p},'codes')
        if sum(DATA.spkcache{e,p}.codes > 0) > 1
            ok = 1;
        elseif ~isempty(DATA.spkcache{e,p}.codes)
            ok = 2;
        end
    end
end
