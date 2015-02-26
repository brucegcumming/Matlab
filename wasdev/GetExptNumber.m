function n = GetExptNumber(name, varargin)
%n = GetExptNumber extracts the number of an experiment from a filename
%
n = 0;


if isfield(name,'Header')
    E = name;
    if isfield(E.Header,'exptno')
        n = E.Header.exptno;
    elseif isfield(E.Header,'Name')
        n = GetExptNumber(E.Header.Name);
    else
        n = NaN;
    end
    return;
end
id = regexp(name, 'Expt[0-9,a]*[ACF][lu].*.mat');
if isempty(id)
    id = regexp(name, 'Expt[0-9]*.p[0-9]*FullV.*.mat');
end    
xid = regexp(name, 'Expt[0-9]*a[ACF][lu].*.mat');
if length(id) == 1
    n = sscanf(name(id+4:end),'%d');
    if length(xid)
        n = n+0.1;
    end
else
    id = regexp(name, '\.[0-9]*.mat');
    if isempty(id)
        id = regexp(name, '\.[0-9]*.lfp.mat');
    end
    if ~isempty(id)
        n = sscanf(name(id+1:end),'%d');
    end
end
