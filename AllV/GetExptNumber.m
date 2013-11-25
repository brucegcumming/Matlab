function n = GetExptNumber(name, varargin)
%n = GetExptNumber extracts the number of an experiment from a filename
%
n = 0;

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
end
