function ClusterLog(name, varargin)
%ClusterLog(name,...) Show contents of ClusterLog
%ClusterLog.mat files are modified each time a save is made
%records space, cluster #, shape, ncut/nevents
%ClusterLog{j}.recluster(1) records whether AllVPcs was called with
%reclassify (2) or recluster (1). If any manual changes were made, it will
%be 0.
%ClusterLog{j}.recluster(2) records whether is an automatic cut.

setprobe = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'probe',5)
        j = j+1;
        setprobe = varargin{j};
    end
    j = j+1;
end

if exist(name)
    ClusterLog = {};
    load(name);
elseif iscell(name)
    ClusterLog = name;
end
probes = CellToMat(ClusterLog,'probe');
id = 1:length(ClusterLog);
if setprobe
    id = find(probes == setprobe);
end
for j = 1:length(id)
    
    C = ClusterLog{id(j)};
    if ~isfield(C,'quick')
        C.quick = NaN;
    end
    fprintf('%s %s E%dP%d %d %d\n',datestr(C.savetime),C.user,C.exptno,C.probe,C.recluster(1),C.quick);
end
