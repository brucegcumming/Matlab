function res = iscluster(C,c,p)
%res = iscluster(Cluster,clid,probe) is Cluster(clid,probe) defined /set
%res = 1 = defined and set. 2 = defined and not set or only auto set
if length(p) > 1
    for j = 1:length(p)
        res(j) = cmb.iscluster(C,c,p(j));
    end
    return;
end
if isfield(C,'Header') 
    if ~isfield(C,'Cluster') %Expt with no clusters defined
    res = 0;
    return;
    else
        C = C.Cluster;
    end
end
if ~iscell(C) | isempty(C) | size(C,1) < c | size(C,2) < p | isempty(C{c,p})
    res = 0;
elseif ~isfield(C{c,p},'params') | ~isfield(C{c,p},'x')
    res = 0;
elseif ~isfield(C{c,p},'nspk') || C{c,p}.nspk == 0
    res = 2;
else
    res = 1;
end

