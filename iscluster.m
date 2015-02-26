function res = iscluster(C,c,p)
%res = iscluster(C) determine if C is a cluster strcutre
%res = iscluster(C,c,p) true if cluster is defined for probe p, cluster c
% if C is  a cell array of cluster cell arrays
%res = iscluster(C,e,p) true if cluster is defined for expt e, cluster c
%
res = 0;
if nargin == 1
    if iscell(C) && iscell(C{end}) && isfield(C{end}{end},'mahal') %Cell array of cluster arrays
        res = 3;
    elseif iscell(C) && isfield(C{end},'mahal')
        res = 2;
    elseif isfield(C,'mahal')
        res = 1;
    end
elseif isempty(c)
    res = 0;
elseif isempty(C)
    res = 0;    
elseif iscell(C{1}) %Clusters{e}{p} like PlotClusters
    if c > length(C) || p > length(C{c})
        res = 0;
    elseif ~isfield(C{c}{p},'shape')
        res = 0;
    else
        res = 1;
    end
elseif iscell(C) && isempty(C{1})
    res = 0;
    for j = 1:length(C)
        if ~isempty(C{j})
            if iscell(C{j}) %Plotclusters styel C{ex}{p}
                if j == c && isfield(C{c}{p},'shape')
                    res = 1;
                    return;
                end
            end
        end
    end
else
    if ~iscell(C) | isempty(C) | size(C,1) < c | size(C,2) < p | isempty(C{c,p})
        res = 0;
    elseif ~isfield(C{c,p},'params') | ~isfield(C{c,p},'x')
        res = 0;
    elseif ~isfield(C{c,p},'nspk') || C{c,p}.nspk == 0
        res = 2;
    else
        res = 1;
    end
end
