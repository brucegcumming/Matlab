function clist = ClusterDef(C, varargin)
%returns list of defined clusters in C

if C.ncut > 0
    clist = 1;
else 
    clist = 0;
end

for j = 1:length(C.next)
    if ~isempty(C.next{j})
        clist(end+1) = j+1;
    end
end