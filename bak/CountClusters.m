function [isset, details] = CountClusterList(C, varargin)


for e = 1:length(C)
    if isfield(C{e},'Cluster') & ~isempty(C{e}.Cluster)
    for j = 1:size(C{e}.Cluster,2);
        if isfield(C{e}.Cluster{1,j},'x')
            isset(e,j) = 1;
            if isfield(C{e}.Cluster{1,j},'autocut') && C{e}.Cluster{1,j}.autocut > 0
                isset(e,j) = 2;
            end
            if isfield(C{e}.Cluster{1,j},'nspk')
            nspks(e,j) = C{e}.Cluster{1,j}.nspk;
            else
            nspks(e,j) = NaN;
            end
        else
            isset(e,j) = 0;
        end            
        if isfield(C{e}.Cluster{1,j},'dprime')
            details.dprime(e,j) = C{e}.Cluster{1,j}.dprime;
        end
    end
    else
        isset(e,:) = 0;
    end
end