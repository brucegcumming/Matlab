function C = StripClusters(Clusters)        for j = 1:length(Clusters)    C{j} = Clusters{j};    if isfield(C{j},'xy')        C{j} = rmfield(C{j},'xy');    endend