function nc = CountClusters(Cluster)
nc = 0;
if isempty(Cluster)
nc = 0;
else
for k = 1:size(Cluster,2) 
ncs(k) = 0;
for j = 1:min([size(Cluster,1) 7])  %Clusters > 7 are artifact
if isfield(Cluster{j,k},'x')
ncs(k) = j;
end
end
end
nc = max(ncs);
end

