function Cluster = SmallCluster(AllClusters, e, p)

Cluster.probe = p;
if length(AllClusters) < e || length(AllClusters{e}) < p
Cluster.mahal = NaN;
return;
end
C = AllClusters{e}(p);

f = {'mahal' 'crit' 'ctime' 'excludetrialids' 'missingtrials' 'clustersign' 'cluster' 'dips' 'dropi' 'savetime' 'user'};
for j = 1:length(f)
if isfield(C,f{j})
Cluster.(f{j}) = C.(f{j});
end
end

