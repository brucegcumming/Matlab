function res = isExptCluster(E,c,p)
if ~isfield(E,'Cluster');
res = 0;
else
res = cmb.iscluster(E.Cluster,c,p);
end


