function AutoClusters = SaveAutoClusters(DATA)

cfile = cmb.ClusterFile(DATA,'auto');
if exist(cfile,'file')
load(cfile);
end
for j = 1:length(DATA.Expts)
for k = 1:length(DATA.probelist)
p = DATA.probelist(k);
if isfield(DATA.Expts{j}.Cluster{1,p},'autocut') && DATA.Expts{j}.Cluster{1,p}.autocut
AutoClusters{j}{1,p} = DATA.Expts{j}.Cluster{1,p};
end
end
end
save(cfile,'AutoClusters');


