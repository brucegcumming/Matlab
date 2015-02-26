function DATA = ClearExptClusters(DATA)

for j = 1:length(DATA.Expts)
DATA.Expts{j}.Cluster = {};
end


