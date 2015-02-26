function DATA = UpdateAllClusters(DATA)
ts = now;
for j = 1:length(DATA.AllClusters)
cfile = sprintf('%s/Expt%dClusterTimes.mat',DATA.datafilename,j);
d = dir(cfile);
if length(d) == 1 && d.datenum > DATA.Clusterfile{j}.datenum
fprintf('ReReading %s\n',cfile);
load(cfile);
for k = 1:length(Clusters)
if DATA.Clusterfile{j}.quick
DATA.AllClusters{j}(k).times = round(Clusters{k}.times' .* 10000);
DATA.AllClusters{j}(k).codes = ones(length(Clusters{k}.times),2);
end
end
end
end
mytoc(ts);

