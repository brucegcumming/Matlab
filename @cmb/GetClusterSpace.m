function [x,y, DATA] = GetClusterSpace(DATA, Expt)


x = DATA.plot.clusterX;
y = DATA.plot.clusterY;
p = DATA.probe;
if isfield(Expt,'Cluster') & ~isempty(Expt.Cluster) & size(Expt.Cluster,2) >= p
j = 1;
while j < size(Expt.Cluster,1) && isempty(Expt.Cluster{j,p})
j = j+1;
end
if j <= size(Expt.Cluster,1) & p <= size(Expt.Cluster,2) & ~isempty(Expt.Cluster{j,p}) & isfield(Expt.Cluster{j,p},'params')
x = Expt.Cluster{j,p}.params(1);
y = Expt.Cluster{j,p}.params(2);
DATA.clusterArange = Expt.Cluster{j,p}.Arange;
DATA.clusterBrange = Expt.Cluster{j,p}.Brange;
if isfield(Expt.Cluster{j,p},'Erange')
DATA.clusterErange = Expt.Cluster{j,p}.Erange;
end
end
else
x = DATA.plot.clusterX;
y = DATA.plot.clusterY;
end




