function DATA = LoadOnlineClusters(DATA, cfile)
if length(DATA.probelist) > 4
minq = 5;
else
minq = 0;
end
if exist(cfile,'file')
load(cfile);
DATA.state.recut = 1;
p = DATA.probe;
for j = 1:length(DATA.Expts)
eid(j,1) = DATA.Expts{j}.Trials(1).id;
eid(j,2) = DATA.Expts{j}.Trials(end).id;
end
for j = 1:length(AllClusters);
if ~isempty(AllClusters{j})
if isfield(AllClusters{j},'ids')
exid = find(eid(:,1) < AllClusters{j}.ids(2) & ...
eid(:,2) > AllClusters{j}.ids(1));
if length(exid) == 1 & isfield(AllClusters{j},'Cluster')
k = exid;
DATA.Expts{k}.OnlineCluster = AllClusters{j}.Cluster;
for ic = 1:length(AllClusters{j}.Cluster(:))
%only copy online cuts that were not automatic
%Feb 20011 Change to apply test to the quality assigned online
%don't copy assigned quality - it needs resetting when offline clsuteres
%are cut
C= AllClusters{j}.Cluster{ic};
if ~isfield(C,'autocut') || C.autocut == 0 
if (isfield(C,'quality') && C.quality >= minq) || minq == 0
C.quality = 0;
DATA.Expts{k}.Cluster{ic} = C;
DATA.Expts{k}.Cluster{ic}.touched = -1;
else
DATA.Expts{k}.OnlineCluster{ic}.touched = -1;
end
end
DATA.Expts{k}.OnlineCluster{ic}.touched = -1;
end
DATA.Expts{k}.gui.clustertype = 2;
if ~isfield(DATA.Expts{k},'Cluster')
DATA.Expts{k}.Cluster = {};
end
%
%Fill in fields that might be missing from old cluster files
%remove firstpk, lastspk refs from the online file - these spike id #s will
%not necessarily match the saved fil
for m = 1:size(DATA.Expts{k}.Cluster,1)
DATA.Expts{k}.Cluster{m,p}.firstspk = NaN;
DATA.Expts{k}.Cluster{m,p}.lastspk = NaN;
if ~isfield(DATA.Expts{k}.Cluster{m,p},'params')
DATA.Expts{k}.Cluster{m,p}.params = [1 2];
end
end
DATA.Expts{k}.gui.classified = 0; %loaded def, but not classified
end
end
end
end
end


