function [RefClusters, rfile] = LoadRefClusters(DATA, varargin)

RefClusters = {};
currentset = 0;

if isdir(DATA.name)
    rfile = [DATA.name '/RefClusters.mat'];
else
    rfile = [fileparts(DATA.name) '/RefClusters.mat'];    
end
if isfield(DATA,'refcut') && isfield(DATA.refcut,'currentset') && DATA.refcut.currentset > 0
    rfile = strrep(rfile,'RefClusters.mat',['Ref' int2str(DATA.refcut.currentset) 'Clusters.mat']);
    currentset = DATA.refcut.currentset;
end

if exist(rfile,'file')
    load(rfile);
    RefClusters = FixCluster(Clusters);
    for j = 1:length(RefClusters)
        RefClusters{j}.refcutset = currentset;
    end
else
    fprintf('Can''t load RefCluster File %s\n',rfile);
end