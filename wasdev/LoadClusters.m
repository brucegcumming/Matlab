function [AllClusters, AllFullVData] = LoadClusters(dirname, eid, varargin)
%[AllClusters, AllFullVData] = LoadClusters(dirname, eid, varargin)
%OBSOLETE
%Now use LoadCluster()

AllClusters = {};

for j = 1:length(eid)
    loadname = [];
    name = [dirname '/Expt' num2str(eid(j)) 'ClusterTimes.mat'];
    dname = [dirname '/Expt' num2str(eid(j)) 'ClusterTimesDetails.mat'];
    daname = [dirname '/Expt' num2str(eid(j)) 'AutoClusterTimesDetails.mat'];
    aname = [dirname '/Expt' num2str(eid(j)) 'AutoClusterTimes.mat'];
    if exist(aname,'file')
        load(aname);
        AutoClusters = Clusters;
        AllFullVData{j} = FullVData;
    end
    if exist(name,'file')
        load(name);
        loadname = name;
        AllClusters{j} = Clusters;
        AllFullVData{j} = FullVData;
    else
        AllClusters{j} = AutoClusters;
        fprintf('Can''t read %s\n',name);
    end
    if exist(daname,'file')
        load(daname);
        AutoClusterDetails = ClusterDetails;
    else
        AutoClusterDetails = [];
    end
    if exist(dname,'file')
        load(dname);
        for k = 1:length(AutoClusterDetails)
            if k > length(ClusterDetails) || isempty(ClusterDetails{k})
                ClusterDetails{k} = AutoClusterDetails{k};
            end
        end
        f = {'Evec' 'clst' 't'};
        for k = 1:length(ClusterDetails)
            for n = 1:length(f)
                if isfield(ClusterDetails{k},f{n})
                    AllClusters{j}{k}.(f{n}) = ClusterDetails{k}.(f{n});
                end
            end
        end
        for k = 1:length(AllClusters{j})
            AllClusters{j}{k}.loadname = loadname;
        end
    end
end

if length(AllClusters) == 1
    AllClusters = AllClusters{1};
    AllFullVData = AllFullVData{1};
end