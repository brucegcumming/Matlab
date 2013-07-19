function ClusterDetails = LoadCluterDetails(name)
%ClusterDetails = LoadCluterDetails(name)
%load ClusterDetails files matching a named ClusterTimes file
%will get relevant data from both ClusterTimesDetails and Auto

    dname = strrep(name,'ClusterTimes', 'ClusterTimesDetails');
    daname = strrep(name,'ClusterTimes', 'AutoClusterTimesDetails');
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
    else
        ClusterDetails = AutoClusterDetails;
    end
