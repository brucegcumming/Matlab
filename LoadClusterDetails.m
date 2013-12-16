function [ClusterDetails, D] = LoadClusterDetails(name, varargin)
%[ClusterDetails, details] = LoadClusterDetails(name)
%load ClusterDetails files matching a named ClusterTimes file
%will get relevant data from both ClusterTimesDetails and Auto
%Converts any ints to doubles for newer versions

getauto = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noauto',5)
        getauto = -1;
    end
    j = j+1;
end
dname = name;
daname = name;
    if isempty(strfind(name,'ClusterTimesDetails'))
        dname = strrep(name,'ClusterTimes', 'ClusterTimesDetails');
    end
    if isempty(strfind(name,'AutoClusterTimes'))
        daname = strrep(name,'ClusterTimes', 'AutoClusterTimesDetails');
    else
        daname = strrep(name,'ClusterTimes', 'ClusterTimesDetails');
    end
    if exist(daname,'file') && (getauto >= 0 || ~exist(dname,'file'))
        tic;
        load(daname);
        D.loaddur(2) = toc;
        AutoClusterDetails = ClusterDetails;
    else
        AutoClusterDetails = [];
    end
    if exist(dname,'file')
        tic;
        load(dname);
        D.loaddur(1) = toc;
        for k = 1:length(AutoClusterDetails)
            if k > length(ClusterDetails) || isempty(ClusterDetails{k})
                ClusterDetails{k} = AutoClusterDetails{k};
            end
        end
    else
        ClusterDetails = AutoClusterDetails;
    end
tic;
for j = 1:length(ClusterDetails)
    if isfield(ClusterDetails{j},'maxv') && isinteger(ClusterDetails{j}.xy)
        ClusterDetails{j} = CondenseDetails(ClusterDetails{j},'reverse');
    end
end
D.loaddur(3) = toc;
    
