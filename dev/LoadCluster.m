function [AllClusters, AllFullVData, details] = LoadCluster(dirname, eid, varargin)
%LoadCluster(dirname, expts, ...)
%Load Cluster info for an expt or list of expts
%if length(expts) = 1 returns a 1 x nprobes cell matrix, else 
%a 1 x length(expts) matrix each containing 1 x nprobes cells
%
%LoadCluster Combines info from Clusters and ClusterDetails files
%so that clst, t and Evec fields are in clusters
%LoadCluster(dirname, expts, 'gextxy') also inlcudes xy
%LoadCluster(dirname, expts, 'rawxy') also inlcudes xy with any roation
%removed (PlotClusters needs it this way
%LoadCluster(dirname, expts, 'alltimes') 
%replaces Clusters{}.times with ClusterDetails{}.t, saving memory
%(Also needed by PlotClusters)

AllClusters = {};
AllFullVData = {};
f = {'Evec' 'clst' 't'}; %fields to copy
rawxy = 0;
alltimes = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'getxy',5)
        f = {f{:} 'xy' 'triggerV'};
    elseif strncmpi(varargin{j},'rawxy',5)
        rawxy = 1;
        f = {f{:} 'xy' 'triggerV'};
    elseif strncmpi(varargin{j},'alltimes',5)
        alltimes = 1;
    end
    j = j+1;
end
ts = now;
details = {};

if ispc && dirname(1) == '/'
    dirname(1) = '\';
end
for j = 1:length(eid)
    e = floor(eid(j));
    if rem(eid(j),e) == 0.1
        xs = 'a';
    else
        xs = '';
    end
    name = [dirname '/Expt' num2str(e) xs 'ClusterTimes.mat'];
    dname = [dirname '/Expt' num2str(e) xs 'ClusterTimesDetails.mat'];
    daname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimesDetails.mat'];
    aname = [dirname '/Expt' num2str(e) xs 'AutoClusterTimes.mat'];
    if exist(aname,'file')
        load(aname);
        AutoClusters = Clusters;
        for p = 1:length(AutoClusters)
            AutoClusters{p}.auto = 1;
        end
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        details{j}.loadname = aname;
    else 
        AutoClusters = {};
    end
    if exist(name,'file')
        load(name);
        details{j}.loadname = name;
        AllClusters{j} = Clusters;
        if exist('FullVData','var')
            AllFullVData{j} = FullVData;
        end
        for k = 1:length(AutoClusters)
            if k > length(Clusters) || isempty(Clusters{k})
                AllClusters{j}{k} = AutoClusters{k};
            end
        end
    elseif ~isempty(AutoClusters)
        AllClusters{j} = AutoClusters;
        fprintf('Can''t read %s\n',name);
    else
        fprintf('Can''t read %s or %s\n',name,aname);
        Clusters = {};
        return;
    end
    details{j}.loadtime = now;
    ClusterDetails = LoadClusterDetails(name);

    for k = 1:length(ClusterDetails)
        if isfield(ClusterDetails{k},'ctime') && ...
                isfield(AllClusters{j}{k},'savetime') && ...
                ClusterDetails{k}.ctime > AllClusters{j}{k}.savetime(end)
            fprintf('Cluster %d Details newer than Cluster!!\n',k);
            AllClusters{j}{k}.needed = 1;
        end
        for n = 1:length(f)
            if isfield(ClusterDetails{k},f{n})
                AllClusters{j}{k}.(f{n}) = ClusterDetails{k}.(f{n});
            end
        end
        if isfield(ClusterDetails{k},'next')
        for c = 1:length(ClusterDetails{k}.next)
            %only load ClusterDetails if Cluster is still defined in next{c}
            if c <= length(AllClusters{j}{k}.next) && isfield(AllClusters{j}{k}.next{c},'space')
            for n = 1:length(f)
                if isfield(ClusterDetails{k}.next{c},f{n})
                    AllClusters{j}{k}.next{c}.(f{n}) = ClusterDetails{k}.next{c}.(f{n});
                end
            end
            end
        end
        end
 %some old files saved as row, not column
        if isfield(AllClusters{j},'clst') && size(AllClusters{j}(k).clst,1) == 1
            AllClusters{j}(k).clst = AllClusters{j}(k).clst';
        end
        if rawxy && isfield(ClusterDetails{k},'xy')
            xy = ClusterDetails{k}.xy;
            C = AllClusters{j}{k};
            if C.shape == 0
                AllClusters{j}{k}.xy = xy;
            else
                AllClusters{j}{k}.xy = xyrotate(xy(:,1),xy(:,2),-C.angle);
            end
        end
        if alltimes && isfield(AllClusters{j}{k},'t')
            AllClusters{j}{k}.times = AllClusters{j}{k}.t;
            AllClusters{j}{k} = rmfield(AllClusters{j}{k},'t');
        end
    end
    for k = 1:length(AllClusters{j})
        if isfield(AllClusters{j}{k},'space')
            if ~isfield(AllClusters{j}{k},'quick')
                AllClusters{j}{k}.quick = 0;
            end
            if ~isfield(AllClusters{j}{k},'dropi')
                AllClusters{j}{k}.dropi = [0 0 0 0];
            end
            if ~isfield(AllClusters{j}{k},'trigdt')
                AllClusters{j}{k}.trigdt = 0;
            end
            if ~isfield(AllClusters{j}{k},'manual')
                AllClusters{j}{k}.manual = 0;
            end
            if ~isfield(AllClusters{j}{k},'next')
                AllClusters{j}{k}.next = {};
            end
            if ~isfield(AllClusters{j}{k},'clusterprog')
                AllClusters{j}{k}.clusterprog = '';
            end
        end
        if isfield(AllClusters{j}{k},'times') && diff(size(AllClusters{j}{k}.times)) < 0
            cprintf('red','%s P%d times is a row vector\n',name,k);
            AllClusters{j}{k}.times = AllClusters{j}{k}.times';
        end
        if isfield(AllClusters{j}{k},'t') && diff(size(AllClusters{j}{k}.t)) < 0
            cprintf('red','%s P%d t is a row vector\n',name,k);
            AllClusters{j}{k}.t = AllClusters{j}{k}.t';
        end
    end
    details{j}.loadur = mytoc(ts);
end

if length(AllClusters) == 1
    AllClusters = AllClusters{1};
    if length(AllFullVData)
    AllFullVData = AllFullVData{1};
    end
    details = details{1};
end

