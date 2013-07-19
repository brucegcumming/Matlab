function Clusters = LoadClusterFiles(name, varargin)
%LoadClusterFiles(dir, ......)
%reads in Cluster files made by AllVPcs and plots up clusters for each

verbose = 1;
cmpdrive = [];
forcewrite = 0;
TAGTOP = 'PlotClusters';
args = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        TOPTAG = varargin{j};
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
initcall =0 ;

if isstruct(name) || iscell(name)  %plot results
elseif isdir(name)
    d = dir(name);
    strings = {};
    xid = [];
    for j  = 1:length(d)
        if strfind(d(j).name,'ClusterTimes.mat') 
            strings = {strings{:} [name '/' d(j).name]};
            if strfind(d(j).name,'OldClusterTimes.mat')
                xid = [xid j];
            elseif strfind(d(j).name,'Copy of')
                xid = [xid j];
            end
        end
    end
    id = setdiff(1:length(strings),xid);
    strings = strings(id);
end

    for j = 1:length(strings)
        eid(j) = GetExptNumber(strings{j});
    end
    expts = unique(eid);
    for j = 1:length(strings)
        dname = strrep(strings{j},'.mat','Details.mat');
        if ~exist(dname,'file')
            dname = strrep(name,'ClusterTimes','ClusterTimesDetails.mat');
        end
        if verbose
        fprintf('Loading %s\n',strings{j});
        end
        X =  LoadClusters(strings{j});
        e = GetExptNumber(strings{j});
        e = find(expts == e);
        if strfind(strings{j},'Auto')
            AutoClusters{e} = X;
            for k = 1:length(AutoClusters{e})
                AutoClusters{e}{k}.auto = 1;
            end
        else
            Clusters{e} = X;
        end
    end
    for j = 1:length(expts)
        if j > length(Clusters) || isempty(Clusters{j})
            Clusters{j} = AutoClusters{j};
        else
            for k = 1:length(AutoClusters{j})
                if AutoClusters{j}{k}.savetime(1) >= Clusters{j}{k}.savetime(1) && Clusters{j}{k}.auto
                    fprintf('AutoCluster Newer for E%.0fP%.0f\n',expts(j),k);
                    Clusters{j}{k} = AutoClusters{j}{k};
                elseif ~isfield(Clusters{j}{k},'clst')
                    if Clusters{j}{k}.auto
                        Clusters{j}{k} = AutoClusters{j}{k};
                    end
                    fprintf('Missing clst for E%.0fP%.0f\n',expts(j),k);
                end
            end
        end
    end
    fprintf('Done\n');
    Clusters{1}{1}.prefix = name;
 
    
    function [Clusters, details] = LoadClusters(name)

exptno = -1;

load(name);
Clusters = Clusters;
for j = 1:length(Clusters)
    if isfield(Clusters{j},'exptno')
        exptno = Clusters{j}.exptno;
    end
end
if exptno < 0
    id = strfind(name,'Expt');
    if length(id)
        [exptno, a,b,c] = sscanf(name(id(1)+4:end),'%d');
        if name(id(1)+3+c) == 'a'
            exptno = exptno+0.1;
        end
    end
end
details.exptno = exptno;
dname = strrep(name,'.mat','Details.mat');
if ~exist(dname,'file')
dname = strrep(name,'ClusterTimes','ClusterTimesDetails.mat');
end    
load(dname);
for j = 1:length(ClusterDetails)
    if isfield(ClusterDetails{j},'xy')
        Clusters{j}.xy = ClusterDetails{j}.xy;
        if isfield(ClusterDetails{j},'clst')
            Clusters{j}.clst = ClusterDetails{j}.clst;
        elseif ~isfield(Clusters,'clst')
            Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
            id = find(ismember(ClusterDetails{j}.t,Clusters{j}.times));
            Clusters{j}.clst = ones(size(ClusterDetails{j}.t));
            Clusters{j}.clst(id) = 2;
        end
        Clusters{j}.times = ClusterDetails{j}.t;
        if exptno >= 0
            Clusters{j}.exptno = exptno;
        end
    end
    if isfield(ClusterDetails{j},'Evec')
        Clusters{j}.Evec = ClusterDetails{j}.Evec;
    end
end
