function ReWriteClusters(name, varargin)

if iscell(name)
    for j = 1:length(name)
        fprintf('%s\n',name{j});
        ReWriteClusters(name{j},varargin{:});
    end
    return;
end

savefiles = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4)
        savefiles = 1;
    end
    j = j+1;
end


exptno = -1;

load(name);
Clusters = Clusters;
dname = strrep(name,'.mat','Details.mat');
if ~exist(dname,'file')
dname = strrep(name,'ClusterTimes','ClusterTimesDetails');
end    
load(dname);

if isempty(strfind(name,'Auto'))
    autoname =  strrep(name,'ClusterTimes','AutoClusterTimesDetails');
    a = load(autoname);
end
for j = 1:max([length(ClusterDetails) length(a.ClusterDetails)])
    if j > length(ClusterDetails) && Clusters{j}.auto == 1
        ClusterDetails{j} = a.ClusterDetails{j};
    elseif isempty(ClusterDetails{j}) && Clusters{j}.auto == 1
        ClusterDetails{j} = a.ClusterDetails{j};
    end
    if isfield(ClusterDetails{j},'clst') && ...
        length(Clusters{j}.times) > Clusters{j}.ncut && ...
            length(Clusters{j}.times) == length(ClusterDetails{j}.clst)
        id = find(ClusterDetails{j}.clst > 1);
        Clusters{j}.times = Clusters{j}.times(id);
        fprintf('Probe %d times %d -> %d\n',j,length(ClusterDetails{j}.t),length(Clusters{j}.times));
    end
    if isfield(Clusters{j},'clst') & isfield(ClusterDetails{j},'clst')
        Clusters{j} = rmfield(Clusters{j},'clst');
        fprintf('Probe %d Removing clst\n',j);
    end
end
if savefiles
    fprintf('Writing %s\n',name);
    save(name,'Clusters');
end
