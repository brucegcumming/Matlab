function res = SetClusters(name, varargin)
%modifies ClusterTimesFiles on disk
%SetClusters(name, 'badexpt')  marks all Clusters in a file a Bad Expts\
%SetClusters(dir, 'goodprobe', p) Removes any bad probe marks for probe p. 
ops = {};
probes = [];
savefix = 0;
res = [];
j = 1; 
while j <= length(varargin)
    if strncmpi(varargin{j},'badexpt',5)
        ops{end+1} = 'badexpt';
    elseif strncmpi(varargin{j},'usearray',5)
        ops{end+1} = 'usearray';
    elseif strncmpi(varargin{j},'save',4)
        savefix = 1;
    end
    j = j+1;
end

if iscellstr(name)
    for j = 1:length(name)
        res{j} = SetClusters(name{j},varargin{:});
    end
    return;
end

res.name = name;
if strcmp(ops{1},'usearray')
    Array = GetArrayConfig(name);
    if ~isfield(Array,'X')
        return;
    end
    probes = 1:length(Array.X);
    res.probes = probes;
end

if isdir(name)
    d = mydir(sprintf('%s/Expt*ClusterTimes.mat',name));
    for j = 1:length(d)
        exn = GetExptNumber(d(j).name);
        load(d(j).name);
        for p = probes(:)';
            if strcmp(ops{1},'goodprobe')
                if isfield(Clusters{p},'marked') && Clusters{p}.marked == 3
                    Clusters{p}.marked = 0;
                end
            elseif strcmp(ops{1},'usearray')
                if isfield(Clusters{p},'exptno')
                    e = Clusters{p}.exptno;
                    if e ~= exn
                        fprintf('Cluster %d exptno mismatch\n',p);
                        fix(exn,p) = 3;
                    end
                else
                    if ~isfield(Clusters{p},'space')
                        fprintf('Cluster %d empty\n',p);
                        res.fix(exn,p) = 4;
                    else
                        fprintf('Cluster %d missing exptno\n',p);
                        res.fix(exn,p) = 2;
                    end
                end
                if isfield(Clusters{p},'marked') && ismember(Clusters{p}.marked,[3, -1000])
                    if (p > length(Array.badprobes) || Array.badprobes(p) == 0) && ...
                            (e > length(Array.badexpts) || Array.badexpts(e) == 0)
                        fprintf('%s Resetting E%dP%d to unmarked\n',d(j).name,e,p);
                        Clusters{p}.marked = 0;
                        res.fix(e,p) = 1;
                    end
                end
            end
        end
        res.expts(j) = e;
        if savefix
            save(d(j).name,'Clusters','FullVData');
        end
    end
elseif exist(name)
    load(name);
    for j = 1:length(ops)
        for c = 1:length(Clusters)
            if strcmp(ops{j},'badexpt')
                Clusters{c}.marked = -1000;
                Clusters{c}.badexpt = 1;
                Clusters{c}.savetime(3) = now;
            end
        end
    end
    save(name,'Clusters','FullVData');
end