function DATA = LoadClusters(DATA, cfile, varargin)
clearold = 1;
allprobes = 0;
nprobes = 1;
probe = DATA.probe;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{1},'allprobes',4)
        allprobes = 1;
    elseif strncmpi(varargin{j},'noclear',5)
        clearold = 0;
    elseif strncmpi(varargin{j},'probe',5)
        j =j+1;
        probe = varargin{j};
    end
    j = j+1;
end
if ~isfield(DATA.AllData,'Spikes')
    return;
end
if isempty(DATA.AllData.Spikes) && isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{probe};
else
    Spks = DATA.AllData.Spikes;
end
if isempty(Spks)
    return;
end

if exist(cfile,'file')
    excludelist = [];
    clid = [];
    classified = 0;
    load(cfile);
    DATA.state.recut = 1;
    firstspk = length(Spks.times);
    lastspk = 0;
    if size(clid,1) == size(Spks.codes,1);
        Spks.codes(:,2) = clid;
        last = max(find(clid > 0));
        maxclasst = Spks.times(last);
        classified = 1;
    end
    for j = 1:min([length(excludelist) length(DATA.Expts)]);
        if max(excludelist{j}) <= length(DATA.Expts{j}.Trials)
            lst = -abs([DATA.Expts{j}.Trials(excludelist{j}).Trial]);
            if lst
                for k = 1:length(lst)
                    [DATA.Expts{j}.Trials(excludelist{j}(k)).Trial] = lst(k);
                end
            end
            
        end
    end
    p = probe;
    %should chnage this to use AllClusters, then find correct expt
    %especially for loading clusters define online, where Expt order can
    %change.
    nx=0;
    for j = 1:length(DATA.Expts)
        eid(j,1) = DATA.Expts{j}.Trials(1).id;
        eid(j,2) = DATA.Expts{j}.Trials(end).id;
    end
    nset = 0;
    for j = 1:length(AllClusters)
        %only load the data for the current probe if this is a single probe file
        if ~isempty(AllClusters{j})
            if isfield(AllClusters{j},'ids')
                exid = find(eid(:,1) < AllClusters{j}.ids(2) & ...
                    eid(:,2) > AllClusters{j}.ids(1));
                if length(exid) == 1
                    nx = exid;
                elseif length(exid) > 1
                    overlap = [];
                    for k = 1:length(exid)
                        overlap(k) = sum(ismember([eid(exid(k),1):eid(exid(k),2)],[AllClusters{j}.ids(1):AllClusters{j}.ids(2)]));
                    end
                    [a,k] = max(overlap);
                    nx = exid(k);
                end
            elseif j <= length(DATA.Expts)
                nx = j;
            end
            if nx > 0 && isfield(AllClusters{j},'Cluster')
                sz = size(AllClusters{j}.Cluster);
                if allprobes
                    DATA.Expts{nx}.Cluster = AllClusters{j}.Cluster;
                    % some old single probe files have AllClusters the other way around
                elseif sz(1) == 1 && sz(2) > 1 && length(DATA.probes) == 1
                    for k = 1:sz(2)
                        DATA.Expts{nx}.Cluster{k,1} = AllClusters{j}.Cluster{k};
                    end
                elseif sz(1) >= 0 && sz(2) >= 1 && length(DATA.probes) == 1
                    DATA.Expts{nx}.Cluster(1:sz(1),1) = AllClusters{j}.Cluster(:,sz(2));
                elseif sz(1) > 0 && sz(2) >= p
                    DATA.Expts{nx}.Cluster(1:sz(1),p) = AllClusters{j}.Cluster(:,p);
                    for k = (sz(1)+1):size(DATA.Expts{nx}.Cluster,1)
                        DATA.Expts{nx}.Cluster{k,p} = [];
                    end
                    %If Allclsuters{j).Cluster is empty, make sure DATA.Expts is too.
                elseif sz(1) == 0 && sz(2) >= p  & isfield(DATA.Expts{nx},'Cluster')
                    for k = 1:size(DATA.Expts{nx}.Cluster,1)
                        DATA.Expts{nx}.Cluster{k,p} = [];
                    end
                end
                if exist('clustertypes','var') %%saved type
                    DATA.Expts{nx}.gui.clustertype = clustertypes(j);
                else
                    DATA.Expts{nx}.gui.clustertype = 1;
                end
                firstspk = length(Spks.times);
                lastspk = 0;
                nset = 0;
                if ~isfield(DATA.cluster{1,1},'Arange')
                    DATA.cluster{1,1}.Arange = [6:8];
                    DATA.cluster{1,1}.Brange = [11:20];
                    DATA.cluster{1,1}.Erange = [1:100];
                end
            end
            
            %only do this for one probe, as AllData.Spikes only has one probes data
            p = probe;
            if nx > 0 && isfield(DATA.Expts{nx},'Cluster')
                nprobes = size(DATA.Expts{nx}.Cluster,2);
                
                for k = 1:size(DATA.Expts{nx}.Cluster,1)
                    if p > nprobes
                        DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                        DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                        nprobes = p;
                    end
                    if isempty(DATA.Expts{nx}.Cluster{k,p}) | ~isfield(DATA.Expts{nx}.Cluster{k,p},'firstspk')
                        DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                        DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                    else
                        nset = nset + 1;
                    end
                    if DATA.Expts{nx}.Cluster{k,p}.firstspk == 0
                        DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                    end
                    if DATA.Expts{nx}.Cluster{k,p}.lastspk == 0
                        DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                    end
                    if DATA.Expts{nx}.Cluster{k,p}.firstspk < firstspk
                        ts = Spks.times(DATA.Expts{nx}.Cluster{k,p}.firstspk);
                        if ts  > DATA.Expts{nx}.Header.trange(1) && ts < DATA.Expts{nx}.Header.trange(2)
                            firstspk = DATA.Expts{nx}.Cluster{k,p}.firstspk;
                        end
                    end
                    if DATA.Expts{nx}.Cluster{k,p}.lastspk > lastspk && ...
                            DATA.Expts{nx}.Cluster{k,p}.lastspk <= length(Spks.times)
                        ls = Spks.times(DATA.Expts{nx}.Cluster{k,p}.lastspk);
                        if ls  > DATA.Expts{nx}.Header.trange(1) && ls < DATA.Expts{nx}.Header.trange(2)
                            lastspk = DATA.Expts{nx}.Cluster{k,p}.lastspk;
                        end
                    end
                    if DATA.Expts{nx}.gui.clustertype == 2 % online cluster - can't use spk coun
                        DATA.Expts{nx}.Cluster{k,p}.firstspk = NaN;
                        DATA.Expts{nx}.Cluster{k,p}.lastspk = NaN;
                    end
                    if ~isfield(DATA.Expts{nx}.Cluster{k,p},'params')
                        DATA.Expts{nx}.Cluster{k,p}.params = [1 2];
                    end
                    if ~isfield(DATA.Expts{nx}.Cluster{k,p},'Arange') & isfield(DATA.cluster{1,1},'Arange')
                        DATA.Expts{nx}.Cluster{k,p}.Arange = DATA.cluster{1,1}.Arange;
                        DATA.Expts{nx}.Cluster{k,p}.Brange = DATA.cluster{1,1}.Brange;
                        DATA.Expts{nx}.Cluster{k,p}.Erange = DATA.cluster{1,1}.Erange;
                    end
                    
                end
            end
            
            if lastspk < firstspk %nevet found these in cluster file
                firstspk = NaN;
                lastspk = NaN;
            end
            if lastspk > size(Spks.codes,1)
                lastspk = size(Spks.codes,1);
            end
            % with multiple probes, empty clsuters can arise because clusters have beed set
            % on other probes.
            if nx > 0
                if nset == 0 && nprobes == 1 % actively set no clusters. Only works if one probe
                    DATA.Expts{nx}.gui.classified = 1;
                    DATA.Expts{nx}.gui.ncluster = 0;
                elseif isnan(firstspk)
                    DATA.Expts{nx}.gui.classified = 0;
                    DATA.Expts{nx}.gui.ncluster = 0;
                elseif classified && lastspk > firstspk && max(Spks.codes(firstspk:lastspk,2)) > 0
                    DATA.Expts{nx}.gui.classified = classified;
                    DATA.Expts{nx}.gui.ncluster = nset;
                else
                    DATA.Expts{nx}.gui.classified = 0;
                    DATA.Expts{nx}.gui.ncluster = 0;
                end
            end
        end
    end
elseif clearold
    p = probe;
    for j = 1:length(DATA.Expts);
        if isfield(DATA.Expts{j},'Cluster')
            for k = 1:size(DATA.Expts{j}.Cluster,1)
                DATA.Expts{j}.Cluster{k,p} = {};
            end
        end
    end
end

if isempty(DATA.AllData.Spikes)
    DATA.AllSpikes{probe}.codes = Spks.codes;
else
    DATA.AllData.Spikes.codes = Spks.codes;
end
