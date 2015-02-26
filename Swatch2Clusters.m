function Swatch2Clusters(name, varargin)
%Swatch2Clusters(name, varargin) Convert old files to PlotClusters
%Not much use unless classification has been run on all probes/expts,
%because leaves holes
%Also not easy for PlotClsuters to read spikes because in Spike2 files
% the spikes file does not line up with expts - arbitraty division into 10

d  = dir([name '/Spikes/*.mat']);
cnames = {};
[a,b,c] = GetMonkeyName(name);
prefix = [b c];
probefile = [name '/' prefix 'probes.mat'];
load(probefile);
probelist = unique([probes.probe]);
matfile = [name '/' prefix '.mat'];
[~, Expts] = APlaySpkFile(matfile,'nospikes');

for j = 1:length(probelist)
    p = probelist(j);
    fprintf('Loading Probe %d\n',p);
    allt = [];
    values = [];
    pid = find([probes.probe] == p);
    cnames{p} = sprintf('%s/%s.p%dcl.mat',name, prefix, p);
    clst = [];
    if exist(cnames{p})
        X{p} = load(cnames{p});
        clst = X{p}.clid;
    else
        X{p} = [];
    end
    for t = 1:length(pid)
        spkname = [name '/Spikes/' probes(pid(t)).filename];
        S = load(spkname);
        Spk = S.(probes(pid(t)).var);
        ranges(t,:) = [1+length(allt) length(allt) + length(Spk.times)];
        allt = [allt; Spk.times];
        values = [values; Spk.values];
    end
    if length(allt) == length(clst)
        ncl = unique(clst);
        for t = 1:length(Expts)
            trange = Expts{t}.Header.trange./10000;
            aid = find(allt >= trange(1) & allt <= trange(2));
            id = find(clst ==1 & allt >= trange(1) & allt <= trange(2));
            nid = find(clst ==0 & allt >= trange(1) & allt <= trange(2));
            C{t,p}.times = allt(id);
            C{t,p}.MeanSpike.ms = mean(values(id,:));
            C{t,p}.MeanSpike.mu = mean(values(nid,:));
            C{t,p}.nspks = length(aid);
            C{t,p}.ncut = length(id);
%            C{t,p}.clst = clst(xid);
            for cl = 2:max(ncl)
                id = find(clst ==cl);
                if ~isempty(id)
                    C{t,p}.next{cl-1}.times = allt(id);
                end
            end
            C{t,p}.mahal = [0 0 0 0];
            C{t,p}.fitdprime = [0 1 1 0 0];
        end
    end
end
for e = 1:size(C,1)
    Clusters = C(e,:);
    for j = 1:length(Clusters)
        if ~isfield(Clusters{j},'nspks')
            Clusters{j}.nspks = 0;
        end
    end
    outname = sprintf('%s/Expt%dClusterTimes.mat',name,e);
    fprintf('Saving %s\n',outname);
    save(outname,'Clusters');
end

%DATA.AllData.Spikes = GetProbeFiles(DATA,DATA.probe,DATA.subprobe,'trange',DATA.Expts{eid}.Header.trange);


function cfile = ClusterFile(DATA,varargin)
getonline = 0;
getauto = 0;
allprobes = 0;
spkcodefile = 0;
probe = DATA.probe;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'auto',4) %Read Clusters from Online mat file
        getauto = 1;
    elseif strncmpi(varargin{j},'allprobes',4) %Read Clusters from Online mat file
        allprobes = 1;
    elseif strncmpi(varargin{j},'getonline',4) %Read Clusters from Online mat file
        getonline = 1;
    elseif strncmpi(varargin{j},'codes',4) %Save codes/times
        spkcodefile = 1;
    elseif strncmpi(varargin{j},'probe',4) %Read Clusters from Online mat file
        j = j+1;
        probe = varargin{j};
    end
    j = j+1;
end
if DATA.state.online == 0 %%not online data
    if getonline
        [a,b] = splitpath(DATA.datafilename);
        if length(b)
        cfile = [b '/OnlineClusters.mat'];
        else
        cfile = 'OnlineClusters.mat';
        end
    elseif getauto
        cfile = strrep(DATA.datafilename,'.mat','.autocl.mat');
    elseif length(DATA.probelist) > 1 && allprobes
        cfile = strrep(DATA.datafilename,'.mat','.allcl.mat');
    elseif spkcodefile 
        cfile = strrep(DATA.datafilename,'.mat','codes.mat');
    elseif length(DATA.probelist) > 1 
        cfile = strrep(DATA.datafilename,'.mat',['.p' num2str(probe) 'cl.mat']);
    else
        cfile = strrep(DATA.datafilename,'.mat','.cl.mat');
    end
else
    if spkcodefile 
        cfile = sprintf('%s/SpikeCodes%d.mat',DATA.datafilename,DATA.Expts{DATA.currentexpt(1)}.Trials(end).id);
    else
        cfile = [DATA.datafilename '/Clusters.mat'];
    end
end
