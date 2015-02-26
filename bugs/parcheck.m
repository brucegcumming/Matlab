function parcheck(name)

%try to replicate parallelization error

prefix = '/b/data/lem/M306/';
filenames{1} = [prefix 'Expt1FullV.mat'];
filenames{2} = [prefix 'Expt2FullV.mat'];
res.exptlist(1) = 1;
res.exptlist(2) = 2;
ns = length(res.exptlist);

args = {'tchan' [1:24] 'quantifyall'};
args = {args{:} 'nowatch' 'nocheck' 'noninteractive' 'verbose'};
workers = [];
parfor  (j = 1:ns)
    t = mygetCurrentTask();
    fprintf('Expt%d worker id is %d\n',res.exptlist(j),t.ID);
    cls{j} =  ProcessSuffix(filenames{j}, res.exptlist(j), args{:});
    workers(j) = t.ID;
end
res.cls = cls;
res.workerid = workers;

function res = ProcessSuffix(path, ex, varargin)
template = [];
nocut = 0;
res = {};
j = 1;
checkexpts = 0;
makesmall = 1;
quantify = 0;
checkspikes = 0;
userefcluster = 0;
probes = [];

args = {};

nid = [];
res.filename = path;
outname = path;
if exist(outname,'file')
    ts = now;
    fprintf('%s Exists\n',outname);
    d = dir(outname);
else
    fprintf('Cant Read %s\n',outname);
    return;
end

    cfile = strrep(outname,'FullV.mat','ClusterTimes.mat');
    if ~exist(cfile)
            fprintf('No Cluster File %s\n',cfile);
            return;
    else
        load(cfile);
        ClusterDetails = {};
        for j = 1:length(Clusters)
            [needc(j), reason{j}] = AllV.NeedToQuantify(Clusters{j}, outname, []);
        end
        if sum(needc) == 0
            fprintf('No Clusters need quantifying in %s\n',cfile);
            return;
        else
            nid = find(needc > 0);
            fprintf('Need %s (%s) in %s\n',sprintf('%d ',nid),sprintf('%d ',needc(nid)),cfile);
        end
    end

