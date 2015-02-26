function parforbug(name)
%The old PC file naming bug fixed in 2014a seems to have
%an unfixed version in parfor loops.  If prefix is
%/b/data/lem/M306 , without the drive name, then dir returns
%an empty list, even though load() works and exist() says true
%if prefix is set to
%Y:/b/data/lem/M306
%then it works.  
% if parallel is set to 0, both forms work. 

prefix = '/b/data/lem/M306/';
prefix = '/bgc/matlab/bugs/';
filenames{1} = [prefix 'Expt1ClusterTimes.mat'];
filenames{2} = [prefix 'Expt2ClusterTimes.mat'];
ns = length(filenames);

parallel = 1;
if parallel
    parfor  (j = 1:ns)
        cls{j} =  ProcessSuffix(filenames{j});
    end
else
    for  (j = 1:ns)
        cls{j} =  ProcessSuffix(filenames{j});
    end
end

function res = ProcessSuffix(cfile)

res.filename = cfile;
d = dir(cfile);
fprintf('%s dir finds %d matches\n',cfile,length(d));

if ~exist(cfile)
    fprintf('No Cluster File %s\n',cfile);
    return;
else
    load(cfile);
    fprintf('Loaded %s. %d Clusters\n',cfile,length(Clusters));
end

