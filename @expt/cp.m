function cp(src, tgt, varargin)
%expt.cp(src, tgt,...) copy expt files
%expt.cp(src, tgt,spikes) copy expt files

dospikes = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'spikes')
        dospikes = 1;
    end
    j = j+1;
end

needfiles = {'ArrayConfig.mat' 'Errors.mat' 'CellList.mat'};

if ~exist(tgt)
    if confirm(sprintf('Create Folder %s?',tgt))
        mkdir(tgt);
    else
        return;
    end
end

d = mydir([src '/*.mat']);
for j = 1:length(d)
    go = 0;
    srcfile = d(j).name;
    if ~isempty(strfind(d(j).name,'ClusterTimes')) ||  ~isempty(strfind(srcfile,'idx.mat')) || sum(strcmp(d(j).filename,needfiles))
        tgtfile = strrep(d(j).name,src,tgt);
        go = 1;
    end
    if go
        fprintf('Copying %s to %s\n',d(j).name,tgtfile);
        copyfile(d(j).name,tgtfile);
    end
end
if dospikes
        copyfile([src '/Spikes' ,[tgt '/Spikes']);
end