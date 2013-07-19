function path = CheckPath(inpath, varargin)

%out = CheckPath(path) tries to ensure that paths for data files are corret
% adding /bgc to paths that lack it - can be missing on windows machines.

inpath = strrep(inpath,'\','/'); %works on all archs
if ~exist(inpath,'file') & ~strncmp(inpath,'/bgc',4)
    path = ['/bgc' inpath];
else
    path = inpath;
end
    