os = computer;
if exist('/b/bgc/matlab','dir')
    addpath('/b/bgc/matlab');
    addpath('/b/bgc/matlab/dev');
elseif strmatch(os,{'MAC' 'MACI' 'GLNX64' 'MACI64'})
    if isempty(strfind(path,'/bgc/matlab;')) %check bgc/matlab is in path first
        path('/bgc/bgc/matlab',path);
    end
    path('/Volumes/bgc/bgc/matlab/dev',path);
elseif strmatch(os,'GLNXA64') %lsr-mc2
    path('/sc/bgc/bgc/matlab/dev',path);
elseif exist('Z:/bgc/matlab','dir')
    if isempty(strfind(path,'Z:\bgc\matlab;')) %check bgc/matlab is in path first
        path('Z:\bgc\matlab',path);
    end
    path('Z:/bgc/matlab/dev',path);
else
    path('C:/bgc/bgc/matlab/dev',path);
end