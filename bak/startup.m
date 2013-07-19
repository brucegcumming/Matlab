
os = computer;
if strmatch(os,{'MAC' 'MACI' 'MACI64'})
    path(path,'/bgc/group/matlab');
    path('/Volumes/bgc/bgc/matlab',path);
    cd /bgc/bgc/anal;
elseif strmatch(os,'GLNXA64')
    path(path,'/sc/bgc/group/matlab');
    path('/sc/bgc/bgc/matlab',path);
    cd /sc/bgc/bgc/anal;
elseif strmatch(os,{'PCWIN' 'PCWIN64'})
    if exist('Z:/bgc/matlab','dir')
        path('Z:/bgc/matlab',path);
        path(path,'Z:/group/matlab');
        cd Z:/bgc/matlab/dev;
    else
        path('/bgc/bgc/matlab',path);
        path(path,'/bgc/group/matlab');
        cd /bgc/bgc/matlab/dev;
    end
end


fprintf('Starting with %s under %s\n',mfilename('fullpath'),os);

