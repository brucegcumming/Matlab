ts = now;
os = computer;
if strmatch(os,{'MAC' 'MACI' 'MACI64'})
    path(path,'/b/group/matlab/agb');
    path(path,'/b/group/matlab');
    path('/b/bgc/matlab',path);
%    path('/b/bgc/matlab/dev',path);
    cd /b/bgc/matlab
elseif strmatch(os,'GLNXA64')
    path(path,'/b/group/matlab');
    path('/b/bgc/matlab',path);
    cd /b/bgc/anal;
elseif strmatch(os,{'PCWIN' 'PCWIN64'})
    if exist('Z:/bgc/matlab','dir')
        path('Z:/bgc/matlab',path);
        path(path,'Z:/group/matlab');
        cd Z:/bgc/matlab;
    else
        path('/b/bgc/matlab',path);
        path(path,'/b/group/matlab');
        cd /b/bgc/matlab;
    end
end


fprintf('Starting with %s under %s took %.2f\n',mfilename('fullpath'),os,mytoc(ts));
dbstop if error;

