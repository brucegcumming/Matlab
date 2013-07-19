
fprintf('Hello Bruce from %s\n',mfilename('fullpath'));

if exist('Z:/bgc/matlab','dir')
    path(['Z:/bgc/matlab;' path]);
    path(path,'Z:/matlab');
    cd Z:/bgc/bgc/anal;
else    
    path(['/bgc/bgc/matlab;' path]);
    path(path,'/bgc/matlab');
    global bgcfileprefix;
    bgcfileprefix = 'C:\bgc';
    cd \bgc\bgc\anal;
end
% Now run startup.m from this path
%
startup
