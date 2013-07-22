os = computer;
if isempty(strfind(path,'/bgc/matlab:')) %check bgc/matlab is in path first
   path('/bgc/bgc/matlab',path); 
end
if strmatch(os,{'MAC' 'MACI' 'GLNX64' 'MACI64'})
    path('/Volumes/bgc/bgc/matlab/dev',path);
elseif strmatch(os,'GLNXA64') %lsr-mc2
    path('/sc/bgc/bgc/matlab/dev',path);
elseif exist('Z:/bgc/matlab','dir')
    path('Z:/bgc/matlab/dev',path);
else
    path('C:/bgc/bgc/matlab/dev',path);
end