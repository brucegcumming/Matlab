if ispc
    home='Z:/bondya/matlab/dev';
    cd(home)
    path('Z:/bgc/matlab/BlackRock',path);
    path('Z:/bgc/matlab',path);
    path('Z:/bgc/anal/orbw/',path);
    devpath
    path('Z:/bondya/matlab',path);
    path(home,path);
elseif ismac
    home='/sc/bgc/bondya/matlab/dev';
    cd(home)
    path('/sc/bgc/bgc/matlab/BlackRock',path);
    path('/sc/bgc/bgc/matlab',path);
    path('/sc/bgc/bgc/anal/orbw',path);
    devpath;
    path('/sc/bgc/bondya/matlab',path);
    path(home,path);
end
fprintf(['Hi Adrian.  Path is set.  Working directory is',' ',home,'\n']);