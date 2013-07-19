if ispc
    home='Z:/bondya/matlab/dev';
    cd(home)
    path('Z:/group/matlab',path);
    path('Z:/bgc/matlab/BlackRock',path);
    path('Z:/bgc/matlab',path);
    path('Z:/bgc/anal/orbw/',path);
    devpath
    path('Z:/bondya/matlab',path);
    path(home,path);
elseif ismac
    home='/sc/bgc/bondya/matlab/dev';
    if isdir(home)
        cd(home)
        path('/sc/bgc/group/matlab',path);
        path('/sc/bgc/bgc/matlab/BlackRock',path);
        path('/sc/bgc/bgc/matlab',path);
        path('/sc/bgc/bgc/anal/orbw',path);
        devpath;
        path('/sc/bgc/bondya/matlab',path);
        path(home,path);
    else
        home='/Volumes/bgc/bondya/matlab/dev';
        if isdir(home)
            path('/Volumes/bgc/group/matlab',path);
            path('/Volumes/bgc/bgc/matlab/BlackRock',path);
            path('/Volumes/bgc/bgc/matlab',path);
            path('/Volumes/bgc/bgc/anal/orbw',path);
            path('/Volumes/bgc/bgc/matlab/dev',path);
            path('/Volumes/bgc/bondya/matlab',path);
            path(home,path);
        end
    end
            
end
fprintf(['Hi Adrian.  Path is set.  Working directory is',' ',home,'\n']);