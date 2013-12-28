root='/b/'; 
home='bondya/matlab/dev';
if ispc
    cd Y:/
end
if isdir(root)
    if isdir([root,home])
        cd([root,home]);
        path([root,'group/matlab'],path);
        path([root,'bgc/matlab/BlackRock'],path);
        path([root,'bgc/matlab'],path);
        path([root,'bgc/anal/orbw'],path);
        path([root,'bgc/matlab/dev'],path);
        path([root,'bondya/matlab/ralf'],path);        
        path([root,'bondya/matlab'],path);
        path([root,home],path);
        fprintf(['Hi Adrian.  Path is set.  Working directory is',' ',root,home,'\n']);
    else
        PrintMsg(0,'Could not set path. %s does not exist.',home);
    end      
else
    PrintMsg(0,'Could not set path. %s does not exist.',root);
end