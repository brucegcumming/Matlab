function LoadErrors(name, varargin)
%Load error files from disk
if ~isdir(name)
    name = fileparts(name);
end


errfile = [name '/Errors.mat'];
cerrfile = [name '/ClusterErrors.mat'];
