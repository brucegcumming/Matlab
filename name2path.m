function [path, dirpath] = name2path(name,varargin)

prefix = '/bgc/data/';
global bgcfileprefix
checkreal = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'smr',3)
        prefix = '/smr/';
    elseif strncmpi(varargin{j},'check',3)
        checkreal = 1;
    elseif strncmpi(varargin{j},'prefix',3)
        j = j+1;
        prefix = varargin{j};
    end
    j = j+1;
end
if strncmpi(name,'duf',3)
    monkey = 'dufus';
    dir = strrep(name,'duf','');
elseif strncmpi(name,'ruf',3)
    monkey = 'rufus';
    dir = strrep(name,'ruf','');
elseif strncmpi(name,'lem',3)
    monkey = 'lem';
    dir = strrep(name,'lem','');
end
dots = strfind(dir,'.');
if isempty(dots)
    dir = dir;
else
    dir = dir(1:dots(1)-1);
end
path = [prefix monkey '/' dir '/' name];
dirpath = [prefix monkey '/' dir];
if exist('bgcfileprefix','var') & length(bgcfileprefix)
    path = [bgcfileprefix path];
    dirpath = [bgcfileprefix dirpath];
end