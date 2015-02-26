function [path, dirpath] = name2path(name,varargin)
%[path, dirpath] = name2path(name,...) Generate full path from name
%uses default prefix '/b/data/'
%   ....,'check')
%   ....,'prefix',P) Force prefix to P

prefix = '/b/data/';
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
id = strfind(name,'/');
if ~isempty(id)
    name = name(id(end)+1:end);
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
    if regexp(dir,'^[S,0-9]')
    dir = ['SE/' dir];
    end
elseif strncmpi(name,'jbe',3)
    monkey = 'jbe';
    dir = strrep(name,'jbe','');
    if regexp(dir,'^[S,0-9]')
    dir = ['SE/' dir];
    end
else
    dir = fileparts(name);
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