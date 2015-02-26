function [name, dir] = splitpath(path, varargin)
%[name, dir] = splitpath(path, varargin) replicate fileparts
%
getcellname = 0;
j = 1;
while j < nargin
   if strncmpi(varargin{j},'cellpref',5)
       getcellname = 2;
   elseif strncmpi(varargin{j},'cell',4)
       getcellname = 1;
   elseif strncmpi(varargin{j},'dir',3)
       getcellname = 3;
   end
    j = j + 1;
end

if iscell(path)
    for j = 1:length(path)
        [dir{j}, name{j}, ext] = fileparts(path{j});
        name{j} = [name{j} ext];
    end
    return;
end

idx = findstr(path,'/');
if isempty(idx)
    idx = findstr(path,'\');
end
if isempty(idx)
    name = path;
    dir = '';
else
    name = path(idx(end)+1:end);
    dir = path(1:idx(end)-1);
end

if getcellname == 1
    idx = regexp(name,'\.[0-9]\.');
    if ~isempty(idx)
        name = name(1:idx+4);
    end
elseif getcellname == 2
    idx = regexp(name,'\.[0-9]\.');
    if ~isempty(idx)
        name = name(1:idx+1);
    end
elseif getcellname == 3
    idx = regexp(path,'[\/\\][0-9]*[\\\/]');
    last = regexp(path,'[\/\\][0-9]*[\\\/]','end');
    if ~isempty(idx)
        name = path(idx+1:last-1);
    end
end