function [name, dir, fitname] = splitname(path)
% [dir, name, fitname] = Splitname(path)
% Splits pathname into directory and filename
% also returns ./fits/name for convenience
%

idx = findstr(path,'/');
name = path(idx(end)+1:end);
dir = path(1:idx(end));
fitname = sprintf('./fits/%s',name);
  