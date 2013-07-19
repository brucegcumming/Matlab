function d = mydir(path, pattern,varargin)
% d = mydir(path)
%works just like dir, but the name field contains full path.
%In order to avoid trouble with wildcards, is path is just hte
%directory name, have '/' at the end
%if path is a cell array of strings, calls mydir for each
pattern = [];

if iscellstr(path)
    for j = 1:length(path)
        if ~isempty(pattern)
            d{j} = mydir([path{j} '/' pattern],varargin{:});
        else
            d{j} = mydir(path{j},varargin{:});
        end
        good(j) = ~isempty(d{j});
    end
    d = cat(1,d{good});
else
root = fileparts(path);
d = dir(path);
for j = 1:length(d)
    d(j).name = [root '/' d(j).name];
end
end
