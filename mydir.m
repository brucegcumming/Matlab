function d = mydir(path, pattern,varargin)
% d = mydir(path)
%works just like dir, but the name field contains full path.
%In order to avoid trouble with wildcards, if path is just hte
%directory name, have '/' at the end
%if path is a cell array of strings, calls mydir for each

if nargin < 2
    pattern = [];
end
useregexp = 0;
if iscellstr(path)
    good = [];
    for j = 1:length(path)
        if ~isempty(pattern)
            d{j} = mydir([path{j} '/' pattern],varargin{:});
        else
            d{j} = mydir(path{j},varargin{:});
        end
        good(j) = ~isempty(d{j});
    end
    d = cat(1,d{find(good)});
else
    if isdir(path)
        root = path;
    else
        [root, name] = fileparts(path);
    end
    d = dir(path);
    if useregexp
    if isempty(d) %may have used regexp
        nx = 0;
        x  = dir(root);
        clear d;
        for j = 1:length(x)
            if regexp(x(j).name,name);
                nx = nx+1;
                x(j).filename = [root '/' x(j).name];
                d(nx) = x(j);
            end
        end
        return;
    end
    end
    if ~isempty(pattern)
        good = [];
        for j = 1:length(d)
            if regexp(d(j).name,pattern)
                good(j) = 1;
            end
        end
        d = d(good);
    end
    for j = 1:length(d)
        d(j).filename = d(j).name;
        if sum(strcmp(d(j).name,{'..' '.'}))
            good(j) = 0;
            else
            d(j).name = [root '/' d(j).name];
            good(j) = 1;
        end
        
    end
    if ~isempty(d)
    d = d(find(good));
    end
end
