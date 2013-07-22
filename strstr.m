function [found, ids] = strstr(str, pattern)
%[found, ids] = strstr(str, pattern)
% if str is a string, calls strfind.
% if str is a cell array of strings, calls strfind for each element.
if iscellstr(pattern)
    for j = 1:length(pattern)
        id = strfind(str,pattern{j});
        if ~isempty(id)
            ids(j) = id(1);
        else
            ids(j) = 0;
        end
    end
    found = sum(ids);
elseif ischar(pattern)
    ids = strfind(str,pattern);
    if ~isempty(ids)
        found = ids(1);
    else
        found = 0;
    end
end