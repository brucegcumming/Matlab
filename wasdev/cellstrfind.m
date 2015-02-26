function [id, pos] = cellstrfind(str, pat)
%cellstrfind(str, pat) calls strfind on each cell in a cell array
%

id = [];
pos = [];
for j = 1:length(str)
    a = strfind(str{j},pat);
    if ~isempty(a)
        id(end+1) = j;
        pos{end+1} = a;
    end
end