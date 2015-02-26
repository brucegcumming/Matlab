function [id, pos] = cellstrfind(str, pat)
%cellstrfind(str, pat) calls strfind on each cell in a cell array str, to
%find pat
%if str is a char, and pat is a cellarray, searches str for each pat

id = [];
pos = [];
if ischar(str) && iscellstr(pat)
for j = 1:length(pat)
    a = strfind(str,pat{j});
    if ~isempty(a)
        id(end+1) = j;
        pos{end+1} = a;
    end
end
else
for j = 1:length(str)
    a = strfind(str{j},pat);
    if ~isempty(a)
        id(end+1) = j;
        pos{end+1} = a;
    end
end
end