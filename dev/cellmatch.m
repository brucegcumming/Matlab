function id = cellmatch(s, C, varargin)
%id = cellmatch(s, C) returns and index of elements in C that match s
%like strmatch, but works on all cell arrays, not just cell string arrays

id = [];
for j = 1:length(C)
    if strncmp(s, C(j),length(s))
        id = [id j];
    end
end
