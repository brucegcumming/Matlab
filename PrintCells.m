function PrintCells(C, varargin)
%PrintCells(C)
%prints a cell string array (more compact than format compact
matchpat = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'match',5)
        j = j+1;
        matchpat = varargin{j};
    end
    j = j+1;
end

for  j = 1:length(C)
    if isempty(matchpat) | regexp(C{j},matchpat)
    fprintf('%s\n',C{j});
    end
end