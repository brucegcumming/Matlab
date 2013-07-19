function [fstrings, details] = ReadCellList(list)
%read a list of .mat files for a list GUI
%will follow a list of lists
%removes lines begninnin with #
%removes commends following a space.

details = [];
fstrings = textread(list,'%s');
if strncmp(fstrings{1},'lists',5)
    allstrings = {};
    for j = 2:length(fstrings)
        newstrings = textread(fstrings{j},'%s');
        allstrings = {allstrings{:} newstrings{:}};
    end
    fstrings = allstrings;
end
laststr = 0;

for j = 1:length(fstrings)
    if fstrings{j}(1) == '#';
        deletestr(j) = 1;
        laststr = 1;
    elseif laststr & isempty(strfind(fstrings{j},'.mat')) %% still in comment
        deletestr(j) = 1;
    else
        goodstr(j) = 1;
        laststr = 0;
    end
end
fstrings = {fstrings{find(goodstr)}};

