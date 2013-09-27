function ListTree(Tree, pattern)

listall = 1;
if length(pattern) > 3
    listall = 1;
end
if ischar(Tree)
    if isdir(Tree)
        return;
    elseif exist(Tree,'file')
    load(Tree);
    else
        return;
    end
elseif iscellstr(Tree)
    for j = 1:length(Tree)
        fprintf('%s\n',Tree{j});
        ListTree(Tree{j},pattern);
    end
end

for j = 1:length(Tree)
    if ~isempty(regexp(Tree(j).name,pattern)) && (~isempty(Tree(j).children) || listall)
        fprintf('%s\n',Tree(j).name);
    end
    ListTree(Tree(j).children,pattern);
end