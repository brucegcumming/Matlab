function ListTree(Tree, pattern)
%Go through AGB disk listings to find files
%These are typically on bgc/group/Disk\ Trees
%
listall = 1;
if length(pattern) > 3
    listall = 1;
end
if ischar(Tree)
    if isdir(Tree)
        d = mydir([Tree '/*.mat']);
        for j = 1:length(d)
            fprintf('Checking %s\n',d(j).name);
            ListTree(d(j).name,pattern);
        end
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
    if isfield(Tree,'parent') && ~isempty(regexp(Tree(j).parent,pattern))
        fprintf('%s/%s\n',Tree(j).parent,Tree(j).name);
    elseif ~isempty(regexp(Tree(j).name,pattern)) && (~isempty(Tree(j).children) || listall)
        fprintf('%s\n',Tree(j).name);
    end
    for k = 1:length(Tree(j).children)
        Tree(j).children(k).parent = Tree(j).name;
    end
    ListTree(Tree(j).children,pattern);
end