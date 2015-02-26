function names = ListTree(Tree, pattern)
%ListTree(Tree, pattern) Go through AGB disk Listings
%to find files
%These are typically on b/group/Disk Trees
%
listall = 1;
if length(pattern) > 3
    listall = 1;
end
if ischar(Tree)
    if isdir(Tree)
        if strcmp(pattern,'summary')
            names = ListTree(pattern,'lem');
            names = ListTree(pattern,'jbe');
        else
            d = mydir([Tree '/*.mat']);
            for j = 1:length(d)
                fprintf('Checking %s\n',d(j).name);
                X = ListTree(d(j).name,pattern);
                names(j) = X;
                names(j).dir = d(j).name;
            end
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

names.names = {};
names.nmatch = [];
nf = 0;
nmatch = 0;
nfullv = 0;
isdatdir = 0;
for j = 1:length(Tree)
    if isfield(Tree,'parent') && ~isempty(regexp(Tree(j).parent,pattern))
        nf = nf+1;
        fprintf('%s/%s\n',Tree(j).parent,Tree(j).name);
        names.names{nf} = sprintf('%s/%s',Tree(j).parent,Tree(j).name);
        fprintf('%s\n',names.names{nf});
        isdatdir = 1;
    elseif ~isempty(regexp(Tree(j).name,pattern)) && (~isempty(Tree(j).children) || listall)
        fprintf('%s\n',Tree(j).name);
        nmatch = nmatch+1;
        parent = Tree(j).parent;
    end
    if ~isempty(strfind(Tree(j).name,'FullV'))
        nfullv = nfullv+1;
    end
    for k = 1:length(Tree(j).children)
        Tree(j).children(k).parent = Tree(j).name;
    end
    newnames = ListTree(Tree(j).children,pattern);
    if ~isempty(newnames.names)
        names.names = {names.names{:} newnames.names{:}};
        names.nmatch = [names.nmatch newnames.nmatch];
    else
        names.nmatch = newnames.nmatch; %Number of matching children
        names.nfullv = newnames.nfullv; %Number of matching children
    end
    if isfield(Tree(j),'parent');
        names.dir = Tree(j).parent;
    end
end
if ~isfield(names,'dir')
    names.dir = [];
end
names.nfullv = nfullv;
if nmatch > 0
    names.nmatch = nmatch;
    fprintf('%d matches in %s\n',nmatch,parent);
end