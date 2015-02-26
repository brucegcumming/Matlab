function CompareLists(lista, listb, varargin)
%Compare contents of two lists of files

if ischar(lista) && exist(lista,'file')
    fid(1) = fopen(lista);
    a = textscan(fid(1),'%s','delimiter','\n');
    atxt = a{1};
    fclose(fid(1));
elseif iscellstr(lista)
    atxt = lista;
end

if ischar(listb) && exist(listb,'file')
    fid(2) = fopen(listb);
    a = textscan(fid(2),'%s','delimiter','\n');
    btxt = a{1};
    fclose(fid(2));
end

anames = lst2names(atxt);
bnames = lst2names(btxt);

for j = 1:length(anames)
    id = find(strcmp(anames{j},bnames));
    if isempty(id)
        missing(j) = 1;
    else
        missing(j) = 0;
    end
end

mid = find(missing);
fprintf('Missing %s\n',anames{mid});

for j = 1:length(bnames)
    id = find(strcmp(bnames{j},anames));
    if isempty(id)
        missing(j) = 1;
    else
        missing(j) = 0;
    end
end
mid = find(missing);
fprintf('A Missing %s\n',bnames{mid});

function names = lst2names(lst)
names = {};

for j = 1:length(lst)
    [a,b,c] = fileparts(strrep(lst{j},'\','/'));
    names{j} = b;
end