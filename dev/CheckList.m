function missing = CheckList(list,root,suff)

fid = fopen(list,'r');
if fid >0 
    ok = 1;
    names = textscan(fid,'%s');
    files = splitpath(names{1});
    fclose(fid);
else
    return;
end

filenames = TreeFind(root,'name',suff);
k = 1;
missing = {};
for j = 1:length(filenames)
    id = strmatch(splitpath(filenames{j}),files);
    if isempty(id)
        missing{k} = filenames{j};
        k = k+1;
    end
end    