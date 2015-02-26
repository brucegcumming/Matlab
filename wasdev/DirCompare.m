function DirCompare(lpath, npath)


[a,b,c,d] = TreeFind(lpath);
gb = (1024 * 1024 * 1024); 
totalbytes = sum(b);
for j = 1:length(a)
    lnames{j} = strrep(a{j},'lpath','');
end
[nname,nsize,c,d] = TreeFind(npath);
for j = 1:length(nname)
    nnames{j} = strrep(nname{j},'npath','');
end

for j = 1:length(a)
    if strfind(a{j},'.smr');
        types{j} = 'smr';
    elseif strfind(a{j},'FullV.mat');
        types{j} = 'FullV';
    elseif strfind(a{j},'FullV.mat');
        types{j} = 'FullV';
    else
        types{j} = 'other';
    end
    
end

alltypes = unique(types);
for j = 1:length(alltypes);
    id = strmatch(alltypes{j},types);
    sz = sum(b(id))./gb;
    fprintf('%s: %d files, %.3fGb\n',alltypes{j},length(id),sz);
end