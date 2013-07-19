function AddToList(list,name)

%[id, ok] = CheckListforName(list,name);
fid = fopen(list,'a');
if fid > 0 
    if iscell(name)
        for j = 1:length(name)
           fprintf(fid,'%s\n',name{j});
        end
    else
    fprinf(fid,'%s\n',name);
    end
    fclose(fid);
end


function [id, ok] = CheckListforName(list,name)

fid = fopen(list,'r');
if fid >0 
    ok = 1;
    names = textscan(fid,'%s');
    files = splitpath(names{1});
    id = strmatch(splitpath(name),files);
    fclose(fid);
else
    ok = 0;
    id = [];
end

