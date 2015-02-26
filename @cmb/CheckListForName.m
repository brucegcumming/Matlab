function [id, ok] = CheckListForName(list,name)

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


