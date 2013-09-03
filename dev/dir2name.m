function name = dir2name(path, type)
% name = dir2name(path, type)
% construct data/analysis filename from directory name

   [a,b] = fileparts(path);
   monkey = GetMonkeyName(path);
if strcmp(type,'rf')
   name = [path '/' monkey b '.rf.mat'];
else
   name = [path '/' monkey b '.mat'];
end