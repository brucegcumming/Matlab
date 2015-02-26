function name = dir2name(path, type)
% name = dir2name(path, type) superceded by BuildPath
% construct data/analysis filename from directory name

if nargin ==1
    type = 'mat';
end

   [a,b] = fileparts(path);
   if isempty(b)
       [a,b] = fileparts(a);
   end
   monkey = GetMonkeyName(path);
if strcmp(type,'rf')
   name = [path '/' monkey b '.rf.mat'];
elseif strcmp(type,'penlog')
    rfname = dir2name(path,'rf');
    if exist(rfname,'file')
        load(rfname);
    else
        ufl = ReadUfl(path);
    end
elseif strcmp(type,'filename')
    if strncmp(monkey,b,length(monkey))
        name = b;
    else
        name = [monkey b];
    end
else
   name = [path '/' monkey b '.mat'];
end