function [missing, filetypes] = CheckList(list,root,suff)
%
% CheckList runs checks on a list of .mat files for a project
% CheckList(list, path, suff) searches in path for all files with suffix
% suff, and sees if they are in the list. returns a list of files missing
% from the list.
filetypes = [];


if iscell(list)
    [files, dirs] = splitpath(list);
else
fid = fopen(list,'r');
if fid >0 
    ok = 1;
    names = textscan(fid,'%s');
    [files, dirs] = splitpath(names{1});
    fclose(fid);
else
    return;
end
end

if nargin == 3
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
return;
elseif strncmpi(root,'sxcx',4)
    for j = 1:length(dirs)
        sxs = dir([dirs{j} '/*IP*']);
        sfs = dir([dirs{j} '/*SF*mat']);
        if length(sxs)
            filetypes(j) = 1;
            missing{j} = sxs(1).name;
        elseif length(sfs)
            msfs = dir([dirs{j} '/*SF*mat']);
            filetypes(j) = 2;
            if length(msfs)
                missing{j} = [dirs{j} '/' msfs(1).name];
            else
                missing{j} = [dirs{j} '/' sfs(1).name];
            end
                
        else
            filetypes(j) = 0;
        end
    end
else
end