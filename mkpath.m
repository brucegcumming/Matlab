function mkpath(pathname)

% mkpath(path)
% makes sure all the directories in path
% exist

slashes = strfind(pathname,'/');

for j = 1:length(slashes)
    dir = pathname(1:slashes(j));
    if ~exist(dir);
        fprintf('Making %s\n',dir);
        mkdir(dir);
    end
end
