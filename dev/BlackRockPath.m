function BlackRockPath()

if isempty(strfind(path,'BlackRock'))
    path(path,'/bgc/bgc/matlab/BlackRock');
end