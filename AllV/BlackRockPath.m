function BlackRockPath()

if isempty(strfind(path,'BlackRock'))
    if exist('./BlackRock','dir')
        path(path,'./BlackRock');
    else
        path(path,'/bgc/bgc/matlab/BlackRock');
    end
end