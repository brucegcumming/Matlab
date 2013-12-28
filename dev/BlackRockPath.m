function BlackRockPath()

if isempty(strfind(path,'BlackRock'))
    addpath([GetFilePath('matlab') '/BlackRock']);
end