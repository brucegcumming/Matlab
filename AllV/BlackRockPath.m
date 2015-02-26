function BlackRockPath()
%Add Blackrock scripts to path
if isempty(strfind(path,'BlackRock'))
    addpath([GetFilePath('bgcmatlab') '/BlackRock']);
end