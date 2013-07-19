function SaveArrayConfig(DATA)

ArrayConfig.X = ones(1,DATA.nprobes);
ArrayConfig.Y = 1:DATA.nprobes;
ArrayConfig.id = 1:DATA.nprobes;

if isfield(DATA,'datadir')
    datdir = DATA.datadir;
elseif isdir(DATA.name)
    datdir = DATA.name;
end
name = [datdir '/ArrayConfig.mat'];
if ~exist(name)
    fprintf('Saving %s\n',name);
    save(name,'ArrayConfig');
end