function DATA = SetSpkLists(DATA)
if isdir(DATA.datafilename)
return;
end
DATA.spklist = [];

for j = 1:length(DATA.Expts)
DATA.Expts{j}.gui.setispk = 0;
end

ifile = strrep(DATA.datafilename,'.mat',sprintf('.p%dispk.mat',DATA.probe));
if exist(ifile,'file')
load(ifile);
if length(ispklen) == length(DATA.Expts)
for j = 1:length(DATA.Expts)
if length(expispk{j})
DATA.Expts{j}.gui.spks = expispk{j};
DATA.Expts{j}.gui.setispk = length(expispk{j});
DATA.Expts{j}.gui.spkrange = [expispk{j}(1) expispk{j}(end)];
end
end
end
end
nset = 0;
for j = 1:length(DATA.Expts)
if DATA.Expts{j}.gui.setispk == 0
DATA = SetExptSpikes(DATA, j, 'setrange');
expispk{j} = DATA.Expts{j}.gui.spks;
ispklen(j) = length(expispk{j});
nset = nset+1;
end
end
if nset > 0
try
save(ifile,'expispk','ispklen');
end
end

