function DATA = ReadGridFile(DATA)

idx = [];
monks = {'lem' 'jbe'};
datdir = 'F:/Utah/dufus/';
datdir = fileparts(DATA.name);
plotsummary = 2;
[a,b,c,d] = fileparts(DATA.name);
BlackRockPath();
if DATA.state.online
id = strfind(DATA.name,'/');
if ~isempty(id)
dname = DATA.name(id(end)+1:end);
monk = dname(1:3);
dname = dname(4:end);
else
dname = test;
end
if strfind(a,'Spike')
datdir = ['F:/Utah/' monk '/' dname '/'];
else
datdir = DATA.name;
end
end
idx = BuildGridIndex(datdir, DATA.Expts, 'reindex');
DATA.grididx = idx;
id = find(idx.expt == DATA.currentexpt(1));
bid = id;
for j = 1:length(id)
nfiles{j} = [datdir '/' idx.names{id(j)}];
end
if DATA.state.usensx == 1 %using ns5 for spikes
elseif length(bid) > 1
DATA = ReadGridFiles(DATA,nfiles, 'toff',idx.toff(bid));
elseif length(bid) == 1
DATA = ReadGridFiles(DATA,nfiles, 'toff',idx.toff(bid));
%       DATA = ReadGridSpikes(DATA,nfiles{j},'toff',idx.toff(bid));
end
fprintf('%d files %d trials\n',length(id),sum(idx.nt(bid)));

