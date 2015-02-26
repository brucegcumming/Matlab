function DATA = ListSubExpts(DATA, id, varargin)

setv = 1;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'relist',4)
setv = get(DATA.elst,'value');
end
j = j+1;
end
suffs = DATA.suffs;
na = 1;
subexplist = {};
expid = [];
exptnos = [];
if id == 0
id = get(DATA.clst,'value');
end
nrp = zeros(size(DATA.explist));
for j = 1:length(DATA.Expts);
if isfield(DATA.Expts{j},'Header')
exptno = GetExptNumber(DATA.Expts{j});
if exptno <1
exptno = j;
end
eid = strmatch(DATA.Expts{j}.Header.expname,{DATA.explist{id}},'exact');
ok = 1;
if  DATA.listbycell == 1 && length(eid)
if strncmp(DATA.filetype,'Grid',4)
rowid = DATA.cellexid(j);
else
cid = regexp(DATA.Expts{j}.Header.Name,'\.[0-9]*\.mat');
en = sscanf(DATA.Expts{j}.Header.Name(cid+1:end),'%d');
rowid = find(DATA.CellDetails.exptids == en);
end
if isempty(rowid) || rowid == 0
ok = 0;
else
pid = find(DATA.CellList(rowid(1),:) == DATA.probe);
if length(pid) ~= 1 %> 1 if cell defined twice
ok = 0;
end
end

else
pid = DATA.probe;
end
else
ok = 0;
eid = [];
end
if  ( id == 1 | ~isempty(eid)) & isfield(DATA.Expts{j},'Trials') & ok
tid = regexp(DATA.Expts{j}.Header.Name,'Expt[0-9]*.mat');
if ~isempty(tid) %online file
label = sprintf(' (%s:%d-%d)',DATA.Expts{j}.Header.Name(tid:end-4),...
DATA.Expts{j}.Trials(1).Trial,...
DATA.Expts{j}.Trials(end).Trial');
elseif DATA.show.times
label = sprintf(' %.1f-%.1f (id %d-%d)',DATA.Expts{j}.Trials(1).Start(1)./10000,...
DATA.Expts{j}.Trials(end).Start(1)./10000,DATA.Expts{j}.Trials(1).id,DATA.Expts{j}.Trials(end).id);
else
label = sprintf(' %d-%d (id %d-%d)',DATA.Expts{j}.Trials(1).Trial,...
DATA.Expts{j}.Trials(end).Trial,DATA.Expts{j}.Trials(1).id,DATA.Expts{j}.Trials(end).id);
end
if DATA.Expts{j}.Header.psych
label = [label ' P'];
end
if DATA.Expts{j}.Header.rc
label = [label ' RC'];
end
if DATA.listbycell
label = [label 'p' num2str(pid)];
end
label = [label cmb.ShowString(DATA,DATA.Expts{j})];
if id == 1
expi = strmatch(DATA.Expts{j}.Header.expname, DATA.explist,'exact');
nrp(expi) = nrp(expi)+1;
subexplist{na} = [DATA.Expts{j}.Header.expname num2str(exptno) label];
else
subexplist{na} = [DATA.explist{id(eid)} suffs(mod(na-1,length(suffs))+1) label];
subexplist{na} = [DATA.explist{id(eid)} num2str(exptno) label];
end
DATA.explabels{j} = subexplist{na};
if DATA.state.online && strncmp(DATA.filetype,'Grid',4)
if ismember(j,DATA.grididx.expt) % have .ns5 files
subexplist{na} = [subexplist{na} '*'];
end
if ismember(j,DATA.fullvidx) % have .ns5 files
subexplist{na} = [subexplist{na} '(V)'];
end
elseif DATA.Expts{j}.gui.clustertype == 0
subexplist{na} = [subexplist{na} '*'];
elseif DATA.Expts{j}.gui.clustertype == 2
subexplist{na} = [subexplist{na} '(O)'];
elseif DATA.Expts{j}.gui.ncluster == 0
subexplist{na} = [subexplist{na} '(Z)'];            
end
expid(na) = j;
exptnos(na) = exptno;
na = na+1;
elseif strmatch('unknown',{DATA.explist{id}})
subexplist{na} = [DATA.explist{id} suffs(na)];
expid(na) = j;
na = na+1;        
end
end
p  = get(DATA.elst,'Listboxtop');
if p > length(subexplist)
set(DATA.elst,'ListboxTop',1);
end
set(DATA.elst,'string',subexplist,'value',setv);
DATA.expid = expid;

DATA.subexplist = subexplist;
if ~isfield(DATA,'currentexpt')
DATA.currentexpt = 1;
end


