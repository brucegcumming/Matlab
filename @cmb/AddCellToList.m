function AddCellToList(a,b, varargin)


DATA = GetDataFromFig(a); 
eid = DATA.exabsid;

itype =  1;
np = length(DATA.probelist);

for j = 1:DATA.state.listlen
cluster(j) = 1;
it = findobj(get(a,'parent'),'Tag',sprintf('CellProbeAdd%d',j));
if ~isempty(it)
probe(j) = get(it(1),'value');
if probe(j) > np
probe(j) = probe(j)-np;
cluster(j) = 2;
end
end
it = findobj(get(a,'parent'),'Tag',sprintf('CellNumber%d',j));
if ~isempty(it)
cell(j) = get(it(1),'value')-1;
end
it = findobj(get(a,'parent'),'Tag',sprintf('CellQuality%d',j));
if ~isempty(it)
quality(j) = get(it(1),'value');
end
end

it = findobj(get(a,'parent'),'Tag','PlotType');
if ~isempty(it)
itype = get(it(1),'value');
end

[a,b] = cmb.TrialRange(DATA);


if max(Counts(cell(cell>0))) > 1  %If list same cell twice, and probes are differnt, use templace scores to assign
[n,c] = Counts(cell);
cid = find(n > 1 & c > 0);
for j = 1:length(cid)
id = find(cell == c(cid(j)));
cmb.SetCellByTemplateScore(DATA, c(cid(j)), probe(id), [a:b],3);
end
return;
end
for j = 1:DATA.state.listlen
if cell(j) > 0
if quality(j) <= 1
DATA.CellList(cell(j),a:b) = 0;
else
DATA.CellList(cell(j),a:b) = probe(j);
end
if quality(j) > 8  %%Use recorded cluster quality
for e = 1:length(eid)
c = DATA.Expts{eid(e)}.Trials(1).Trial;
if eid(e) < length(DATA.Expts)
d = DATA.Expts{eid(e)+1}.Trials(1).Trial-1; %fill in gaps
else
d = DATA.Expts{eid(e)}.Trials(end).Trial; %fill in gaps
end
starts(e) = c;
ends(e) = d;
C = DATA.Expts{eid(e)}.Cluster{1,probe(j)};
if isfield(C,'quality') & C.quality > 0;
Q(e) = C.quality;
elseif isfield(C,'autocut') & C.autocut == 0
Q(e) =  10; %tag. Sort later
else
Q(e) = 0;
end
end
id = find(Q < 10);
if length(id) % at least one is set
last = id(end);
id = find(Q == 10)
isset = 0;

for e = 1:length(eid)
%if have set qual in this block, and it is set again later, fill in with the
%last set value
if Q(e) == 10 & isset & e < last
Q(e) = Q(e-1);
else
isset = 1;
end
if Q(e) == 0
DATA.CellQuality(cell(j),starts(e):ends(e)) =  Q(e);
DATA.CellList(cell(j),starts(e):ends(e)) =  0;

elseif Q(e) < 10
DATA.CellQuality(cell(j),starts(e):ends(e)) =  Q(e);
end
end
else
fprintf('No Cluster Quality data for Cell%d, Probe %d\n',cell(j),probe(j));
end
else
DATA.CellQuality(cell(j),a:b) = quality(j);
end
DATA.CellListCluster(cell(j),a:b) = cluster(j);
end
end

len = max([length(DATA.CellQuality)  length(DATA.CellList) length(DATA.CellListCluster)]);
if length(DATA.CellQuality) < len
lq = length(DATA.CellQuality);
DATA.CellQuality(:,lq+1:len)= 0;
end
if length(DATA.CellList) < len
lq = length(DATA.CellList);
DATA.CellList(:,lq+1:len) = 0;
end
if length(DATA.CellListCluster) < len
DATA.CellListCluster(:,len) = 0;
end

DATA.cellid = cell(1);
set(DATA.toplevel,'UserData',DATA);
%cmb.SaveCellList(DATA);
cmb.PlotCellList(DATA,'plotshapes');
cmb.SaveCellList(DATA,0,'temp');


