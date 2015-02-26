function ReCombineByName(DATA, name, backup)

load(name);
if isfield(Expt.Header,'Combineids')
DATA.combineids = Expt.Header.Combineids;
else
DATA.combineids = [];
end
id = find(ismember(DATA.expid,DATA.combineids));
if ~isempty(id)
DATA.extype = 1;
DATA.spikelist = Expt.Header.Spikelist;
DATA.outname = name;
cmb.SetClusterCheck(DATA);

nstr = length(get(DATA.elst,'String'));
if max(id) > nstr
fprintf('%d is longer than expt list (%d)\n',max(id),nstr);
else
set(DATA.elst,'value',id);
end
[nExpt, DATA] = cmb.CombinePlot(DATA, 0,'ids',DATA.combineids);
end
CompareExpts(nExpt,Expt);
if backup
bname = strrep(name,'.mat','bak.mat');
if strcmp(name,bname)
fprintf('backup (%s) same as original(%s)\n');
return;
else
try save(bname,'Expt'); catch fprintf('ERROR saving %s\n',bname); return; end
end
end
Expt = nExpt;
try save(name,'Expt'); catch fprintf('ERROR saving %s\n',name); end

