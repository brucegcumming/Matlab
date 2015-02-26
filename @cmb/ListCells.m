function ListCells(a,b)
onoff = {'off' 'on'};
DATA = GetDataFromFig(a);
DATA.listbycell = ~DATA.listbycell;
set(a,'Checked',onoff{DATA.listbycell+1});
eid = get(DATA.clst,'value');
cmb.ListSubExpts(DATA,eid);
pit = findobj(DATA.toplevel,'Tag','ProbeId');
pn = get(pit,'value');
if DATA.listbycell
DATA.state.recount = 1;
cells = unique(DATA.CellList(:));
cells = cells(cells > 0);
for j = 1:length(cells)
end
set(pit,'string',num2str(cells));
setappdata(pit,'probelist',cells);
if pn > length(cells)
set(pit,'value',1);
end
%        DATA = cmb.cListExpts(DATA, DATA.Expts);
else
set(pit,'string',DATA.probenames);
setappdata(pit,'probelist',DATA.probelist);
eid = get(DATA.clst,'value');
DATA = cmb.ListSubExpts(DATA,eid);
DATA = cmb.SetProbe(DATA,1);
set(pit,'value',DATA.probe);
end
cmb.SetGui(DATA);
set(DATA.toplevel,'UserData',DATA);

