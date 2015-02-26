function AddOneCellToList(caller,b,varargin)
DATA = GetDataFromFig(caller); 
p = DATA.probe;
cl = DATA.currentcluster;
dp = DATA.cluster{cl,p}.dprime;
it = findobj(get(caller,'Parent'),'Tag','ClusterQuality');
q = get(it,'value');
if q == 10
if dp > 5
q = 6;
elseif dp > 3.5
q = 5;
elseif dp > 2.8
q = 4;
elseif dp > 2.2
q = 3;
else
q = 2;
end
end
it = findobj(get(caller,'Parent'),'Tag','AddOneCellToList');
cell = get(it,'value');
[a,b] = cmb.TrialRange(DATA);
if DATA.state.online
a = a + DATA.Expts{DATA.currentexpt(1)}.gui.firsttrial-1;
b = b + DATA.Expts{DATA.currentexpt(1)}.gui.firsttrial-1;
end
DATA.CellList(cell,a:b) = p;
DATA.CellQuality(cell,a:b) = q;
DATA.CellListCluster(cell,a:b) = DATA.currentcluster;
DATA.cellid = cell(1);
set(DATA.toplevel,'UserData',DATA);
%cmb.SaveCellList(DATA);
GetFigure(DATA.tag.celllist);
cmb.PlotCellList(DATA,'qual');
cmb.SaveCellList(DATA,0,'temp');


