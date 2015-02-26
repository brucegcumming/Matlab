function RePlotCellList(a,b)
if isstruct(a)
DATA = a;
if ~isempty(findobj('Tag',DATA.tag.celllist))
GetFigure(DATA.tag.celllist);
end
else
DATA = GetDataFromFig(a);
itype =  1;
it = findobj(get(a,'parent'),'Tag','PlotType');
if ~isempty(it)
itype = get(it(1),'value');
end
DATA.plot.cellplot = itype;
end

if isempty(findobj('Tag',DATA.tag.celllist))
return;
end
cmb.PlotCellList(DATA);
set(DATA.toplevel,'UserData',DATA);


