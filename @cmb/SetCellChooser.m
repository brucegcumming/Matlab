function SetCellChooser(DATA)
F = findobj('type','figure','Tag',DATA.tag.allexpts);
if ~isempty(F)
    it = findobj(F,'Tag','ExptCellChoose');
else
    it = findobj('Tag','ExptCellChoose');
end
if isempty(it)
    return;
end
%length(it) can be > if Keep has been used.  

idn = 1;
if length(it) > 1 %Find Right Figure
    
end
delete(get(it(idn),'Children'));
AllCellRes = getappdata(DATA.toplevel,'AllCellRes');

for j = 1:length(AllCellRes)
if AllCellRes(j).cellid > 0
uimenu(it(idn),'label',sprintf('Cell%d (%.1f)',AllCellRes(j).cellid,AllCellRes(j).Header.probe),'callback',{@cmb.AllCellPlots, j});
else
uimenu(it(idn),'label',sprintf('P%dmu',AllCellRes(j).Header.probe),'callback',{@cmb.AllCellPlots, j});
end
end

