function SetTrialRange(a,b, type)
DATA = GetDataFromFig(a);

if type == 2
DATA.firsttrial = DATA.currenttrial;
elseif type ==3
DATA.firsttrial = 0;
DATA.lasttrial = 0;
elseif type ==1
DATA.lasttrial = DATA.currenttrial;
end

set(DATA.toplevel,'UserData',DATA);
cmb.CheckSpoolButton(DATA);
cmb.Recmb.PlotCellList(DATA);



