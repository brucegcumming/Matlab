function SetCellCompare(a,b, cellid)DATA = GetDataFromFig(a);D = get(gcf,'UserData');if ~CellIsEmpty(DATA.CellDetails,'MeanSpike', cellid)    DATA.comparecell(cellid) = ~DATA.comparecell(cellid);    if DATA.comparecell(cellid)        AllV.AddCellMean(DATA,cellid);    elseif strcmp(D.plottype,'QuickSpikes')s %redraw without        AllV.QuickSpks(DATA,1000);    endendset(DATA.toplevel,'UserData',DATA);    