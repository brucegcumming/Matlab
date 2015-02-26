function AllCellKeyPressed(src,ks, op)
DATA = GetDataFromFig(src);
AllCellRes = getappdata(DATA.toplevel,'AllCellRes');
P =  getappdata(DATA.toplevel,'AllCellPlot');
if strcmp(ks.Key,'downarrow') && P.currentcell < P.nplots
    P = cmb.NextPlot(P,AllCellRes,0);
elseif strcmp(ks.Key,'uparrow') && P.currentcell > 1
    P = cmb.NextPlot(P,AllCellRes,-1);
end

if isfield(DATA,'fig') && isfield(DATA.fig,'CombinerAllExpts') && gcf == DATA.fig.CombinerAllExpts;
    setappdata(gcf,'ParentFigure',DATA.toplevel);
    set(gcf,'KeyPressFcn',@cmb.AllCellKeyPressed);
end

setappdata(DATA.toplevel,'AllCellPlot',P);
set(0,'CurrentFigure',gcbf);

