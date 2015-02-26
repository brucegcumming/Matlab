function AddComment(a,b)
DATA = GetDataFromFig(a);
n = 1;
if ~isfield(DATA,'Comments')
DATA.Comments = PlotComments(DATA.datadir,'parent',DATA.toplevel,'noninteractive');
end
DATA.Comments = PlotComments(DATA.Comments, 'addhidden', get(a,'string'),DATA);




