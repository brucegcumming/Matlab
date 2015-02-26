function PlotSubspace(a,b, varargin)
if isfield(a, 'state') && isfield(a,'AllData')
DATA = a;
else
DATA = GetDataFromFig(a);
end
oldrc = DATA.plot.condenseRC;
DATA.plot.condenseRC = 0;
cmb.PlotCombined(DATA, DATA.Expt);
DATA.plot.condenseRC = oldrc;

