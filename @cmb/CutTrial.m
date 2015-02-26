function CutTrial(a,b)
%DATA = cmb.combine('getstate');
DATA = GetDataFromFig(a);

t = DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).Trial;
if length(DATA.probelist) > 1 & isfield(DATA,'CellList')
cells = cmb.CellsSelected(DATA);
for j = 1:length(cells)
DATA.CellQuality(cells(j),t) = 3;
end
else
DATA.Expts{DATA.currentexpt(1)}.Trials(DATA.currenttrial).Trial = -abs(t);
end
it = findobj(DATA.svfig,'Tag','ChooseTrial');
set(it,'string',sprintf('%d|',[DATA.Expts{DATA.currentexpt(1)}.Trials.Trial]),'value',1);
set(DATA.toplevel,'UserData',DATA);
cmb.PlayOneTrial(DATA,a,1);

