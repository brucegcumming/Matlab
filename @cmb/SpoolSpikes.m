function SpoolSpikes(a,b)
DATA = GetDataFromFig(a);
DATA.spooling = 1;
if isfigure(a) & Bstrmatch(get(a,'Tag'),'SpoolSpikes') & DATA.currenttrial > 1
DATA.spooling = 2; %spool from current spike to end
end
if DATA.state.nospikes == 2
DATA = cmb.LoadSpikes(DATA,DATA.currentexpt(1));
end
DATA = cmb.PlaySpikes(DATA,DATA.currentexpt(1));
set(DATA.toplevel,'UserData',DATA);

