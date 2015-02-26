function  DATA = MainMenu(a,b,type)
DATA = GetDataFromFig(a);


if type == 1
DATA.plot.useprobe = ~DATA.plot.useprobe;
for j = 1:length(DATA.plot.useprobe)
it = findobj(DATA.toplevel,'Tag', sprintf('UseProbe%d',j));
set(it,'value',DATA.plot.useprobe(j));
end
end

set(DATA.toplevel,'UserData',DATA);


