function MarkProbes(a,b)
DATA = GetDataFromFig(a);
for j = 1:length(DATA.probelist)
it = findobj(get(a,'parent'),'Tag',sprintf('MarkProbe%d',j));
if length(it) == 1 && get(it,'value') > 0
DATA.markexp(DATA.currentexpt(1),j) = 1;
DATA.plot.useprobe(j) = 1;
else
DATA.plot.useprobe(j) = 0;
end
end
imagesc(DATA.markexp);
set(DATA.toplevel,'UserData',DATA);

