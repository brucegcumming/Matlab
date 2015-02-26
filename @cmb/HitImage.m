function HitImage(a,b, type)
DATA = GetDataFromFig(a);
ax = get(a,'Parent');
xy = get(ax,'currentpoint');
l = get(ax,'Children');
tag = get(get(ax,'Parent'),'Tag');
ex = round(xy(1,2));
p = round(xy(1,1));
bt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});

zval = NaN;
xdata = [];
if length(xdata)  %%need to map data space into image pixel value
for j = 1:length(l)
a = get(l(j));
if isfield(a,'CData')
Z = get(l(j),'CData');
zval = Z(ex,p);
end
end
fprintf('Hit %.0f,%.0f %.3f\n',ex,p,zval);
end
if strcmp(type,'allexpt')
f = gcf;
X = get(gcf,'UserData');
id = X.sequence(ex);
Expt = DATA.AllExpts{id};
GetFigure(DATA.tag.dataplot);
hold off;
if Expt.Header.rc
[a,b,c] = PlotRC(Expt.plotres,'sdf');
if isfield(c,'figb')
figure(c.figb);
title(cmb.ProbeLabel(Expt));
end
if isfield(c,'figa')
figure(c.figa);
title(cmb.ProbeLabel(Expt));
end
else
PlotExpt(Expt);
end
figure(f);
end
offset = 0;
if ismember(type, [1 3]) % 3 = hit cell image - set cell#
C = DATA.AllClusters{ex}(p);
if type == 3
it = findobj('Tag','CellNumberId');
if DATA.CellList(ex,p,offset+1) > 0
DATA.currentcell = DATA.CellList(ex,p,offset+1);
set(it,'value',DATA.currentcell);
end
end
end
DATA.currentpoint = [ex p];
DATA = cmb.MarkCurrentCell(DATA);
set(DATA.toplevel,'UserData',DATA);


