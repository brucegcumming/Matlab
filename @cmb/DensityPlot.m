function DensityPlot(a,b)
global mousept;

%DATA = cmb.combine('getstate');
DATA = GetDataFromFig(a);

if isfield(DATA,'spklist') & length(DATA.spklist) > 10
expspks = DATA.spklist;
else
expspks = DATA.spkrange(1):DATA.spkrange(2);
end


figure(DATA.xyfig);
if DATA.densityplot
DATA.densityplot = 0;
hold off;
DATA = cmb.DrawXYPlot(DATA,expspks);
%    ClearMouse;
set(DATA.toplevel,'UserData',DATA);
cmb.SetGui(DATA);
return;
end

[x,y,z] = CombineCalcDensity(DATA, expspks, 2);

erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
toc
pcolor(x,y,z);
shading('interp')
hold on;
if 0 %code for trackin peak ridge. May help track clusters/auto cut
[a,b] = max(z);
py = smooth(y(b,1),3,'gauss');
plot(x(1,:), py);
base = 0.1;
for j = length(py):-1:1
dip(j) = a(j)./(base+max(a(j:end)).^2);
end
plot(dip);
end
DATA = cmb.DrawClusters(DATA, DATA.cluster, 0);
DATA.densityplot = 1;
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
x = caxis;
set(it,'string',sprintf('%.2f',x(2)));
cmb.SetGui(DATA);
set(DATA.toplevel,'UserData',DATA);

