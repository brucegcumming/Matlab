function lProbeClusters(a,b)
if isstruct(a)
DATA = a;
else
DATA = GetDataFromFig(a);
end
for j = 1:length(DATA.Expts)
DATA.allexp = j;
cmb.PlotAllProbeXY(DATA);
drawnow;
end

