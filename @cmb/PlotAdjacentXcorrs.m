function PlotAdjacentXcorrs(DATA, probes, type)
GetFigure('CrossCorrelation');
np = length(probes);
[nr, nc] = NSubplots(np-1);
for j = 1:np-1
subplot(nr,nc,j);
cmb.CalcXcorr(DATA, DATA.currentexpt(1),probes(j),probes(j+1));
set(gca,'xtick',[],'ytick',[])
title(sprintf('%d,%d',j,j+1));
end


