function PlotISIHist(DATA, c)

Clusters = PC.CheckClusterLoaded(DATA,c(1));
C = Clusters{c(1)}{c(2)};
PC.SetFigure(DATA, 'ISIHist');
PlotISI(C);