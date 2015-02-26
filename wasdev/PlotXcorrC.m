function PlotCXcorr(C)
%Plots cross-correlations for a Clusters set
%wrapper that calls PlotAllXcorr from PlotClusters

DATA.name = 'Test';
DATA.toplevel = GetFigure('Xcorr');
for j = 1:length(C)
    cells(j).p = C{j}.probe(1);
    cells(j).cl = 1;
end

PlotAllXcorr(DATA, C, cells);
