function PlotEig(DATA)
if isfield(DATA.AllData,'EigVec')
GetFigure('EigenVectors');
plot(DATA.AllData.EigVec(:,end-5:end));
end

