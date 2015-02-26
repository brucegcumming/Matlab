function PlotTemplates(a,b, varargin)
DATA = GetDataFromFig(a);
GetFigure('CombinerSpikeV');
plot(DATA.Templates');
for j = 1:min([size(DATA.Templates,1) size(DATA.TemplateInfo,2)])
if isfield(DATA.TemplateInfo(j),'probe')
labels{j} = sprintf('P%dE%d',DATA.TemplateInfo(j).probe,DATA.TemplateInfo(j).exid);
end
end
legend(labels,'position',[0.8 0.2 0.2 0.2]);

