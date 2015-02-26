function SaveCellList(a,b, varargin)


if isstruct(a)
DATA = a;
else
DATA = GetDataFromFig(a);
end

if isfield(DATA,'CellList')
CellList = DATA.CellList;
CellQuality = DATA.CellQuality;
CellListCluster = DATA.CellListCluster;
else
CellList = [];
CellQuality = [];
CellListCluster = [];
end

Templates = DATA.Templates;
TemplateInfo = DATA.TemplateInfo;
if ~isfield(DATA,'cellfile')
DATA.cellfile = strrep(DATA.datafilename,'.mat','.cells.mat');
end

cellfile = DATA.cellfile;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'temp',4)
cellfile = strrep(DATA.cellfile,'cells','celltmp');
end
j = j+1;
end

save(cellfile,'CellList','CellQuality', 'Templates','CellListCluster','TemplateInfo');
set(DATA.toplevel,'UserData',DATA);

