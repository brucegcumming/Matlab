function xcorrhit(a,b, type, varargin);
DATA = GetDataFromFig(a);
GetFigure('CrossCorrelation');
if type == 1
if isfield(DATA.plot,'useprobe')
id = find(DATA.plot.useprobe)
if length(id)
cmb.PlotXcorr(DATA, id,3); %or [DATA.probes(id).probe] if not 1:nprobes
end
end
elseif type == 2
cmb.PlotAdjacentXcorrs(DATA,DATA.probelist,1);
end


