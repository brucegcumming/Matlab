function CombineRFData(varargin)
%rflist = CombineRFData(rflist, fits)
%Combines together RF date from ulf files (rflist) with fitted RFs
%(for multiple probes in array data) to produce one rflist for PlotMat

fits = []
rfs = [];
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        if length(CellToMat(varargin{j},'proberf')) > 0
            fits = varargin{j};
        elseif length(CellToMat(varargin{j},'depth')) > 0
            rfs = varargin{j};
        end
    end
    j = j+1;
end

for j = 1:length(rfs)
    [a, rfnames{j}] = fileparts(rfs{j}.name);
end
for j = 1:length(fits)
    if isfield(fits{j},'dirname');
        fitnames{j} = fileparts(fits{j}.dirname);
        id = find(strcmp(fitnames{j},rfnames));
    end
end

