function DATA = CheckCellExptList(DATA)
if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'exptids') && isfield(DATA,'exabsid')
for j = 1:length(DATA.exabsid)
a = find(DATA.CellDetails.exptids == DATA.exptnos(DATA.exabsid(j)));
if length(a) == 1
DATA.cellexid(j) = a;
else
DATA.cellexid(j) = 0;
end
end
end


