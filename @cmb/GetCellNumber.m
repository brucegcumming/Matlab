function GetCellNumber(DATA, eid, probe)
%was used by combine. now in cmb.GetCellNumber
if isfield(DATA,'CellList')
t = [DATA.Expts{eid}.Trials.Trial];
cid = mean(DATA.CellList(:,t));
[a,b] = min(diff(cid-probe));
end


