function [eid, cid, shifts] = FindNoCell(DATA, probe, expts)
%find expt/probes wiht no cell defined
%arg 3 is ignored at the moment
if isfield(DATA,'CellList')
    cellid = sum(DATA.CellList,3);
    eid = [];
    cid = [];
    if isfield(DATA.CellDetails,'probedrift') && DATA.state.fixdrift
        drift = round(DATA.CellDetails.probedrift);
        np = size(cellid,2);
        for j = 1:size(cellid,1)
            p = probe + drift(j);
            if p > 0 && p <= np && cellid(j,p) == 0
                eid(end+1) = j;
                cid(length(eid)) = p;
            end
        end
    else
        [eid] = find(cellid(:,probe) == 0);
        cid = ones(size(eid)).*probe;
    end
else
    eid = DATA.expid;
    cid = DATA.probelist;
end


