function [eid, cid, clid] = FindCell(DATA, cellid)
%[eid, cid, clid] = FindCell(DATA, cellid)
%find rows, columns, clusters in CellList that are cell cellid
id = find(DATA.CellList == cellid);
[eid, cid, clid] = ind2sub(size(DATA.CellList),id);
[eid, id] = sort(eid);
cid = cid(id);
clid = clid(id);




