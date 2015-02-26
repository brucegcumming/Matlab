function [CellId, details] = MakeCellId(cellps)     %cellps = Nstarts X N expts array, each cell being the best matching probe %for that expt and that starting point. for j = 1:size(cellps,2)    cim(:,j) = hist(cellps(:,j),[0:24]);end imagesc(cim(2:end,:));%cmid is frequency with which probes are selected, for each exptid = find(cim > 0);[a,cid] = sort(cim(id),'descend');cid = id(cid); %nonzero elements in descending order[a,b] = ind2sub(size(cim),cid);CellId = zeros(size(cim,2),size(cim,1)-1);CellCounts = CellId;%b is list of expts%a is list of probes +1 - if a == 1 means probe was 0  nc = 0;for j = 1:length(a)    if a(j) > 1%find all runs that have matching probe at this expt    id = find(cellps(:,b(j)) == a(j)-1);%c = hist(cellps(id,:),[0:24]);[d,e] = max(c(2:end,:));did = find(d>0);    x = sub2ind(size(CellId),did,e(did));    eid = find(d(did) > CellCounts(x));    x = x(eid);    if length(x) > 0%check to see if this is already defined        if sum(CellId(x)) == 0        oldcell = 0;    else        oldcell =  mode(CellId(x(CellId(x) > 0)));        if sum(CellId(x) > 0 & CellId(x) ~= oldcell) > sum(CellId(x) == oldcell)/3            oldcell = 0;        end    end    if oldcell == 0         nc = nc+1;        CellId(x) = nc;    else        CellId(x) = oldcell;    end    CellCounts(x) = d(did(eid));    end    endendfor j = 1:size(CellId,1)    [a,b] = Counts(CellId(j,:));    id = find(a >1 & b > 0);    for k = 1:length(id)        aid = find(CellId(j,:) == b(id(k)));        [a,c] = max(CellCounts(CellId(j,aid)));        CellId(j,aid) = 0;        CellId(j,aid(c)) = b(id(k));    endendhold off;imagesc(CellId);fprintf('%d Cells\n',length(unique(CellId(:)))-1);    