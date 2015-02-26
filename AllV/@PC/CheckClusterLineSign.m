function DATA = CheckClusterLineSign(DATA)    Clusters = getappdata(DATA.toplevel,'Clusters');Expts = getappdata(DATA.toplevel,'Expts');sgn = CellToMat(Clusters,'sign');shp = CellToMat(Clusters,'shape');angles = CellToMat(Clusters,'angle');if isfield(DATA,'tagged')    oldsum = sum(DATA.tagged(:) ==2);    DATA.tagged(DATA.tagged ==2) = 0;else    oldsum = 0;endp = find(sum(sgn < 0) & sum(sgn > 0));for j =1:length(p)    [a,b] = Counts(sgn(:,p(j)));    if a(b == -1) > a(b == 1)        id = find(sgn(:,p(j)) == 1);        gid = find(sgn(:,p(j)) == -1);    else        id = find(sgn(:,p(j)) == -1);        gid = find(sgn(:,p(j)) == 1);    end    shape = prctile(shp(gid,p(j)),50);    angle = prctile(angles(gid,p(j)),50);    for k = 1:length(id)        if shp(id(k),p(j)) == shape && cos(angle - angles(id(k),p(j))) > 0            mycprintf('blue','E%d P%d Cluster Line Inverted\n',id(k),p(j));            DATA.tagged(id(k),p(j)) = 2;        end    end     endif isfield(DATA,'tagged') && sum(DATA.tagged(:) == 2) || oldsum > 0    DATA.markcell.tagged = 1;    DATA = PC.PlotCellList(DATA);    set(DATA.toplevel,'UserData',DATA);end