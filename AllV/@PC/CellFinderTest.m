function CellFinderTest(DATA, type)    X = sum(DATA.CellList,3);X = isnan(X) + X > 0;[a,b] = find(X > 0);for j  = 1:length(a)    DATA = PC.TrackCluster(DATA,b(j),a(j),'UsePts',X);    cellps(j,:) = DATA.usepeaks;endplot(cellps);PC.MakeCellId(cellps);