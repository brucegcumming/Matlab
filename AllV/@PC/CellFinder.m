function CellFinder(DATA)        mid = DATA.mahaltype;%First make a smoothed map of mahal distances, and find local maximaif DATA.mahaltype == 5    X = squeeze(DATA.mahal(:,:,1));else    X = squeeze(DATA.mahal(:,:,mid));end[a,b,G] = Gauss2D(0.8,[-3:1:3]);Y = conv2(X,G,'same');xm = diff(sign(diff(Y)));ym = diff(sign(diff(Y,1,2)),1,2);[a,b] = find(xm(:,2:end-1) < 0 & ym(2:end-1,:) < 0);mxs = xm(:,2:end-1) + ym(2:end-1,:);%add in any squares where mahal distance > 3[c, d] = find(DATA.mahal(:,:,1) > 3 | DATA.mahal(:,:,2) > 3);a = cat(1,a+1 ,c);b = cat(1, b+1, d);for j  = 1:length(a)    DATA = PC.TrackCluster(DATA,b(j),a(j));    cellps(j,:) = DATA.usepeaks;endplot(cellps);PC.MakeCellId(cellps);