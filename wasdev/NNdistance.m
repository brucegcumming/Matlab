function NNdistance(C, varargin)
%estimate cluster separation based on distance to neighbors
%
plottype = [1 1];
xytag = 'ClusterXY';
ftag = 'NearestNeighbor';
j = 1;
while j <= length(varargin)
    j = j+1;
end

aid = find(C.clst ==1);
bid = find(C.clst ==2);

A = C.xy(aid,1) + i.* C.xy(aid,2);
B = C.xy(bid,1) + i.* C.xy(bid,2);
AA = repmat(A,1,length(bid));
BB = repmat(B.',length(aid),1);
D = abs(AA-BB); %mutual distance matrix
[ab, iab] = min(D); %for each element of B, distance to nearest A
[ds, did] = sort(ab);
for j = 1:100
    bd = sort(D(:,did(j))); %distances to A for this shortest
    x = B(did(j));
    ax = abs(x-A);
    ax = sort(abs(x-B))';
    pts = 1:min([length(bd) length(ax) 1000]);
    y(j,:) = cumsum(bd(pts)')./pts;
    z(j,:) = cumsum(ax(pts))./pts;
end
if plottype(1)
GetFigure(ftag);
hold off;
plot(mean(y));
hold on;
plot(mean(z),'r');
end

if plottype(2)
    GetFigure(xytag);
    PlotND(C.xy,[],'.','idlist',C.clst)
    plot(C.xy(bid(did(1:100)),1),C.xy(bid(did(1:100)),2),'rx');
end
