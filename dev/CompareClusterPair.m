function res = CompareClusterPair(A,B, varargin)

aid = find(A.clst == 2);
bid = find(B.clst == 2);
at = A.times(aid);
bt = B.times(bid);

[a, ida, idb] = intersect(round(at.*1000),round(bt.*1000));
res.overlap = length(ida);
res.nspka = length(at);
res.nspkb = length(bt);
if A.auto == 1 && B.auto == 1
    res.auto = 4;
elseif A.auto == 1
    res.auto = 3;
elseif B.auto == 1
    res.auto = 2;
else
    res.auto = 1;
end
res.ctimea = A.savetime(1);
res.ctimeb = B.savetime(1);

xc = corrcoef(A.MeanSpike.ms,B.MeanSpike.ms);
res.spikexc = xc(1,2);
if size(A.MeanSpike.vdprime,1) == size(B.MeanSpike.vdprime,1) 
    xc = corrcoef(A.MeanSpike.vdprime,B.MeanSpike.vdprime);
    res.dpxc = xc(1,2);
else
    res.dpxc = NaN;
end

