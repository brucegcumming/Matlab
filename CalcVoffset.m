function [V, minV, maxV] = CalcVoffset(AllV, chspk, overlap)
%V = CalcVoffset(AllV, chspk, overlap)
%Cacluate voltage offsets so that channels don't overlap

if nargin < 3 || overlap < 0 || overlap > 1
    crit = 2;
else
    crit = 100.*overlap/2;
end
chspk = sort(chspk);
if max(chspk) > size(AllV,1)
    V(chspk) = 0;
    maxV = 0;
    minV = 0;
    return;
end
V(1) = 0;
for j = 1:length(chspk)
    yl(2) = prctile(max(AllV(chspk(j),:,:)),100-crit);
    yl(1) = prctile(min(AllV(chspk(j),:,:)),crit);
    maxV(chspk(j)) = yl(2);
    minV(chspk(j)) = yl(1);
end
pv = 0;
for j = 1:max(chspk)
    V(j) = pv-minV(j);    
    pv = maxV(j);
end
V = cumsum(V);
V(end:size(AllV,1)) = V(end);