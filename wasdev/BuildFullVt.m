function vt = BuildFullVt(FullV)
% builds vector of sample times from blkstart/blklen in a FullV struct
% vt = BuildFullVt(FullV)

    first = 1;
    vt = zeros(size(FullV.V));
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        vt(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        first = last+1;
    end
