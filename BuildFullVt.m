function vt = BuildFullVt(FullV)
% builds vector of sample times from blkstart/blklen in a FullV struct
% vt = BuildFullVt(FullV)

    first = 1;
    if isfield(FullV, 'V')
    vt = zeros(size(FullV.V));
    else
        vt = zeros(1,sum(FullV.blklen));
    end
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        vt(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        first = last+1;
    end
