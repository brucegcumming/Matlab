function DATA = AddTimesToFullV(DATA)
% DATA = AddTimesToFullV(DATA)
%calculate timestamp values for each point in a FullV structre
%takes longer thanyou think, so worth precalculating. 

if isfield(DATA,'blklen')
    first = 1;
    vt(sum(length(DATA.blklen))) = 0;
    for j = 1:length(DATA.blklen)
        last = first+DATA.blklen(j)-1;
        vt(first:last) = DATA.blkstart(j)+[1:DATA.blklen(j)].*DATA.samper;
        first = last+1;
    end
else
    vt= [1:size(Vall.V,2)] .* DATA.interval;
end
DATA.t = vt;