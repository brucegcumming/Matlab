function id = FullVTimes2id(Vall, t)
% id = VallTimes2id(Vall, t) 
%return sample points matching list of times, t, in a FullV structure
id = [];
for j = 1:length(Vall.blkstart)
    tid = find(t >= Vall.blkstart(j) & t <= Vall.t(sum(Vall.blklen(1:j))));
    tid = round((t(tid)-Vall.blkstart(j))./Vall.samper);
    id = [id tid];
end