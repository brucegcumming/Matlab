function missedtrials = FindMissingTrials(DATA, Vall, ei)

missedtrials = [];
if isempty(Vall)
fprintf('Missing FullV Data for Expt %d\n',ei);
return;
end
blkend = (Vall.blkstart + Vall.blklen .* Vall.samper).*10000;
blkstart = Vall.blkstart .* 10000;

for j = 1:length(DATA.Expts{ei}.Trials)
tid = find(blkstart < DATA.Expts{ei}.Trials(j).Start(1) & blkend > DATA.Expts{ei}.Trials(j).End(end));
if isempty(tid)
missedtrials = [missedtrials DATA.Expts{ei}.Trials(j).id];
end
end


