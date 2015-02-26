function FullVData = MakeFullVInfo(Vall)
%FullVData = MakeFullVInfo(Vall)
%Builds struct with summary infor about fullV, but none of the big stuff

f = {'submean', 'usealltrials' 'builddate' 'samper' 'missedtrials' 'blklen'...
    'blkstart' 'highpass' 'coilnoiseratio'};
for j = 1:length(f)
    if isfield(Vall,f{j})
        FullVData.(f{j}) = Vall.(f{j});
    end
end
