function probesource = FindProbeSources(DATA)
if ~isfield(DATA.probes,'source')
    probesource = ones(size(DATA.probelist));
else
    for j = 1:length(DATA.probelist)
        id = find([DATA.probes.probe] == DATA.probelist(j));
        probesource(j) = DATA.probes(id(1)).source;
    end
end

