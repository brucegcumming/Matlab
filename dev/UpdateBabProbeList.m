function Array = UpdateBadProbeList(dirname, Clusters)
%FInds probes markded as bad in Clusters, and Adds these to the ArrayConfig
%File

Array = GetArrayConfig(dirname);
if iscell(Clusters)
    marked = CellToMat(Clusters,'marked');
    marked(isnan(marked)) = 0;
    ismark = find(sum(marked)>0);
    for j = 1:length(ismark)
        p = ismark(j);
        mark = median(marked(:,p));
        if ismember(mark,[3 5]) %bad or duplicate
        if p > length(Array.badprobes) || mark ~= Array.badprobes(p)
            if mark == 3
                Array = GetArrayConfig(dirname,'markbad',p);
                str = 'Bad';
            elseif mark == 5
                Array = GetArrayConfig(dirname,'markdup',p);
                str = 'Duplicate';
            end
            fprintf('Marking probe %d %s\n',p,str)
        end
        end
    end
end