function CompareIdx(a,b)

f = fields(a.Trials);
for j = 1:length(f)
    if iscell(a.Trials.(f{j}))
        for k = 1:length(a.Trials.(f{j}))
            dl(k) = length(a.Trials.(f{j})(k))-length(b.Trials.(f{j})(k));
        end
        if sum(abs(dl))
            fprintf('Size mismatch for %s\n',f{j});
        end
    else
        id = find(~isnan(a.Trials.(f{j})));
        diffs = a.Trials.(f{j})(id)-b.Trials.(f{j})(id);
        if sum(abs(diffs(:)))
            fprintf('Mismatch for %s\n',f{j});
        end
    end
end