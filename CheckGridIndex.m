function errs = CheckGridIndex(idx)

errs = {};
nerr = 0;
id = find(idx.expt == 0);
[a,b] = max(idx.nt(id));
if a > 1
    nerr = nerr+1;
    errs{nerr} = sprintf('%d unused trials in %s',a,idx.names{id(b)});
end

for j = 1:length(errs)
    fprintf('%s\n',errs{j});
end