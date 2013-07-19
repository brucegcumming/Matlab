function ediffs = CompareGridIdx(aidx,bidx)

if length(aidx.expt) == length(bidx.expt)
    id = find(aidx.expt ~= bidx.expt);
    ediffs = id;
    for j = 1:length(id)
        fprintf('%s %d %s %d\n',aidx.names{id(j)},aidx.expt(id(j)),bidx.names{id(j)},bidx.expt(id(j)));
    end
    gid = find(aidx.expt == bidx.expt);
    fprintf('Tdiffs%s\n',sprintf(' %.1f',aidx.tdiff(gid)-bidx.tdiff(gid)));
else
    if length(aidx.expt) > length(bidx.expt)
        ediffs = find(aidx.expt(1:length(bidx.expt)) ~= bidx.expt);
        ediffs = [ediffs aidx.expt(length(bidx.expt)+1:end)];
    else
        ediffs = find(aidx.expt) ~= bidx.expt(1:length(aidx.expt));
        ediffs = [ediffs bidx.expt(length(aidx.expt)+1:end)];
    end
end