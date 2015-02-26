function ClearCell(name, expt, cell)

cellfile = BuildFileName(name,'celllist');
mod = 0;
if exist(cellfile)
    X = load(cellfile);
    CD = X.CellDetails;
    eid = find(ismember(CD.exptids,expt));
    for j = 1:length(eid)
        [a,b] = find(X.CellList(eid(j),:,:) == cell);
        if ~isempty(a)
            X.CellList(eid(j),b,a) = 0;
            mod = a;
        end
    end
    if mod
        BackupFile(cellfile);
        save(cellfile,'-struct',X);
    end
end
