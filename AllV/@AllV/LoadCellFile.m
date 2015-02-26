function DATA = LoadCellFile(DATA)    cellfile = BuildFileName(DATA.name,'celllist');    if exist(cellfile,'file')        load(cellfile);        DATA.CellList = CellList;        DATA.CellDetails = CellDetails;        DATA.CellChanges = CellChanges;        DATA = AllV.LoadComments(DATA);        if isfield(DATA,'tagged')            [a,b] = find(DATA.tagged > 0);            id = find(DATA.CellDetails.exptids(a) == DATA.exptno);            if ~isempty(id)                DATA.TaggedProbes = DATA.tagged(a(id(1)),:);                for j = 1:length(id)                    fprintf('Expt %d Probe %d Tagged %d\n',DATA.exptno,b(id),DATA.TaggedProbes(b(id(j))));                end                AllV.ShowTaggedProbes(DATA);            end        end        nc = max(DATA.CellList(:));        n = length(DATA.comparecell);        if nc > n            DATA.comparecell(n+1,nc) = 0;        end    end                