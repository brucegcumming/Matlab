function exitallv(src, evnt)    DATA = get(src,'UserData');    if 0    if isfield(DATA,'tag')        f = fields(DATA.tag);        for j = 1:length(f)            if ~strcmp(f{j},'top')                CloseTag(DATA.tag.(f{j}));            end        end    end    end    AllV.MiscMenu(DATA,[],'savelastconfig');    AllV.MiscMenu(DATA,[],'savelastlayout');    CloseChildren(src);    try        fclose(DATA.logfid);    end    fids = fopen('all');    for j = 1:length(fids)        fname = fopen(fids(j));        if regexp(fname,'ClusterLog.*.txt')            fclose(fids(j));        end    end    delete(src);    