function [x, details] = FindExcludedTrials(DATA,e,p, cluster, C)%FindExcludedTrials(DATA,e,p, cluster, C) e and p are indicesdetails.xid = [];    x = [];    if ~isfield(DATA,'trialids') || cluster < 1        return;    end    if isfield(C,'missingtrials')        x = find(ismember(DATA.trialids{e},C.missingtrials));    end    if isfield(C,'excludetrialids')        y = find(ismember(DATA.trialids{e},C.excludetrialids));        x = [x y];        details.xid = C.excludetrialids;    end    details.nt = length(DATA.trialids{e});    details.fraction = 0;        if ~isfield(DATA.CellDetails,'excludetrials')        DATA.CellDetails.excludetrials=[];    end    XC = DATA.CellDetails.excludetrials;              if size(XC,1) < e || size(XC,2) < p || size(XC,3) < cluster        return;    end    y = find(ismember(DATA.trialids{e},XC{e,p, cluster}));    x = [x y];    details.xid = DATA.trialids{e}(x);    details.fraction = length(x)./details.nt;    return;        %Old stuff    if p <= 0 || size(DATA.CellList,2) < p || size(DATA.CellList,1) < e        return;    end    cellid = DATA.CellList(e,p, cluster);    if cellid > 0 && length(XC) >= cellid && length(XC{cellid}) >= e        x = [x XC{e, cellid, cluster}];        details.fraction = length(x)./details.nt;    end     