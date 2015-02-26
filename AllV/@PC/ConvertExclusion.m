function DATA = ConvertExclusion(DATA)    %excludetrials and trialids are indexed by exptid.%i.e. match rows, not expts    Expts = getappdata(DATA.toplevel,'Expts');    if isempty(Expts) && size(DATA.CellDetails.excludetrials,1) < size(DATA.CellList,1)        DATA.CellDetails.excludetrials{size(DATA.CellList,1),DATA.nprobes,6} = [];        return;    end    exlist = GetExptNumber(Expts);    expid = DATA.CellDetails.exptids;    nxid = length(expid);    if size(DATA.CellDetails.excludetrials,1) < nxid        DATA.CellDetails.excludetrials{nxid,DATA.nprobes,6} = [];    end    if isfield(DATA.CellDetails,'version') & DATA.CellDetails.version > 1.1        PC.CheckExclusion(DATA);        return;    end    XC = DATA.CellDetails.excludetrials;    for j = 1:size(DATA.CellDetails.excludetrials,1)        e = find(exlist == expid(j));         for k = 1:size(DATA.CellDetails.excludetrials,2)            for c = 1:size(DATA.CellDetails.excludetrials,3)                X = DATA.CellDetails.excludetrials{j,k,c};                if length(X)                    if max(X) > length(Expts{e}.Trials)                        fprintf('Too Many Excluded E%dP%dC%d: %d vs %d\n',j,k,c,length(X),length(Expts{e}.Trials));                        X = X(X < length(Expts{e}.Trials));                        ids = [Expts{e}.Trials(X).id];                        XC{k,j,c} = ids;                    else                        ids = [Expts{e}.Trials(X).id];                        XC{k,j,c} = ids;                    end                end            end        end    end    n =1;    DATA.CellDetails.version = DATA.version;    if size(DATA.CellDetails.excludetrials,1)    DATA.CellDetails.excludetrials = XC;    end%    set(DATA.toplevel,'UserData',DATA);    