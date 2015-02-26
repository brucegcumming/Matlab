 function X = CheckClusters(DATA, type)    Clusters = getappdata(DATA.toplevel,'Clusters');Expts = getappdata(DATA.toplevel,'Expts');if iscell(type)    for j = 1:length(type)        X{j} = PC.CheckClusters(DATA,type{j});    end    return;endX.errs = {}; for j = 1:length(Clusters)     if ~isempty(Clusters{j})         e = floor(Clusters{j}{1}.exptno);         eid = Clusters{j}{1}.exptid;     if isempty(Expts{e})         nt = 0;     else        nt = length(Expts{e}.Trials);     end     for k = 1:length(Clusters{j})         C = Clusters{j}{k};         if strcmp(type,'exclusions') && isfield(C,'excludetrialids')             if isfield(C,'restricttimerange')                 fprintf('Row%dP%d %d/%d excluded Trials (%.1f-%.1f)\n',j,k,length(C.excludetrialids),nt,C.restricttimerange(1),C.restricttimerange(2));             else                 fprintf('Row%dP%d %d/%d excluded Trials\n',j,k,length(C.excludetrialids),nt);             end         elseif strcmp(type,'nclusters')             if isfield(C,'next') & length(C.next) > 1                 fprintf('Row%dP%d %d clusters\n',j,k,length(C.next)+1);             end         elseif strcmp(type,'fitspace')             PC.CheckFitDim(C);         elseif strcmp(type,'fittimes')             if isfield(C,'savetime') && C.savetime(end)-C.savetime(1) < -0.1                 cprintf('blue','Expt%d(Row%d)P%d Fit is old\n',DATA.exptid(j),j,k);             end         elseif strcmp(type,'celldefined')             [a, cells, clid] = PC.isacell(DATA,eid,k);              for c = 1:length(clid)                 if clid(c) > 1 && (length(C.next) < clid(c)-1 || isempty(C.next{clid(c)-1}))                    X =  AddError(X,'-show','Expt%d(Row%d)P%d Cluster %d is Cell %d but its empty\n',DATA.exptid(j),j,k,clid(c),cells(c));                 end             end         end     end     end end