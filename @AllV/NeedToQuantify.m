function [needc, reason, C]  = NeedToQuantify(C, fullvname, varargin)


RefClusters = {};
CD = {};
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'clst')
        CD = varargin{j};
    elseif isfield(varargin{j},'RefClusters')
        RefClusters = varargin{j}.RefClusters;
    end
    j = j+1;
end

reason = 'Not set';
needc = 0;
if C.quick || C.trigdt == 4 || ...
        C.manual == 2 ||  C.dropi(3) == 0
    needc = 1;
    reason = 'Missing dropi';
end
if C.auto == 1 && length(RefClusters) >= C.probe(1) && isfield(RefClusters{C.probe(1)},'space')
    needc = 5;
    reason = 'Automatic Cut but RefCluster Exists';
end
if isfield(C,'needed') && C.needed
    needc = 1;
    reason = sprintf('Cluster has reason %d',C.needed);
elseif isfield(C,'manual') && C.manual ==2
    needc = 4;
    reason = 'Manual = 2';
elseif C.dropi(3) == 0
    needc  = 3;
    reason = 'Missing dropi';
elseif CheckSpkFile(C, fullvname) == 0
    needc = 2;
    reason = 'spkfile date mismatch';
end


for k = 1:length(C.next)
    if ~isfield(C.next{k},'space')
        ; %do nothing
    elseif (isfield(C.next{k},'quick') && C.next{k}.quick) ...
            || (isfield(C.next{k},'manual') && C.next{k}.manual == 2)...
            || C.next{k}.dropi(3) == 0
        needc = 1;
        reason = sprintf('Cluster %d need redo',k+1);
    elseif isfield(C.next{k},'space') %other checks
       if isfield(CD,'xy') && ~isfield(CD,'next') || length(CD.next) < k || ~isfield(CD.next{k},'xy')
           needc = 1;
           reason = sprintf('Cluster %d Details missing next or next.xy',k+1);
       end
    end
end
if C.trigdt == 4
    C.trigdt = 0;
end
