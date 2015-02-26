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

CheckSpkFile(C, fullvname);

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
