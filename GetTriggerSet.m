function [T, details] = GetTriggerSet(C, varargin)
%List  which cluster belongs to which triggerset

verbose = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'warn',4)
        verbose =1;
    end
    j = j+1;
end
T(1,1) = C.triggerset;
details.duplicates = [];
for c = 1:length(C.next)
    if isfield(C.next{c},'triggerset')
        T(1,c+1) = C.next{c}.triggerset;
    else
        T(1,c+1) = NaN;
    end
end
if ~isfield(C,'trigset')
   return;
end
nc = size(T,2);
for k = 1:length(C.trigset)
B = C.trigset{k};
if B.triggerset > 0
    details.duplicates(k,c+1) = k;
    if verbose
        cprintf('blue','Cluster 1/%d is in Trigger Sets %d and %d\n',nc, T(1,1),k);
    end
end
for c = 1:length(B.next)
    if isfield(B.next{c},'triggerset')
        if T(1,c+1) >= 0
            details.duplicates(k,c+1) = 1;
            if verbose
                cprintf('blue','Cluster %d/%d is in Trigger Sets %d and %d\n',c+1,nc,T(1,c+1),k);
            end
        end
        T(1,c+1) = B.next{c}.triggerset;
        nc = size(T,2);
    end
end
end