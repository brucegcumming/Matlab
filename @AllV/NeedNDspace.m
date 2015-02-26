function yes = NeedNDspace(C, varargin)
% yes = NeedNDspace(C, varargin) determine if need to make ND auto space
% for classifying
yes = false;
if C.space(1) == 6
    yes = true;
end
if isfield(C,'next')
for j = 1:length(C.next)
    if isfield(C.next{j},'space') && C.next{j}.space(1) == 6
        yes = true;
    end
end
end