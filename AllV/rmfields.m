function S = rmfields(S, varargin)
%rmfields(S, a,b,c....)    removes a list of fields from S
%   just calls rmfield, building a cell string array if necessary.
%   Also checks that the field exits to avoid errors
% if S is a cell array, does it for each element of S

f = {};
if iscell(S)
    for j = 1:length(S)
        S{j} = rmfields(S{j},varargin{:});
    end
    return;
end

j = 1;
while j <= length(varargin)
    if iscellstr(varargin{j})
        for k = 1:length(varargin{j})
            if isfield(S, varargin{j}{k})
                f= {f{:}, varargin{j}{k}};
            end
        end
    elseif isfield(S, varargin{j})
        f= {f{:}, varargin{j}};
    end
    j = j+1;
end
    
if ~isempty(f)
    S = rmfield(S, f);
end