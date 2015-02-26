function S  = AddField(S, f, default)
%S  = AddField(S, f, default) Add a field to structure, IF its not already
%there. Assign to default or zero if nargin  == 2
if nargin == 2
    default = 0;
end

if iscell(f)
    for j = 1:length(f)
        if ~isfield(S,f{j})
            S.(f{j}) = default;
        end
    end
elseif ~isfield(S,f)
    S.(f) = default;
end