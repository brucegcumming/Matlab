function result  = GetStructVal(DATA, str, default)
%Get a value from a structure, checking fields exist.
if nargin > 2
    result = default;
else
    result = NaN;
end
    
f = split(str, '.')
X = DATA;
for j = 1:length(f)
    if isfield(X,f{j})
        X = X.(f{j});
    else
        return;
    end
end