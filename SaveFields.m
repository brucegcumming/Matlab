function SaveFields(name, X, varargin)
%SaveFields(name, X) save the fields of X in file name
% So X = load(name); SaveFields(name, X) writes out only the variables
% loaded

f = fields(X);
for j = 1:length(f)
    eval(sprintf('%s=X.%s',f{j},f{j}));
end
save(name,f);