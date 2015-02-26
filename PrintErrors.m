function PrintErrors(X, varargin)


errs = [];

if isfield(X,'errmsg')
    errs = X.errmsg;
elseif isfield(X,'errs')
    errs = X.errs;
end

for j = 1:length(errs)
    fprintf('%s\n'errs{j});
end
