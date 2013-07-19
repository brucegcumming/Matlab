function res = iscluster(C,c,p)
    if ~iscell(C) | isempty(C) | size(C,1) < c | size(C,2) < p | isempty(C{c,p})
        res = 0;
    elseif ~isfield(C{c,p},'params') | ~isfield(C{c,p},'x')
        res = 0;
    elseif ~isfield(C{c,p},'nspk') || C{c,p}.nspk == 0
        res = 2;
    else
        res = 1;
    end
