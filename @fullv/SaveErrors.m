function errs = SaveErrors(FullV, varargin)
%take any errors in a FullV file and save them into a single
%file in each folder logging all FullV Errros

X.Errors= {};
silent = 0;
if ~isfield(FullV,'errmsg') || isempty(FullV.errmsg) %nothing to do
    return;
end
strs = {};

name = GetName(FullV,'folder');
if isdir(name)
    ename = [name '/Errors.mat'];
    if exist(ename)
        X = load(ename);
        strs = CellToMat(X.Errors,'s');
    end
end

newerrs = 0;
if length(FullV.errdata) < length(FullV.errmsg)
    FullV.errdata{length(FullV.errmsg)} = [];
end
for j = 1:length(FullV.errmsg)
    if isempty(X.Errors)
        X.Errors(1).s = FullV.errmsg{j};
        X.Errors(1) = CopyFields(X.Errors(1),FullV.errdata{j});
        newerrs = newerrs+1;
    elseif sum(strcmp(FullV.errmsg{j},strs)) == 0 %new error
        X.Errors(end+1).s = FullV.errmsg{j};
        X.Errors = CopySFields(X.Errors,length(X.Errors),FullV.errdata{j});
        if ~silent
            fprintf('%s\n',FullV.errmsg{j});
        end
        newerrs = newerrs+1;
    end
end

if newerrs
    save(ename,'-struct','X');
end