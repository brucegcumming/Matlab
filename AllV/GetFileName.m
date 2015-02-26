function name = GetFileName(fpath)
%name = GetFileName(fpath) Retuns name part of a string (after directory)
%For stuctures, finds the name they were loaded from where possible
%See also GetFilePath
if ischar(fpath)
    [a,b,c] = fileparts(fpath);
    name = [b c];
elseif isstruct(fpath)
    E = fpath;
    if isfield(E,'loadname')
        name = E.loadname;
    elseif isfield(E,'Name')
        name = E.Name;
    elseif isfield(E,'name')
        name = E.name;
    elseif isfield(E,'Header')
        name = GetFileName(E.Header);
    end
end

