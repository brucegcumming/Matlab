function date = CreationDate(name, varargin)
%Find the creation date of a file (E.g. smr file)
%from its contents, not from the directory

date = NaN;
if exist(name,'dir')
elseif ~isempty(strfind(name,'.mat')) && exist(name)
    F = load(name);
    if isfield(F,'Ch30'); %.mat made from smr
        id = strmatch('Created',F.Ch30.text);
        if isempty(id)
            id = strmatch('uf',F.Ch30.text);
        end
        if ~isempty(id)
            offset = 8;
            s = strfind(F.Ch30.text(id(1),:),'Created');
            if isempty(s)
                s = strfind(F.Ch30.text(id(1),:),'smr ');
                offset = 4;
            end
            if ~isempty(s)
            date = datenum(F.Ch30.text(id(1),s(1)+offset:end));
            end
        end
    end
end
    