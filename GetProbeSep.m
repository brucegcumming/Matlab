function probesep = GetProbeSep(Header, varargin)

probesep = NaN;
if isfield(Header,'Peninfo')
id = strfind(Header.Peninfo.trode,'Contact');
if length(id)
    x = id(1);
    id = strfind(Header.Peninfo.trode(id:end),' ');
    probesep = sscanf(Header.Peninfo.trode(id+x:end),'%d');
end
end
