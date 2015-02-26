function name = ExptName(Expt, varargin)
shortname = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'long',4)
        shortname = 0;
    end
    j = j+1;
end

if isfield(Expt,'Header') 
    if isfield(Expt.Header,'Name')
        str = Expt.Header.Name;
    elseif isfield(Expt.Header,'name')
        str = Expt.Header.name;
    else
        str = 'Header Missing Name';
    end
    if shortname
        [dpath, name] = fileparts(str);
    else
        name = str;
    end
end