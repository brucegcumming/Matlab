function ArrayConfig = GetArrayConfig(name,varargin)
% ArrayConfig = GetArrayConfig(name,varargin)
%name can be a directory containing .ns5 files, a .ns5 file name, or a
%struct returned from openNEV
%if a config is already in the directory that is loaded. 
%use GetArrayConfig(dirname, 'rebuild') to force a rebuild
rebuild = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'rebuild',6)
        rebuild =1;
    end
    j = j+1;
end

BlackRockPath();

aname = [];
if isfield(name,'MetaTags')
    nsx = name;
elseif isdir(name)
    aname = [name '/ArrayConfig.mat'];
    if exist(aname,'file') && rebuild == 0
        load(aname);
        if ~isfield(ArrayConfig,'id')
            rebuild = 1;
        else
        return;
        end
    end
    d = dir([name '/*.ns5']);
    if isempty(d)
        ArrayConfig = [];
        fprintf('No .ns5 files in %s\n',name);
        return;
    end
    nsx =openNSx([name '/' d(1).name]);
elseif ischar(name)
    [a,b,c] = fileparts(name);
    if strcmp(c,'.ns5')
    nsx = openNSx(name);
    elseif strcmp(c,'.mat')
        ArrayConfig = GetArrayConfig(a);
        return;
    end
end

ArrayConfig = ReadArrayConfig(nsx.MetaTags);
if ~isempty(aname) && ~isempty(ArrayConfig)
    save(aname,'ArrayConfig');
end

function Array = ReadArrayConfig(MetaTags)
Array = [];
if isempty(MetaTags.ElecLabel)
    return;
end
np = size(MetaTags.ElecLabel,1);
if np > 96 && MetaTags.CreateDateTime(1) < 2014
    np = 96;
end
for j = 1:np
    E(j) = sscanf(MetaTags.ElecLabel(j,:),'elec%d'); 
end

Array.Y = 1+mod(E-1,10);
Array.X = ceil(E/10);
Array.id = E;