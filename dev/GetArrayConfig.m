function ArrayConfig = GetArrayConfig(name,varargin)
% ArrayConfig = GetArrayConfig(name,varargin)
%name can be a directory containing .ns5 files, a .ns5 file name, or a
%struct returned from openNEV
%if a config is already in the directory that is loaded. 
%use GetArrayConfig(dirname, 'rebuild') to force a rebuild
%
%GetArrayConfig(dirname,'markbad',p)  
%  records that probe p is bad. 

rebuild = 0;
markprobe = 0;
setprobe=0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'rebuild',6)
        rebuild =1;
    elseif strncmpi(varargin{j},'mark',4)
        markprobe = find(strncmp(varargin{j},{'markbad', 'markdup'},7));
        j = j+1;
        setprobe = varargin{j};
    end
    j = j+1;
end

BlackRockPath();

aname = [];
badprobes = [];
if isfield(name,'MetaTags')
    nsx = name;
elseif isdir(name)
    aname = [name '/ArrayConfig.mat'];
    if exist(aname,'file') && rebuild == 0
        load(aname);
        ArrayConfig.type = SetArrayType(ArrayConfig);
        if setprobe && markprobe
            ArrayConfig.badprobes(setprobe) = markprobe;
            save(aname,'ArrayConfig');            
        end
        if ~isfield(ArrayConfig,'badprobes')
            ArrayConfig.badprobes = [];
        end
        if ~isfield(ArrayConfig,'id')
            rebuild = 1;
            badprobes = ArrayConfig.badprobes;
        else
            return;
        end
    end
    d = dir([name '/*.ns5']);
    if isempty(d)
        ArrayConfig = GetArrayByName(name);
        if isempty(ArrayConfig)
            fprintf('No .ns5 files in %s\n',name);
        else %used to be elseif rebuild. But getting here means need to save. 
            ArrayConfig.badprobes = badprobes; %keep list if rebuild
            save(aname,'ArrayConfig');
        end
        return;
    end
    nsx =openNSx([name '/' d(1).name]);
elseif ischar(name)
    [a,b,c] = fileparts(name);
    if strcmp(c,'.ns5')
    nsx = openNSx(name);
    elseif strcmp(c,'.mat')
        ArrayConfig = GetArrayConfig(a,varargin{:});
        return;
    end
end

ArrayConfig = ReadArrayConfig(nsx);
if ~isempty(aname) && ~isempty(ArrayConfig)
    save(aname,'ArrayConfig');
end

function type = SetArrayType(Array)
reset = 0;
type = 'unknown';

if isfield(Array,'type')
    type = Array.type;
elseif length(Array.X) == 24 && length(unique(Array.X) == 2)
    type = '12x2';
elseif length(Array.X) == 96
    type = 'Utah1';
else
    type = '24';
end

function Array = ReadArrayConfig(nsx)
Array = [];

if ~isfield(nsx.MetaTags,'ElecLabel') & isfield(nsx,'ElectrodesInfo') %V2.3
    for j = 1:length(nsx.ElectrodesInfo); 
        E(j) = sscanf(nsx.ElectrodesInfo(j).Label,'elec%d');
    end
else
if isempty(nsx.MetaTags.ElecLabel)
    return;
end
for j = 1:size(nsx.MetaTags.ElecLabel,1); 
    E(j) = sscanf(nsx.MetaTags.ElecLabel(j,:),'elec%d'); 
end
end

Array.Y = 1+mod(E-1,10);
Array.X = ceil(E/10);
Array.id = E;


function Array = GetArrayByName(name)

Array = [];

d = mydir('/bgc/group/arrays/*.mat');
for j = 1:length(d)
    load(d(j).name);
    Arrays{j} = ArrayConfig;
    labels{j} = ArrayConfig.label;
end

if isdir(name)
    d = dir([name '/*idx.mat']); %These at least need to be built
    for j = 1:length(d)
        load([name '/' d(j).name]);
        if exist('Expt','var') && isfield(Expt,'Comments') && isfield(Expt.Comments,'Peninfo')
           arrayname = Expt.Comments.Peninfo.trode; 
           for k = 1:length(labels)
               match(j,k) = sum(strfind(arrayname,labels{k}));
           end
           fprintf('%s:%s\n',d(j).name,arrayname);
        end
    end
    [a,b] = find(match);
    if ~isempty(b)
        [x,d] = Counts(b,'descend');
        Array = Arrays{d(1)};
    end
    
end


