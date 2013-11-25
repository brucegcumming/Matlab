function ArrayConfig = GetArrayConfig(name,varargin)
% ArrayConfig = GetArrayConfig(name,varargin)
%name can be a directory containing .ns5 files, a .ns5 file name, or a
%struct returned from openNEV
%if a config is already in the directory that is loaded. 
%use GetArrayConfig(dirname, 'rebuild') to force a rebuild
%
%GetArrayConfig(dirname,'markbad',p)  cd ../../matlab
%  records that probe p is bad. 

rebuild = 0;
markprobe = 0;
markbadexpt = 0;
setprobe=0;
setlabel = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'badexpt',6)
        j = j+1;
        markbadexpt = varargin{j};
    elseif strncmpi(varargin{j},'rebuild',6)
        rebuild =1;
    elseif strncmpi(varargin{j},'set',3)
        j = j+1;
        setlabel = varargin{j};
    elseif strncmpi(varargin{j},'useufl',6)
        j = j+1;
        ArrayConfig = GetArrayFromUfl(name);
        return;
    elseif strncmpi(varargin{j},'mark',4)
        markprobe = find(strncmp(varargin{j},{'markbad', 'markdup'},7));
        j = j+1;
        setprobe = varargin{j};
    end
    j = j+1;
end



BlackRockPath();
ArrayConfig = [];
aname = [];
badprobes = [];
badexpts = [];
if isfield(name,'MetaTags')
    nsx = name;
elseif iscellstr(name) %Cell arary of strings - do all
    for j = 1:length(name)
        ArrayConfig{j} = GetArrayConfig(name{j});
    end
    return;
elseif iscell(name) && iscluster(name)
    if isfield(name{end},'loadname')
        ArrayConfig = GetArrayConfig(name{end}.loadname, varargin{:});
        return;
    elseif isfield(name{end},'spkfile')
    end
elseif isdir(name)
    aname = [name '/ArrayConfig.mat'];
    if ~isempty(setlabel)
        if exist(aname,'file')
            yn = questdlg(sprintf('Overwrite Array %s?',aname),'Overwrite','Yes','No','Yes');
        else
            yn = 'yes';
        end
        if strcmp(yn,'yes')
            ArrayConfig = GetArrayByName(setlabel);
            if isfield(ArrayConfig,'label')
                fprintf('Setting Array Config to %s\n',ArrayConfig.label);
                save(aname,'ArrayConfig');
            end
        else
            load(aname);
        end
        return;
    end
    if exist(aname,'file') && rebuild == 0
        load(aname);
        ArrayConfig.type = SetArrayType(ArrayConfig);
        if setprobe && markprobe
            ArrayConfig.badprobes(setprobe) = markprobe;
            BackupFile(aname);
            save(aname,'ArrayConfig');            
        elseif markbadexpt
            ArrayConfig.badexpts(markbadexpt) = 1;
            BackupFile(aname);
            save(aname,'ArrayConfig');                        
        end
        if ~isfield(ArrayConfig,'badprobes')
            ArrayConfig.badprobes = [];
        end
        if ~isfield(ArrayConfig,'badexpts')
            ArrayConfig.badexpts = [];
        end
        if ~isfield(ArrayConfig,'id')
            rebuild = 1;
            badprobes = ArrayConfig.badprobes;
            badexpts = ArrayConfig.badexpts;
        else
            if length(ArrayConfig.X) > 96
                if ~isfield(ArrayConfig,'CreationDate') || ArrayConfig.CreationDate < datenum('1/1/2014')
                    cprintf('red','%s has %d probes\n',aname,length(ArrayConfig.X));
                end
                    
            end
            return;
        end
    end
    d = dir([name '/*.ns5']);
    if isempty(d)
        ArrayConfig = GetArrayByName(name);
        if isempty(ArrayConfig)
            fprintf('No .ns5 files in %s\n',name);
            if exist(aname,'file')
                load(aname);
                if length(ArrayConfig.X) > 96
                    d = mydir([name '/*idx.mat']);
                    fix = 0;
                    if ~isempty(d)
                        load(d(1).name);
                        if Expt.Header.CreationDate < datenum('01/01/2014')
                            fix = 1;
                        end
                    else
                        fix = 1;
                    end
                    if fix
                        ArrayConfig.X = ArrayConfig.X(1:96);
                        ArrayConfig.Y = ArrayConfig.Y(1:96);
                        ArrayConfig.id = ArrayConfig.id(1:96);
                        ArrayConfig.fixed = 1;
                        save(aname,'ArrayConfig');
                    end
                end
            end
        else %used to be elseif rebuild. But getting here means need to save. 
            ArrayConfig.badprobes = badprobes; %keep list if rebuild
            ArrayConfig.badexpts = badexpts; %keep list if rebuild
            save(aname,'ArrayConfig');
        end
        return;
    end
    nsx =openNSx([name '/' d(1).name]);
elseif ~exist(name)
    cprintf('red','No file or directory %s\n',name);
    return;
elseif ischar(name)
    [a,b,c] = fileparts(name);
    if strcmp(c,'.ns5')
        nsx = openNSx(name);
    elseif strcmp(c,'.mat')
        ArrayConfig = GetArrayConfig(a,varargin{:});
        return;
    else
        cprintf('red','No Array config data in %s\n',name);
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
elseif ~isfield(Array,'X')
    type = 'unknown';
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
    np = size(nsx.MetaTags.ElecLabel,1);
    x = nsx.MetaTags.CreateDateTime;
    d = sprintf('%d/%d/%d',x(4),x(2),x(1));
    Array.CreationDate = datenum(d);    
    if np > 96 && nsx.MetaTags.CreateDateTime(1) < 2014
        np = 96;
    end
    for j = 1:np;
        if strncmp(nsx.MetaTags.ElecLabel(j,:),'chan',4)
            E(j) = sscanf(nsx.MetaTags.ElecLabel(j,:),'chan%d');
        else
            E(j) = sscanf(nsx.MetaTags.ElecLabel(j,:),'elec%d');
        end
    end
end

Array.Y = 1+mod(E-1,10);
Array.X = ceil(E/10);
Array.id = E;
ArrayConfig.src = 'ns5';

function Array = GetArrayFromUfl(name)
Array = [];
    d = mydir([name '/*.ufl']); %These at least need to be built
    if isempty(d)
        return;
    end
    etext = [];
    allpen = [];
    for j = 1:length(d)
        pen = [];
        txt = scanlines(d(j).name);
        for k = 1:length(txt)
            id = regexp(txt{k},' pe[0-9]');
            if ~isempty(id)
                a = sscanf(txt{k}(id(1)+3:end),'%d');
                pen(k) = a(1);
            end
            if strncmp('Electrode',txt{k},7)
                etext = txt{k};
            end
        end
        allpen = [allpen unique(pen(pen>0))];
    end
    if isempty(allpen)
        cprintf('red','No pen in UFL files 9n %s\n',name);
    end
    penid = mode(allpen);
    allpen = unique(allpen);
    pentxt = sprintf('/bgc/bgc/anal/%s/pens/pen%d.log', GetMonkeyName(name),penid);
    txt = scanlines(pentxt);
    id = find(strncmp('Electr',txt,6));
    if ~isempty(id)
        fprintf('Array Based on "%s" in penlog\n',txt{id(end)});
    elseif ~isempty(etext)
        txt{end+1} = etext;
        id = length(txt);
        fprintf('Array Based on "%s" in ufl\n',txt{id(end)});
    else
        return;
    end
    if strfind(txt{id(end)},'24Contact')
        Array = GetArrayByName('24probe');
        Array.idstr = txt{id(end)};
    elseif strfind(txt{id(end)},'8Contact')
        Array = GetArrayByName('8probe');
        Array.idstr = txt{id(end)};
    end
    if strfind(txt{id(end)},' 100u')
        Array.spacing = 100;
    elseif strfind(txt{id(end)},' 150u')
        Array.spacing = 150;
    elseif strfind(txt{id(end)},' 60u')
        Array.spacing = 60;
    elseif strfind(txt{id(end)},' 75u')
        Array.spacing = 75;
    elseif strfind(txt{id(end)},' 50u')
        Array.spacing = 50;
    end
    Array.src = 'penlog';

function Array = GetArrayByName(name)

Array = [];
labels = {};
match = [];
Arrays = {};
d = mydir('/bgc/group/arrays/*.mat');
for j = 1:length(d)
    try
        load(d(j).name);
        Arrays{j} = ArrayConfig;
        labels{j} = ArrayConfig.label;
        if strncmp(labels{j},'24 probe',6)
            aid(24) = j;
        elseif strncmp(labels{j},'8 probe',6)
            aid(8) = j;
        elseif strncmp(labels{j},'Utah96',6)
            aid(96) = j;
        end 
    catch
        mycprintf('red','Error Reading %s\n',d(j).name);
    end
end
if isempty(Arrays)
    return;
end
if isdir(name)
    d = dir([name '/*idx.mat']); %These at least need to be built
    if isempty(d)
        Array = GetArrayFromUfl(name);
        if isempty(Array)
        cprintf('blue','No idx.mat or .ufl files in %s\n',name);
        end
        return;
    end
    match = [];
    id = find(~(strcmp('FileIdx.mat',{d.name})));
    d = d(id);
    for j = 1:length(d)
        load([name '/' d(j).name]);
        if exist('Expt','var') && isfield(Expt,'Comments') && isfield(Expt.Comments,'Peninfo')
           arrayname = Expt.Comments.Peninfo.trode; 
           for k = 1:length(labels)
               match(j,k) = sum(strfind(arrayname,labels{k}));
           end
           if sum(match(j,:))              
               fprintf('%s:%s\n',d(j).name,arrayname);
           end
        elseif isfield(Expt,'Probes');
            p = [Expt.Probes.probe];
            if max(p) > 90
                Array = Arrays{aid(96)};
            elseif max(p) > 20
                Array = Arrays{aid(24)};
            elseif max(p) > 7
                Array = Arrays{aid(8)};
            end
            Array.src = 'guess'
            Array.idstr = num2str(max(p));
        end
    end
    [a,b] = find(match);
    if ~isempty(b)
        [x,d] = Counts(b,'descend');
        Array = Arrays{d(1)};
        Array.src = 'peninfo';
        Array.idstr = arrayname;
    elseif isempty(Array) || strcmp(Array.src,'guess')
        X = GetArrayFromUfl(name);
        if ~isempty(X);
            Array = X;
        end
    end
elseif ischar(name)
    for k = 1:length(labels)
        match(k) = sum(strfind(regexprep(labels{k},'\s+',''),regexprep(name,'\s+','')));
    end
    id = find(match > 0);
    if isempty(id)
        warndlg(sprintf('No Arrays Match name: %s',name));
    else
        Array = Arrays{id(1)};
        Array.src = 'label';
    end
end


