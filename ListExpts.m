function [names, details] = ListExpts(E, varargin)
%[names, details] = ListExpts(E, ...)   lists names of Expts in  E
% E can be a cell array of expts, or a directory name
% If its  directory name, only that directorly is listed
% ListExpts(E, 'depth',1) lists all directoris one layer down 
% For Directories, ListExpts Lookd in .idx files for experiments performed.
% ListExpts can search its own output to match specific expt names, and 
% then to search for combined Expt files with a given suffix
%
%ListExpts(...,'show','xx')
%      shows the value for field 'xx' in each expt. If this varies, the
%      mean value is printed
%
% [names, details] = ListExpts(dir);
% list = ListExpts(name, 'rds.dxXce', details);
% list = ListExpts(details, 'rds.dxXce');
% exlist = ListExpts(list,'AC');
%  
%         returns a structure with a list of directoies that contain this
%         expt

depth = 0;
psych = 0;
showvals = {};
j = 1;
argson = {};
while j <= length(varargin)
    if strncmp(varargin{j},'depth',5)
        j = j+1;
        depth = varargin{j};
    elseif strncmp(varargin{j},'psych',5)
        psych = 1;
    elseif strncmp(varargin{j},'show',5)
        j = j+1;
        showvals{end+1} = varargin{j};
    end
    j = j + 1;
end

if iscellstr(E) %list of names
    if isdir(E{1})
        for j = 1:length(E)
            [names{j}, details{j}] = ListExpts(E{j}, varargin{:});
        end
        return;
    else
        [names, details] = FindExpts(E, varargin{:});
    end
elseif psych
    [names,details] = ListPsychDir(E, 0);
elseif iscell(E)
    if iscellstr(E{1})
        [names, details] = FindExpts(E, varargin{:});
    elseif isfield(E{1},'dirpath')
        [names, details] = FindExpts(E, varargin{:});
    else
        [names, details] = ListExptCells(E, varargin{:});
    end
elseif isstruct(E)  % a details struct with a list of dirs containing expts. 
    names = CheckForExpts(E, varargin{:}); %first arg should be suffix to check
elseif isdir(E)
    if depth > 0
    [names, details] = ListExptDirs(E, depth, varargin);
    else
    [names, details] = ListExptDir(E, varargin);
    end
elseif ischar(E)
    [names, details] = ListExptDirs(E, depth, varargin);
end

function E = CheckForExpts(E, suffix, varargin)
 for j = 1:length(E.dirpath)
     s = [E.dirpath{j} '/*' suffix '*'];
     d = dir(s);
     E.nexpts(j) = length(d);
     if E.nexpts(j) == 0
         fprintf('%s no %s Expt Files\n',E.dirpath{j},suffix);
     end
     cells = [];
     for k = 1:length(d)
         id  =findstr(d(k).name,'cell');
         if ~isempty(id)
             cid = sscanf(d(k).name(id(1)+4:end),'%d');
             cells = [cells cid];
         end
     end
     E.cells{j} = cells;
     E.ncells(j) = length(cells);
 end


function [names, ids]= FindExpts(namelist, varargin)
names = {};
details = [];

if isempty(namelist)
    ids = [];
    return;
end

j = 1;
while j <= length(varargin)
    if iscell(varargin{j}) && isfield(varargin{j}{1},'dirpath')
        exdetails = varargin{j};
    end
    j = j+1;
end

if isstruct(namelist{1}) && isfield(namelist{1},'dirpath') %just details
    exdetails = namelist;
    namelist = {};
    for j = 1:length(exdetails)
        if isfield(exdetails{j},'names') %actually has data
            namelist{j} = exdetails{j}.names;
        end
    end
end

if iscellstr(namelist) && isdir(namelist{1}) %cell array of dir names
    for j = 1:length(namelist)
        [names{j}, ids{j}] = ListExptDir(namelist{j},varargin{:});
    end
    return;

elseif iscell(namelist{1}) %not a cell array of strings = array of cellstrings
    ids = {};
    for j = 1:length(namelist)
        if ~isempty(namelist{j})
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
    namelist = namelist(find(good)); %need to remove empties
    exdetails = exdetails(find(good));
    goodlist = find(good);
    for j = 1:length(namelist)
        [names{j}, id{j}] = FindExpts(namelist{j},varargin{:});
        details.matches(j) = length(names{j});
        if j <= length(exdetails) && isfield(exdetails{j},'ntrials')
        details.ntrials(j) = sum(exdetails{j}.ntrials(id{j}));
        end
        details.matchnames{j} = unique(names{j});
    end
    details.ids = id;
    clear id;
    id = find(details.matches);
    details.ids = details.ids(id);
    details.matches = details.matches(id);
    details.ntrials = details.ntrials(id);
    details.matchnames = details.matchnames(id);
    details.matchid = goodlist(id);
    allnames = {};
    for j = 1:length(details.matchnames)
        allnames = {allnames{:} details.matchnames{j}{:}};
    end
    [details.namecounts,details.allnames] = Counts(allnames);
    if ~isempty(details)
        for j = 1:length(id)
            details.dirpath{j} = exdetails{id(j)}.dirpath;
        end
    end
    ids = names;
    names = details;
    return;
end


%get here if its a cell array of strings that are not dir names
%= should be a list of expt names

findstr = {};
ids = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'test',4)
    elseif ischar(varargin{j})
        findstr = {findstr{:} varargin{j}};
    end
    j = j+1;
end
nf = 0;
for k = 1:length(findstr)
for j = 1:length(namelist)
    if strfind(namelist{j},findstr{k})
        nf = nf+1;
        names{nf} = namelist{j};
        ids(nf) = j;
    end
end
end


function [a,b] = CountReps(E)
et = E.Stimvals.et;
e2 = E.Stimvals.e2;
e3 = E.Stimvals.e3;
E.Stimvals.e0 = 0;
E = FillTrials(E,'e0');
x = [E.Trials.(et)];
y = [E.Trials.(e2)];
z = [E.Trials.(e3)];

xv = unique(x);
yv = unique(y);
zv = unique(z);
b = length(xv).*length(yv).*length(zv);
a = length(E.Trials)./b;


function [names, details] = ListExptCells(E, varargin)
showvals = [];
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'show',5)
        j = j+1;
        showvals{end+1} = varargin{j};
    end
    j = j+1;
end

for j = 1:length(E)
    if ~isfield(E{j}.Header,'expname')
        E{j}.Header.expname = Expt2Name(E{j});
    end
    names{j} = E{j}.Header.expname;
    if ~isfield(E{j}.Stimvals,'ei')
        E{j}.Stimvals.ei = 0;
    end
    if isfield(E{j},'Trials')
        details.ntrials(j) = length(E{j}.Trials);
    elseif isfield(E{j},'ntrials')
        details.ntrials(j) = E{j}.ntrials;
    end
    fn = fieldnames(E{j}.Trials);
    fn = setdiff(fn,{'Start' 'TrialStart' 'End' 'dur' 'id' 'TrueEnd' 'bstimes' 'delay' 'Trial' 'exvals' 'endevent' 'Spikes' 'count' 'OptionCode' 'rwtimes'});
    fn = setdiff(fn,{'ch' 'op' 'rw' 'se' });
    str = [];
    idstr = [];
    if isfield(E{j}.Trials,'id')
        idstr = sprintf(' Id%d-%d',minmax([E{j}.Trials.id]));
    end
    
    for k = 1:length(fn)
        str = [str fn{k} ','];
    end
    xstr = ExptString(E{j},showvals);
    for k = 1:length(showvals)
        details.(showvals{k})(j) = GetEval(E{j},showvals{k});
    end
    if isfield(E{j}.Header,'human')
        [a,b] = CountReps(E{j});
        fprintf('%d: (%d * %.1frep = %dTrials) %s\n',j,b,a,details.ntrials(j),Expt2Name(E{j}));
    else
        fprintf('%d: (%dTrials%s) %s %s %s. Fields %s\n',j,details.ntrials(j),idstr,names{j},xstr,Expt2Name(E{j}),str);
    end
end


function [names, details] = ListExptDir(E, varargin)


relist = 0;
listfullv = 0;
j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'fullv',5)
        listfullv = 1;
    elseif strncmpi(varargin{j},'relist',5)
        relist = 1;
    end
    j = j+1;
end

details.dirpath = E;
if listfullv
    expts = [];
    probes = [];
    sizes = [];
    d = dir([E '/Expt*FullV.mat']);
    for j = 1:length(d)
        expts(j) = GetExptNumber(d(j).name);
        probes(j) = GetProbeFromName(d(j).name);
        sizes(j) = d(j).bytes;
    end
    names = unique(expts);
    details.expts = expts;
    details.probes = probes;
    details.sizes = sizes;
    return;
end

outname = [E '/AllExptList.mat'];
if exist(outname) && relist == 0
    load(outname);
    names = details.names;
    if ~isfield(details,'type') && isfield(details,'filename')
        nf = length(unique(details.filename));
        if nf == 1
            details.type = 'WaveMark';
        elseif nf > length(details.filename)/2
            details.type = 'Continuous';
        else
            details.type = 'Unknown';
        end
    end
    if ~isfield(details,'dirpath')
        details.dirpath = E;
    end
    return;
end

d = dir([E '/*idx.mat']);
nx = 0;
names = {};
details = [];
details.dirpath = E;
for j = 1:length(d)
    if sum(strcmp(d(j).name,{'FileIdx.mat' 'OldIdx.mat'}))
    else
    name = [E '/' d(j).name];
    clear ExptList;
    load(name);
    id = regexp(d(j).name,'\.[0-9]*idx');
    if isempty(id)
        fileno = 0;
    else
        fileno = sscanf(d(j).name(id(1)+1:end),'%d');
    end
    if ~exist('ExptList','var') && exist('Expt','var')
        ExptList = BuildExptList(Expt, Expts);
    end
    if exist('ExptList','var')
        for k = 1:length(ExptList)
            e = find([Expts.start] == ExptList(k).start);
            if isempty(Expts(e).result)
                Expts(e).result = 2;
            end
            if length(e) && Expts(e).result == 2
                nx = nx+1;
                names{nx} = ExptList(k).expname;
                details.ntrials(nx) = Expts(e).lasttrial-Expts(e).firsttrial;
                details.filename{nx} = d(j).name;
                details.exptno(nx) = k+fileno;
            end
        end
    end
    end
end
if nx > 0
    [a,b] = sort(details.exptno);
    names = names(b);
    if sum(fileno) == 0
        details.type = 'WaveMark';
    elseif sum(fileno > 0) > length(a)/2
        details.type = 'Continuous';
    else
        details.type = 'Unknown';
    end
    details.ntrials = details.ntrials(b);
    details.filename = details.filename(b);
    details.exptno = details.exptno(b);
    details.names = names;
    try
        save(outname,'details');
    end
elseif ~exist(outname) || relist
    details.names = {};
    details.exptno = [];
    details.type = 'Empty';
    save(outname,'details');
end


function ExptList = BuildExptList(Expt, Expts)

ExptList = [];

nx = 0;
for j = 1:length(Expts)
    if Expts(j).lasttrial - Expts(j).firsttrial > 10
        nx = nx+1;
        ts = Expts(j).firsttrial;
        if isfield(Expt.Trials,'et')
            ExptList(nx).et = Expt.Trials.et{ts};
            ExptList(nx).e2 = Expt.Trials.e2{ts};
            ExptList(nx).e3 = Expt.Trials.e3{ts};
        elseif isfield(Expts,'et')
            ExptList(nx).et = Expts(j).et;
            ExptList(nx).e2 = Expts(j).e2;
            ExptList(nx).e3 = Expts(j).e3;
        end
        ExptList(nx).start = Expts(j).start;
        ExptList(nx).expname = [ExptList(nx).et 'X' ExptList(nx).e2 'X' ExptList(nx).e3];
        if isfield(Expt.Trials,'OptionCode')
        if strfind(Expt.Trials.OptionCode{ts},'+fS')
            ExptList(nx).expname = [ExptList(nx).expname 'FS'];
        end
        end
    end
end

function [names, details] = ListPsychDir(E, depth, varargin)
names = {};
details = {};
d = dir(E);
if isdir(E)
    root = [E '/'];
else
    id = strfind(E,'/');
    if length(id)
        root = E(1:id(end));
    end
end

nd = 0;
ngood = 0;

for j = 1:length(d)
     path = [root d(j).name];
     [a,b,c] = fileparts(d(j).name);
     if ~d(j).isdir && isempty(c)
     expts = PsychMon(path,'getexpts');
     if ~isempty(expts)
     [names{j}, details{j}] = ListExpts(expts);
     end
     end
end


function [names, details] = ListExptDirs(E, depth, varargin)

names = {};
details = {};
d = dir(E);
if isdir(E)
    root = [E '/'];
else
    id = strfind(E,'/');
    if length(id)
        root = E(1:id(end));
    end
end

nd = 0;
ngood = 0;

for j = 1:length(d)
    if d(j).isdir && d(j).name(1) ~= '.'
        if depth > 1
            path = [root d(j).name];
            [a, b] = ListExptDirs(path, depth-1, varargin);
            for k = 1:length(a)
                nd = nd+1;
                names{nd} = a{k};
                details{nd} = b{k};
                ngood = ngood+1;
            end
        else
            nd = nd+1;
            path = [root d(j).name];
            [names{nd}, details{nd}] = ListExptDir(path, varargin);
            details{nd}.dirpath = path;
            if ~isempty(names{nd})
                ngood = ngood+1;
            end
        end
    end
end

if ngood == 0
    names = {};
    details = {};
end