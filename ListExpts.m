function [names, details] = ListExpts(E, varargin)
%[names, details] = ListExpts(E, ...)   lists names of Expts in  E
% E can be a cell array of expts, or a directory name
% If its  directory name, only that directorly is listed
% ListExpts(E, 'depth',1) lists all directoris one layer down 
% For Directories, ListExpts Lookd in .idx files for experiments performed.
% ListExpts can search its own output to match specific expt names, and 
% then to search for combined Expt files with a given suffix
%
% [names, details] = ListExpts(dir);
% list = ListExpts(name, 'rds.dxXce', details);
% exlist = ListExpts(list,'AC');
%  
%         returns a structure with a list of directoies that contain this
%         expt

depth = 0;
j = 1;
argson = {};
while j <= length(varargin)
    if strncmp(varargin{j},'depth',5)
        j = j+1;
        depth = varargin{j};
    end
    j = j + 1;
end

if iscellstr(E) %list of names
    FindExpts(E, varargin);
elseif iscell(E)
    if iscellstr(E{1})
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
if iscell(namelist{1})
    ids = {};
    for j = 1:length(namelist)
        if ~isempty(namelist{j})
            good(j) = 1;
        else
            good(j) = 0;
        end
    end
   % namelist = namelist(find(good));
    for j = 1:length(namelist)
        [names{j}, id{j}] = FindExpts(namelist{j},varargin{:});
        details.matches(j) = length(names{j});
        if j <= length(exdetails) && isfield(exdetails{j},'ntrials')
        details.ntrials(j) = sum(exdetails{j}.ntrials(id{j}));
        end
    end
    details.ids = id;
    clear id;
    id = find(details.matches);
    details.ids = details.ids(id);
    details.matches = details.matches(id);
    details.ntrials = details.ntrials(id);
    details.matchid = id;
    if ~isempty(details)
        for j = 1:length(id)
            details.dirpath{j} = exdetails{id(j)}.dirpath;
        end
    end
    ids = names;
    names = details;
    return;
end


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

for j = 1:length(E)
    names{j} = E{j}.Header.expname;
    if isfield(E{j},'Trials')
        details.ntrials(j) = length(E{j}.Trials);
    elseif isfield(E{j},'ntrials')
        details.ntrials(j) = E{j}.ntrials;
    end
    if isfield(E{j}.Header,'human')
        [a,b] = CountReps(E{j});
        fprintf('%d: (%d * %.1frep = %dTrials) %s\n',j,b,a,details.ntrials(j),Expt2Name(E{j}));
    else
        fprintf('%d: (%dTrials) %s\n',j,details.ntrials(j),names{j});
    end
end


function [names, details] = ListExptDir(E, varargin)

relist = 0;
j = 1;
while j <=length(varargin)
    if strncmpi(varargin{j},'relist',5)
        relist = 1;
    end
    j = j+1;
end

outname = [E '/AllExptList.mat'];
if exist(outname) && relist == 0
    load(outname);
    names = details.names;
    return;
end

d = dir([E '/*idx.mat']);
nx = 0;
names = {};
details = [];
for j = 1:length(d)
    if strcmp(d(j).name,'FileIdx.mat')
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
    if ~exist('ExptList','var')
        ExptList = BuildExptList(Expt, Expts);
    end
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
if nx > 1
    [a,b] = sort(details.exptno);
    names = names(b);
    details.ntrials = details.ntrials(b);
    details.filename = details.filename(b);
    details.exptno = details.exptno(b);
    details.names = names;
    try
        save(outname,'details');
    end
elseif ~exist(outname)
    details.names = {};
    details.exptno = [];
    save(outname,'details');
end


function ExptList = BuildExptList(Expt, Expts)

ExptList = [];

nx = 0;
for j = 1:length(Expts)
    if Expts(j).lasttrial - Expts(j).firsttrial > 10
        nx = nx+1;
        ts = Expts(j).firsttrial;
        ExptList(nx).et = Expt.Trials.et{ts};
        ExptList(nx).e2 = Expt.Trials.e2{ts};
        ExptList(nx).e3 = Expt.Trials.e3{ts};
        ExptList(nx).start = Expts(j).start;
        ExptList(nx).expname = [ExptList(nx).et 'X' ExptList(nx).e2 'X' ExptList(nx).e3];
        if strfind(Expt.Trials.OptionCode{ts},'+fS')
            ExptList(nx).expname = [ExptList(nx).expname 'FS'];
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