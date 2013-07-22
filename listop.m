function details = listop(list, op, varargin)
%
% listop(list, 'copy', ....)
% copies all files in list, on in the current drive to C:
% listop(list, 'copy', 'cpem',....)
%   copies over any .em.mat files that match files in list
%
%
% copies files
% Do useful things to files named in list.

cprefix = 'C:';
skip = 1;
doac = 1;
doorfiles = 1;
doemfiles = 0;
funcfcn = [];
funcargs = {};
if iscellstr(list)
    fstrings = list;
elseif ischar(list)
    fstrings = textread(list,'%s','delimiter','\n');
end
Expts = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sizestrings',7)
        OTTF.sizestrings = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'cpem',4)
        doemfiles = 1;
    elseif strncmpi(varargin{j},'function',3)
        j = j+1;
        funcfcn = varargin{j};
        if length(varargin) > j & iscell(varargin{j+1})
            j = j+1;
            funcargs = varargin{j};
        end
    elseif strncmpi(varargin{j},'size',2)
        j = j+1;
        wsize = varargin{j};
    elseif strncmpi(varargin{j},'skip',2)
        j = j+1;
        skip = varargin{j};
    end
    j = j+1;
end

details.ndone = 0;
ndone = 0;
names = [];

rid = strmatch('#',fstrings);
id = setdiff(1:length(fstrings),rid);
fstrings = fstrings(id);


for j = skip:length(fstrings)
if length(fstrings{j}) == 0;
elseif length(funcfcn)
    pathname = fstrings{j};
    if isempty(funcargs)
        good = eval([funcfcn '(''' pathname ''')']);
    else
        str = [];
        for k = 1:length(funcargs)
            str = [str ',''' funcargs{k} ''''];
        end
        good = eval([funcfcn '(''' pathname '''' str ')']);
    end
elseif strncmpi(op,'recombineorbw',12)
    name = deblank(regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat'));
    name = regexprep(name,'\.cell[0-9]*\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombineorbw','recount');
    end
elseif strncmpi(op,'allprobes',6)
    if regexp(fstrings{j},'\.cell[0-9]*\..*.mat')
    name = regexprep(fstrings{j},'\.cell[0-9]*\..*.mat','.mat');
    if ndone == 0  | isempty(strmatch(name,names))
        PlotAllProbes(fileparts(name),'save');
    ndone = ndone+1;
    names{ndone} = name;
    end
    end
    details.names = names;
    details.ndone = ndone;
elseif strncmpi(op,'load',4)
    ndone = ndone+1;
    Expts{ndone} = LoadExpt(fstrings{j});
elseif strncmpi(op,'recombineall',11)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombinequit','recount');
    end
elseif strncmpi(op,'recombine',6)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) % not a brainwave file
        combine(name,'recombinename',fstrings{j});
    ndone = ndone+1;
    end
elseif strncmpi(op,'chkcluster',6)
    f = load(fstrings{j});
    if isfield(f,'cExpt')
        f.cExpt.Header.Filename = fstrings{j};
        [dp(j), d] = ExptCellQuality(f.cExpt,'verbose');
        ndone = ndone+1;
    elseif isfield(f,'Expt')
        Expt.Header.Filename = fstrings{j};
        [dp(j), d] = ExptCellQuality(Expt,'verbose');
        ndone = ndone+1;
    end
    details.errs(j) = d.errs;
        
elseif strncmpi(op,'relist',6)
    name = regexprep(fstrings{j},'\.c[0-9]\..*.mat','.mat');
    name = regexprep(name,'\.cell[0-9]\..*.mat','.mat');
    if isempty(strfind(name,'.0.mat')) & (ndone == 0  | isempty(strmatch(name,names)))
        APlaySpkFile(name,'relist');
        ndone = ndone+1;
        names{ndone} = name;
    end
    details.names = names;
    details.ndone = ndone;
elseif strncmpi(op,'copy',4)
    tgt = [cprefix fstrings{j}];

    a = CopyFile(fstrings{j},tgt);
    ndone = ndone+a;
    details.ndone = ndone;
    if doemfiles
        emfile = regexprep(fstrings{j},'\..*.mat','.em.mat');
        if ~strcmp(emfile,fstrings{j}) & exist(emfile,'file');
        tgt = [cprefix emfile];
        CopyFile(emfile,tgt);
        end
    end
    if doorfiles
        otfile = regexprep(fstrings{j},'image.ORBW.','image.OT.');
        if ~strcmp(otfile,fstrings{j}) & exist(otfile,'file');
        tgt = [cprefix otfile];
        CopyFile(otfile,tgt);
        end
        otfile = regexprep(fstrings{j},'image.ORBW.','grating.OXM.');
        if ~strcmp(otfile,fstrings{j}) & exist(otfile,'file');
        tgt = [cprefix otfile];
        CopyFile(otfile,tgt);
        end
    end
    if doac
        id = regexp(fstrings{j},'[0-9,F][A-Z][A-Z]*\.mat');
        suffix = fstrings{j}(id:end-4);
        [a,b] = regexp(fstrings{j},'\.[a-z][a-z]*\.','start','end');
        stim = fstrings{j}(a+1:b-1);
        acfile = regexprep(fstrings{j},suffix,'AC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},suffix,'OXAC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.nsines.DP');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.rds.AC');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
        acfile = regexprep(fstrings{j},['\.' stim '\.' suffix],'.rds.DT');
        if ~strcmp(acfile,fstrings{j}) & exist(acfile,'file');
        tgt = [cprefix acfile];
        CopyFile(acfile,tgt);
        end
    end
end
end
if length(Expts)
    details.Expts = Expts;
end
details.ndone = ndone;

function done = CopyFile(src, tgt, varargin)
%copy src to tgt, if tgt doesn't exist or is older than src.
done = 0;
sd = dir(src);
if isempty(sd)
    fprintf('Cant Find %s\n',src);
    return;
end
    d = dir(tgt);
    dd = dir(fileparts(tgt));
    if length(d) == 0  || d(1).datenum < sd(1).datenum
        if length(dd) == 0
           mkdir(fileparts(tgt));
        end
        copyfile(src, tgt);
        cmd = ['copy ' src ' ' tgt];
        cmd = strrep(cmd,'/','\');
        fprintf('%s\n',cmd);
        done = 1;
    end

