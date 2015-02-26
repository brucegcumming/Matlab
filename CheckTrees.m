function [missing, details] = CheckTrees(src,tgt, varargin)
%[missing, details] = CheckTrees(src,tgt) check that files on src exist somewhere on tgt
fargs = {};
missing = {};
details = {};
parallel = 1;
dironly = 0;
hidelist = {'openMotif'};
if nargin == 1
    tgt = {};
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dironly',5)
        dironly = 1;
    elseif strncmpi(varargin{j},'hide',5)
        j = j+1;
        hidelist = {hidelist{:} varargin{j}};
    elseif strncmpi(varargin{j},'parallel',5)
        parallel = 1;
    elseif strncmpi(varargin{j},'serial',5)
        parallel = 0;
    end
    j = j+1;
end

if ischar(tgt) && strcmp(tgt,'serial')
    parallel = 0;
end
ts = now;
errs = {};
if iscellstr(src) && iscellstr(tgt)
    if parallel
        parfor j = 1:length(src)
            [newmiss{j}, details{j}] = CheckTrees(src{j},tgt{j});
        end
    else
        for j = 1:length(src)
            [newmiss{j}, details{j}] = CheckTrees(src{j},tgt{j});
        end
    end
    for j = 1:length(src)
        missing = {missing{:} newmiss{j}{:}};
    end
    return;
elseif iscellstr(src) && sum(strncmp(tgt,{'list' 'size'},4))
    CheckTreeResult(src, tgt, hidelist);
    return;
elseif sum(strcmp(src,{'oldservers' 'bgc2' 'bgc3' 'bgc' 'bgc4' 'bgc5' 'bgc6' 'bgc7'}))
    search = src;
    clear src;
    src = {};
    tgt = {};
    ns = 0;
    if sum(strcmp(search,{'oldservers' 'bgc3'}))
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/uffs/dufus';
    tgt{ns} = '/b/data/dufus';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/uffs/rufus';
    tgt{ns} = '/b/data/rufus';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/uffs/icarus';
    tgt{ns} = '/b/data/icarus';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/cicmiln';
    tgt{ns} = '/Volumes/bgcpriv/cicmiln';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/krugk';
    tgt{ns} = '/Volumes/bgcpriv/krugk';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/ralf';
    tgt{ns} = '/Volumes/bgcpriv/ralf';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/group';
    tgt{ns} = '/b/group';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/data/ruby';
    tgt{ns} = '/b/data/ruby';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/data/psych';
    tgt{ns} = '/b/data/psych';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/smr/dae';
    tgt{ns} = '/b/data/dae';
    ns = ns+1;
    src{ns} = '/Users/bgc/tmp/smr/icarus';
    tgt{ns} = '/b/data/icarus';
    end
    if sum(strcmp(search,{'oldservers' 'bgc4'}))
    [a, b, c] = GetSubDirs('/Volumes/bgc4/ox', '/Volumes/bgcdiv/files/ox','two');
    a = setdiff(a,{'/Volumes/bgc4/ox/usr/data'});
    b = setdiff(a,{'/Volumes/bgcdiv/files/ox/usr/data'});
    src = {src{:} a{:}};
    tgt = {tgt{:} b{:}};
    errs = {errs{:} c{:}};
    
%    [a, b, errs] = GetSubDirs('/Volumes/bgc4/machines', '/Volumes/bgcdiv/files/machines');
%    src = {src{:} a{:}};
%    tgt = {tgt{:} b{:}};
    ns = length(src)+1;
    src{ns} = '/Volumes/bgc4/uffs/hogan';
    tgt{ns} = '/b/data/hogan';
    end
    if sum(strcmp(search,{'oldservers' 'bgc2'}))
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/anal';
    tgt{ns} = '/b/bgc/anal';
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/data/hogan';
    tgt{ns} = '/b/data/hogan';
    [a, b, c] = GetSubDirs('/Volumes/bgc2/hn', '/Volumes/bgcpriv/hn');
    src = {src{:} a{:}};
    tgt = {tgt{:} b{:}};
    errs = {errs{:} c{:}};
    [a, b, c] = GetSubDirs('/Volumes/bgc2/smr/lem', '/b/data/lem');
    src = {src{:} a{:}};
    tgt = {tgt{:} b{:}};
    errs = {errs{:} c{:}};
    
    [a, b, c] = GetSubDirs('/Volumes/bgc2/smr/lem', '/b/data/lem/SE');
    src = {src{:} a{:}};
    tgt = {tgt{:} b{:}};
    errs = {errs{:} c{:}};

    ns = length(src)+1;
    src{ns} = '/Volumes/bgc2/smr/dae';
    tgt{ns} = '/b/data/dae';
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/smr/ica';
    tgt{ns} = '/b/data/ica';
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/smr/icarus';
    tgt{ns} = '/b/data/icarus';
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/smr/jbe';
    tgt{ns} = '/b/data/jbe';
    ns = ns+1;
    src{ns} = '/Volumes/bgc2/smr/rufus';
    tgt{ns} = '/b/data/rufus';
    end

    if sum(strcmp(search,{'oldservers' 'bgc'}))
        [a, b, c] = GetSubDirs('/Volumes/bgc/bgc', '/Volumes/bgcpriv/bgc');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};

        ns = length(src)+1;
        src{ns} = '/Volumes/bgc/group';
        tgt{ns} = '/b/group';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/aidan';
        tgt{ns} = '/Volumes/bgcpriv/aidan';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/ali';
        tgt{ns} = '/Volumes/bgcpriv/ali';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/bonyda';
        tgt{ns} = '/Volumes/bgcpriv/bondya';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/ceb';
        tgt{ns} = '/Volumes/bgcpriv/ceb';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/tanabes';
        tgt{ns} = '/Volumes/bgcpriv/tanabes';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/cicmiln';
        tgt{ns} = '/Volumes/bgcpriv/cicmiln';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/jcr';
        tgt{ns} = '/Volumes/bgcpriv/jcr';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/kangi';
        tgt{ns} = '/Volumes/bgcpriv/kangi';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/krugk';
        tgt{ns} = '/Volumes/bgcpriv/krugk';
        ns = ns+1;
        src{ns} = '/Volumes/uffs/ruby';
        tgt{ns} = '/b/data/ruby';
        ns = ns+1;
        src{ns} = '/Volumes/bgc/files';
        tgt{ns} = '/b/group';
    end

    if sum(strcmp(search,{'oldservers' 'bgc5'}))
        ns = ns+1;
        src{ns} = '/Volumes/bgc5/MRI';
        tgt{ns} = '/b/data/jbe/MRI';
        ns = ns+1;
        src{ns} = '/Volumes/bgc5/Psych Overflow';
        tgt{ns} = '/b/data/psych';
        
        
        [a, b, c] = GetSubDirs('/Volumes/bgc5/smr/lem', '/b/data/lem');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};
        [a, b, c] = GetSubDirs('/Volumes/bgc5/smr/jbe', '/b/data/jbe');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};
        [a, b, c] = GetSubDirs('/Volumes/bgc5/smr/rufus', '/b/data/rufus');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};

        
        ns = length(src)+1;;
        src{ns} = '/Volumes/bgc5/files';
        tgt{ns} = '/Volumes/bgcdiv/files';
    end
    if sum(strcmp(search,{'oldservers' 'bgc6'}))
        ns = ns+1;
        src{ns} = '/Volumes/bgc6/ORBWRC';
        tgt{ns} = '/b/data/psych/ORBWRC';
        ns = ns+1;

        [a, b, c] = GetSubDirs('/Volumes/bgc6/smr/lem', '/b/data/lem');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};
        [a, b, c] = GetSubDirs('/Volumes/bgc6/smr/jbe', '/b/data/jbe');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        [a, b, c] = GetSubDirs('/Volumes/bgc6/smr/ica', '/b/data/ica');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = {errs{:} c{:}};
        [a, b, c] = GetSubDirs('/Volumes/bgc6/smr/dae', '/b/data/dae');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};

    end
    if sum(strcmp(search,{'oldservers' 'bgc7'}))
        [a, b, c] = GetSubDirs('/Volumes/bgc7/smr/lem', '/b/data/lem');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
        errs = c;
        [a, b, c] = GetSubDirs('/Volumes/bgc7/Utah/jbe', '/b/data/jbe');
        src = {src{:} a{:}};
        tgt = {tgt{:} b{:}};
    end

    [missing, details] = CheckTrees(src, tgt);
    missing = {errs{:} missing{:}};
    return;
elseif iscell(tgt) && length(tgt) > 10 %a saved TreeFind
    tfiles = tgt;
elseif iscell(tgt) %build list from several folters;
else
    [tfiles, tsz, tdates, tex] = TreeFind(tgt,fargs{:});
end
[sfiles, ssz, sdates, sex] = TreeFind(src,fargs{:});
worker = mygetCurrentTask('num');
fprintf('Worker %d Listing %s took %.2f. %d and %d files At %s\n',worker,src,mytoc(ts),length(sfiles),length(tfiles),datestr(now));
if dironly
    details.n = [length(sfiles) length(tfiles)];
    return;
end
xtime = now;
nmatch = 0;

for j = 1:length(sfiles)
    id = strfind(sfiles{j},'/');
    c = '';
    nmatch = 0;
    if length(id) > 1 && ~strncmp('/Volumes/bgc/bgc/bgc',sfiles{j},20)
        name = sfiles{j}(id(end-1):end);
        found = 0;
        matches = [];
        nf = 0;
        if ~sum(strcmp(name,{'..' '.'}))
            for k = 1:length(tfiles)
                id = strfind(tfiles{k},name);
                if ~isempty(id)
                    if ssz(j) == tsz(k) %name and size match
                        found = found+1;
                    elseif strcmp(regexprep(sfiles{j},'.*/',''),regexprep(tfiles{k},'.*/',''))
                        nf = nf+1;
                        matches(nf) = k;
                    else %haven't found file yet
                    end
                end
            end
            t = now;
            if t - xtime > 1./(24 * 10); %1 min elapsed
                xtime = t;
                fprintf('Worker%d: %s: %d/%d files checked in %s\n',worker,datestr(now),j, length(sfiles),src);
            end
            if found == 0
                xtime = t;
                for m = 1:nf
                    %if size does not match, make sure name is really identical
                    if strcmp(strrep(sfiles{j},src,''),strrep(tfiles{matches(m)},tgt,''))
                        fprintf('Size mismatch %d for %s %d:',m,sfiles{j},ssz(j)); %identical path up to prefix
                    else
                        fprintf('Possible match %d for %s %d:',m,sfiles{j},ssz(j));
                    end
                    fprintf('%s %d',tfiles{matches(m)},tsz(matches(m)));
                    if tdates(matches(m)) > sdates(j)
                        fprintf('Newer');
                        c = '!';
                    else
                        c = '*';
                    end
                    fprintf('\n');
                    details.matches{length(missing)+1} = matches;
                end
                if isempty(strfind(sfiles{j},'Images/bpnoise')) %don't care about these
                    fprintf('No match for %s\n',sfiles{j});
                    missing{end+1} = [sfiles{j} c];
                end
            else
                nmatch = nmatch+1;
            end
        end
    end
end
details.nmatch = nmatch;
fprintf('Checking %s took %.2f. At %s. %d matches %d missing\n',src,mytoc(ts),datestr(now),nmatch,length(missing));


function [srcsub, tgtsub, errs] = GetSubDirs(src, tgt, varargin)

recurse = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'two',3)
        recurse = 1;        
    elseif strncmpi(varargin{j},'recurse',5)
        j = j+1;
        recurse = varargin{j};        
    end
    j = j+1;
end
a = dir(src);
b = dir(tgt);
a = a([a.isdir] > 0);
b = b([b.isdir] > 0);
srcsub = {};
tgtsub = {};
errs = {};
for j = 1:length(a)
    id = find(strcmp(a(j).name,{b.name}));
    if sum(strcmp(a(j).name,{'.' '..' 'bgc8'}))
    elseif length(id) ==1
        srcsub{end+1} = [src '/' a(j).name];
        tgtsub{end+1} = [tgt '/' b(id).name];
    elseif isempty(id)
        errs{end+1} = sprintf('Mo Match for subdir %s/%s in %s',src,a(j).name,tgt);
        if strncmp(a(j).name,'M1',2)
        fprintf('%s\n',errs{end});
        end
        fprintf('%s\n',errs{end});
    end
end

if recurse
    newsrc = {};
    newtgt = {};
    newerrs = {};
    for j = 1:length(srcsub)
        [a,b, c] = GetSubDirs(srcsub{j},tgtsub{j},'recurse',recurse-1);
        newsrc = {newsrc{:} a{:}};
        newtgt = {newtgt{:} b{:}};
        newerrs = {newerrs{:} c{:}};
    end
    srcsub = newsrc;
    tgtsub = newtgt;
end
    

function CheckTreeResult(missing, type,hidelist)


if strcmp(type,'list')
for j = 1:length(missing)
    if missing{j}(end) == '*' || missing{j}(end) == '!'; %size mismatch
    elseif strfind(missing{j},'ruby')
    elseif strfind(missing{j},'spkblk')
    elseif regexp(missing{j},'/[\.]*$')
    elseif strfind(missing{j},'/._')
    elseif sum(cellstrfind(missing{j},hidelist))
    else
        
        fprintf('%s\n',missing{j});
    end
end
elseif strncmp(type,'size',4)
    if strcmp(type,'sizeall')
        showall = 1;
    else
        showall = 0;
    end
    for j = 1:length(missing)
        if missing{j}(end) == '*' || missing{j}(end) == '!'
            filename = missing{j}(1:end-1);
            if strncmp(filename,'/Volumes/bgc5/smr/lem/Done',15)
                tgt = strrep(filename,'/Volumes/bgc5/smr/lem/Done','/b/data/lem');
            elseif strncmp(filename,'/Volumes/bgc6/smr/lem/Done',15)
                tgt = strrep(filename,'/Volumes/bgc6/smr/lem/Done','/b/data/lem/');
            elseif strncmp(filename,'/Volumes/bgc6/smr/lem/',15)
                tgt = strrep(filename,'/Volumes/bgc6/smr/lem/','/b/data/lem/');
            else
            tgt = strrep(filename,'/Users/bgc/tmp/uffs','/b/data');
            end
            d = dir(filename);
            if exist(tgt)
                c = dir(tgt);
                if c.datenum > d.datenum
                    ok = 1;
                    str = 'ds1 Newer';
                elseif c.datenum < d.datenum
                    str = sprintf('%.2f days Older',c.datenum-d.datenum);
                    ok = 0;
                else 
                    ok = 2; %date and size match.  No worries
                end
                if d.bytes == c.bytes
                    if showall >= ok
                        fprintf('%s matches %s (%s) %.1fkB %s\n',missing{j},tgt, str,d.bytes,str);
                    end
                elseif ok ==0 || (showall && ok ==1)
                    fprintf('%s %.1fkB  %.1fkB %d %s\n',missing{j},d.bytes./1024,c.bytes./1024,d.bytes-c.bytes,str);
                end
            elseif isempty(d)
                fprintf('%s No match\n',missing{j});
            else
                fprintf('%s %.1fkB\n',missing{j},d.bytes./1024);
            end
        end
    end
end

