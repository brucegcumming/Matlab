function res = cpfiles(src, tgt, varargin)
%cpfiles (srcdir, tgtdir,...)
%finds and lists all files in srdir that are not in tgtdir, or are older in tgtdir
%cpfiles (srcdir, tgtdir,'copynew')
%then copies the files from srcdir to tgtdir that are not in tgtdir
%cpfiles (srcdir, tgtdir,'copynewer')
%then copies the files from srcdir to tgtdir that are newer in srcdir
%
%cpfiles (srcdir, tgtdir,'hidenewer')
%lists new files, but not newer

go = [0 0 0];
show = [1 1 0];
fargs = {};
res = [];
excludes = {};
includes = {};
tgts = {};
mindiff = 1./24;  %1 hour by defauls
recurse = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'copynewer',8)
        go(2) = 1;
    elseif strncmpi(varargin{j},'copynew',7)
        go(1) = 1;
    elseif strncmpi(varargin{j},'exclude',7)
        j = j+1;
        excludes = {excludes{:} varargin{j}};
    elseif strncmpi(varargin{j},'include',7)
        j = j+1;
        includes = {includes{:} varargin{j}};
    elseif strncmpi(varargin{j},'smrdir',6) %defaults for syncing an SMR data directory
%        excludes = {excludes{:} '\.smr' '\.S2R' '\.[0-9]*\.mat' '\.lfp\.mat'};
        excludes = {excludes{:} '\.smr' '\.S2R' '\.lfp\.mat'};
    elseif strncmpi(varargin{j},'mindiff',4)
           j = j+1;
        mindiff = varargin{j};
    elseif strncmpi(varargin{j},'dest',4)
        j = j+1;
        tgts = {tgts{:} varargin{j}};
    elseif strncmpi(varargin{j},'nonrecursive',5)
        fargs = {fargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'name',4)
        fargs = [fargs varargin(j) varargin(j+1)]
    elseif strncmpi(varargin{j},'shownewer',8)
        show(2) = 1;
    elseif strncmpi(varargin{j},'hidenewer',10)
        show(2) = 0;
    elseif strncmpi(varargin{j},'hidenew',7)
        show(1) = 0;
    end
    j = j+1;
end

if isfield(src,'srcdir')  %a strcut with previous list
    res = src;
    sfiles = res.srcfiles;
    tfiles = res.tgtfiles;
    src = res.srcdir;
    tgt = res.tgtdir;
    tdates = res.tgtdates;
    sdates = res.srcdates;
else
    [sfiles, ssz, sdates, sex] = TreeFind(src,fargs{:});
    [tfiles, tsz, tdates, tex] = TreeFind(tgt,fargs{:});
    for j = 1:length(tfiles)
        tfiles{j} = strrep(tfiles{j},tgt,'');
    end
    for j = 1:length(tgts)
        [atfiles, atsz, atdates] = TreeFind(tgts{j},fargs{:});
        for k = 1:length(atfiles)
            atfiles{k} = strrep(atfiles{k},tgts{j},'');
        end
        tfiles = {tfiles{:} atfiles{:}};
        tsz = [tsz atsz];
        tdates = [tdates atdates];
    end

    res.srcfiles = sfiles;
    res.tgtfiles = tfiles;
    res.srcdir = src;
    res.tgtdir = tgt;
    res.srcdates = sdates;
    res.tgtdates = tdates;
    res.excludes = excludes;
end

if length(includes)
    gid = [];
    for k = 1:length(includes)
        for j = 1:length(sfiles)
            if regexp(sfiles{j},includes{k})
                gid = [gid j];
            end
        end
    end
    sfiles = sfiles(gid);
    sdates = sdates(gid);
end
exclid = [];
for j = 1:length(sfiles)
    if sex(j).isdir
            exclid = [exclid j];
    end
end
for k = 1:length(excludes)
    for j = 1:length(sfiles)
        if regexp(sfiles{j},excludes{k})
            exclid = [exclid j];
            fprintf('Excluding %s\n',sfiles{j});
        end
    end
end
if length(exclid)
    exclid = unique(exclid);
    gid = setdiff(1:length(sfiles),exclid);
    sfiles = sfiles(gid);
    sdates = sdates(gid);
end


for j = 1:length(sfiles)
    sfiles{j} = strrep(sfiles{j},src,'');
    ti = strmatch(sfiles{j},tfiles,'exact');
    if isempty(ti)
        if show(1)
            fprintf('%s is new\n',sfiles{j});
        end
        if go(1)
            fprintf('Copying %s to %s\n',[src '/' sfiles{j}],[tgt '/' sfiles{j}]);
            copyfile([src '/' sfiles{j}],[tgt '/' sfiles{j}]);
        end
    else
        if sdates(j) > tdates(ti) +mindiff
            if show(2)
                fprintf('%s is newer (%s vs %s)\n',sfiles{j},datestr(sdates(j)),datestr(tdates(ti)));
            end
            if go(2)
                fprintf('Copying %s to %s\n',[src '/' sfiles{j}],[tgt '/' sfiles{j}]);
                copyfile([src '/' sfiles{j}],[tgt '/' sfiles{j}]);
            end
        elseif sdates(j) <= tdates(ti) && show(3)
            fprintf('%s is older/same\n',sfiles{j});
        end
    end
end
            
