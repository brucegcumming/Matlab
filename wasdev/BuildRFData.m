function rfs = BuildRFData(dirname, varargin)
%Reads .ufl/rf.mat  files and builds table for PlotMatp
%rfs = BuildRFData(dirname, ...)
%BuildRFData(dirname, 'save') appends any new rsf to
% /bgc/anal/monkey/newrfs.mat
% if dirname is a cell string array, searches each directory


recurse = 0;
saverfs = 0;
rebuild = 0;
verbose = 0;
rfs = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'find',4)
        recurse = 1;
    elseif strncmpi(varargin{j},'rebuild',4)
        rebuild = 1;
    elseif strncmpi(varargin{j},'save',4)
        saverfs = 1;
    elseif strncmpi(varargin{j},'verbose',4)
        verbose = 1;
    end
    j = j+1;
end

if iscellstr(dirname)
    nc = 0;
    for j= 1:length(dirname)
        if verbose
            fprintf('Reading %s\n',dirname{j});
        end
        rf = GetRFFromDir(dirname{j}, rebuild);
        if ~isempty(rf)
            nc = nc+1;
            rfs{nc} = rf;
        end
    end
elseif recurse
    names = {};
    d = mydir(dirname);
    if d(j).isdir
        names = {names{:} d(j).name};
    end
    rfs = BuildRFData(dirname,varargin{:});
else
    rfs = GetRFFromDir(dirname, rebuild);
end

if saverfs
    if length(rfs) ==1 && ~iscell(rfs)
        X{1} = rfs;
        rfs = X;
    end
    for j = 1:length(rfs)
        monk{j} = GetMonkeyName(rfs{j}.name);
        [a,b,c,newrfs{j}] = GetMonkeyName(rfs{j}.name);
    end    
    monkeys = unique(monk);
    for j = 1:length(monkeys)
        outfile = [GetFilePath('anal') '/' monkeys{j} '/allrfs.mat'];
        if exist(outfile)
            updated = [];
            X=load(outfile);
            if isfield(X,'allrfs')
                allrfs = X.allrfs;
            else
                allrfs = X.rfs;
            end
            for k = 1:length(allrfs)
                [a,b,c, oldnames{k}] = GetMonkeyName(allrfs{k}.name);
                id = find(strcmp(oldnames{k},newrfs));
                if ~isempty(id) && allrfs{k}.readtime(1) < rfs{id(1)}.readtime(1) %new file
                    fprintf('Updating %s (%s)\n',oldnames{k},newrfs{id(1)});
                    allrfs{k} = rfs{id(1)};
                    updated(end+1) = id(1);
                end
            end
            [newnames, newid] = setdiff(newrfs,oldnames);
            for k = 1:length(newnames)
                fprintf('Adding %s to %s\n',newnames{k},outfile);
            end
            allrfs = {allrfs{:} rfs{newid}};
        else 
            allrfs = rfs;
            newid = 1;
        end
        if ~isempty(newid) || ~isempty(updated)
            save(outfile,'allrfs');
        end
    end
end

function therf = GetRFFromDir(dirname,rebuild)
therf = [];

rffile = dir2name(dirname, 'rf');
d = dir(rffile);
if exist(rffile,'file') && ~rebuild    
    try
    load(rffile);
    return;
    catch ME
        cprintf('red','Error Loading %s\n',rffile);
    end
end

d = mydir([dirname '/*.ufl']);
%remove backup files etc
good = [];
for j = 1:length(d)
    if ismember(d(j).filename(1),'#.')
        good(j) = false;
    else
        good(j) = true;
    end
end
if isempty(good)
    return;
end
d = d(find(good));

filename = dir2name(dirname,'filename');
ts = now;
nrf = 0;
dates = [];
depths = [];
StartDepth = [];
rfs = [];

for j = 1:length(d);
    txt = scanlines(d(j).name);
    vid = find(strncmp('VisualArea',txt,7));
    if ~isempty(vid)
        areas{j} = txt{vid(1)}(12:end);
    else
        areas{j} = '';
    end
    rid = find(strncmp('cm=rf',txt,5));
    rfn(j,1) = size(rfs,1)+1;
    for k = 1:length(rid)
        nrf = nrf+1;
        rf{nrf} = sscanf(txt{rid(k)},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(nrf,1:length(rf{nrf})) = rf{nrf};
    end
    rfn(j,2) = size(rfs,1);
    if ~isempty(regexp(d(j).name,'[M][0-9][0-9]'))
        types(j) = 1;
    else
        types(j) = 2;
    end
    did = find(strncmp('Created',txt,7));
    if ~isempty(did) && length(txt{did(end)}) > 13
        dates(nrf) = datenum(txt{did(end)}(9:end));
    end
    idxfile = strrep(d(j).name,'.ufl','idx.mat');
    if exist(idxfile)
        load(idxfile);
        if isfield(Expt.Trials,'ed')
            depths(nrf) = median(Expt.Trials.ed(Expt.Trials.ed > 0));
        end
        if isfield(Expt.Trials,'StartDepth')
            sd = Expt.Trials.StartDepth;
            StartDepth = median(sd(sd>0));
        end
    end
    if isempty(dates > 0)
        matfile = strrep(d(j).name,'.ufl','.mat');
        dates(nrf) = CreationDate(matfile);
    end
end
uflva = unique(areas);
if length(uflva) > 1
    fid = find(strcmp(uflva{1},areas));
    grfs = rfn(fid,1):rfn(fid,2);
    rfs = rfs(grfs,:);
end

d = mydir([dirname '/*idx.mat']);
areas = {};
for j = 1:length(d)
    load(d(j).name);
    if isfield(Expt,'Comments')
    id = find(strncmp('cm=VisualArea',Expt.Comments.text,10));
    if ~isempty(id)
        areas{j} = strrep(Expt.Comments.text{id(end)}(15:end),' .*','');
    end
    end
end
if isempty(areas)
    va = '';
else
    va = unique(areas);
end
gid = find(rfs(:,6)) > 0;
pe = median(rfs(gid,6));
[monkey, mname, mdir] = GetMonkeyName(dirname);
pen = ReadPen(FindPenLog(monkey,pe),'noplot');
if isfield(pen,'files') && isfield(pen,'visualarea')
    areas = {};
    nf=0;
    for j = 1:length(pen.files)
        if strfind(pen.files{j},filename)
            nf = nf+1;
            areas{nf} = deblank(pen.visualarea{j});
        end
    end
    if nf
        vb = unique(areas);
        if isempty(va)
            va = vb;
        elseif ~strcmp(va,vb)
            mycprintf('blue','Area Assignment mismatch for %s: %s in idx, %s in pen log\n',filename, va{:}, vb{:});
        end
    end
else
    fprintf('No Visual area data for %s (pe %.0f)\n',filename,pe);
end
if isfield(pen,'enterdepth') && isempty(StartDepth)
    StartDepth = pen.enterdepth;
end


addfile = [dirname '/' monkey mdir 'Add.txt'];
addtxt = scanlines(addfile,'silent');
nva = 0;
for j = 1:length(addtxt)
    if ~isempty(strfind(addtxt{j},'VisualArea'))
        fprintf('Using %s from %s\n',addtxt{j},addfile);
        nva = nva+1;
        va{nva} = regexprep(addtxt{j},'.*VisualArea=','');
        va{nva} = regexprep(va{nva},'\s.*','');
    end
end

if isempty(va)
    therf.area = 'unknown';
else
    therf.area = va{1};
end
dates = dates(dates>0);
if isempty(dates)
    d = mydir([dirname '/*.online']);
    if ~isempty(d)
        txt = scanlines(d(1).name,'bufsize',9192);
        id = find(strncmp('Reopened',txt,8));
        if ~isempty(id)
            dates(1) = datenum(txt{id(1)}(14:end));
        end
        id = find(strncmp('ed',txt,2));
        id = find(strncmp('StartDepth',txt,8));
        if ~isempty(id)
            a = sscanf(txt{id(1)},'StartDepth%d');
            StartDepth = [StartDepth a];
        end
    end
end
if ~isempty(dates)
    therf.date = median(dates(dates>0));
end
depths = depths(depths > 0);
if ~isempty(depths)
    therf.depth = median(depths);
else
    therf.depth = NaN;
end
StartDepth = unique(StartDepth);
if ~isempty(StartDepth)
    therf.StartDepth = StartDepth;
end
type = prctile(types,50);
if type == 1
    therf.electrode= 'uProbe';
else
    therf.electrode = 'Normal';
end
if size(rfs,1) > 2
    therf.rf = prctile(rfs,50);
else
    therf.rf = rfs(end,:);
end

if isfield(pen,'pos') && sum(abs(pen.pos)) > 0
%pen location setting in log overrides what was written at run time,
%bur record old value;
    if therf.rf(7) ~= pen.pos(1) || therf.rf(7) ~= pen.pos(1)  
        thetf.peninsmr = therf.rf(7:8);
        therf.rf(7:8) = pen.pos(1:2);
    end
end
therf.name = dirname;
therf.readtime = [ts now];
try
    fprintf('Saving %s\n',rffile);
    save(rffile,'therf');
catch
    cprintf('errors','Error Writing %s\n',rffile);
    delete(rffile); %don't leave corrupt/empty file on disk
    therf.saveerror = 1;    
end
