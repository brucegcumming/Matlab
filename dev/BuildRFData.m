function rfs = BuildRFData(dirname, varargin)
%rfs = BuildRFData(dirname, ...)
%Reads .ufl/rf.mat  files and builds table for PlotMatp
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
    for j = 1:length(rfs)
        monk{j} = GetMonkeyName(rfs{j}.name);
        newrfs{j} = rfs{j}.name;
    end    
    monkeys = unique(monk);
    for j = 1:length(monkeys)
        outfile = ['/bgc/anal/' monkeys{j} '/allrfs.mat'];
        if exist(outfile)
            load(outfile);
            for k = 1:length(allrfs)
                oldnames{k} = allrfs{k}.name;
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
        if ~isempty(newid)
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
if isempty(d)
    return;
end

filename = dir2name(dirname,'filename');
ts = now;
nrf = 0;
dates = [];
depths = [];
StartDepth = [];
for j = 1:length(d);
    txt = scanlines(d(j).name);
    rid = find(strncmp('cm=rf',txt,5));
    for k = 1:length(rid)
        nrf = nrf+1;
        rf{nrf} = sscanf(txt{rid(k)},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(nrf,1:length(rf{nrf})) = rf{nrf};
    end
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
monkey = GetMonkeyName(dirname);
pen = ReadPen(FindPenLog(monkey,pe),'noplot');
if isfield(pen,'files') && isfield(pen,'visualarea')
    areas = {};
    nf=0;
    for j = 1:length(pen.files)
        if strfind(pen.files{j},filename)
            nf = nf+1;
            areas{nf} = pen.visualarea{j};
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
therf.name = dirname;
therf.readtime = [ts now];
try
save(rffile,'therf');
catch
    cprintf('errors','Error Writing %s\n',rffile);
    delete(rffile); %don't leave corrupt/empty file on disk
    therf.saveerror = 1;    
end
