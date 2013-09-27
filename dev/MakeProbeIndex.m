function [probes, errs]  = MakeProbeIndex(dname, varargin)
%probes = MakeProbeIndex(path, ...
%read all spkblk files in directory 'path', and build an index file
%listing start/end times.
%
%     files = MakeProbeIndex(probes,'prefix',path,'readfiles',[t1 t2]) reads in
%
%all files needed span times between t1 and t2. With no time range, all
%files are read.
%
%If the index file has already been built, MakeProbeIndex just reads it off
%the disk.  To force rebuiling use
%
%     MakeProbeIndex(path,'recalc')
probes = [];
errs = [];
nc = 0;
readfiles = 0;
recalc = 0;
prefix = [];
trange = [];
res = [];
dosuffs = [];
readsuffix =0;
plottype = 0;
checkprobes  = 1;
readnames = 0;
rescan = 0;
nosave = 0;
verbose = 2;

j =1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'blkstart')
        res = varargin{j};
    elseif strncmpi(varargin{j},'nocheck',4)
        checkprobes = 0;
    elseif strncmpi(varargin{j},'noplot',6)
        plottype = 0;
    elseif strncmpi(varargin{j},'newprobes',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'prefix',4)
        j = j+1;
        prefix = varargin{j};
        nc = length(probes);
    elseif strncmpi(varargin{j},'quiet',4)
        verbose = 1;
    elseif strncmpi(varargin{j},'readfiles',8)
        readfiles = 1;
        plottype = 0;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            trange = varargin{j};
        end
    elseif strncmpi(varargin{j},'readnames',8)
        readnames = 1;
    elseif strncmpi(varargin{j},'suffix',6)
        if strncmpi(varargin{j},'suffixrecalc',10)
            recalc = 1;
            checkprobes = 1;
        end
        j = j+1;
        dosuffs = varargin{j};
        nosave = 1;
    elseif strncmpi(varargin{j},'readsuff',8)
        readfiles = 1;
        plottype = 0;
        checkprobes = 0;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            readsuffix = varargin{j};
            if readsuffix == 0
                fprintf('Suffix must be > 0\n');
                return;
            end
        end
    elseif strncmpi(varargin{j},'recalc',5)
        recalc = 1;
    elseif strncmpi(varargin{j},'reindex',5)
        recalc = 1;
    elseif strncmpi(varargin{j},'rescan',5)
        recalc = 1;
        rescan = 1;
        plottype = 0;
    elseif strncmpi(varargin{j},'update',5) %only read new or modified spkblks
        recalc = 0;
        rescan = 2;
        plottype = 0;
    end
    j = j+1;
end

if ischar(dname) &&  readfiles
    outfile = [dname '/probes.mat'];
    if exist(dname,'file')
        load(outfile);
        prefix = dname;
        dname = probes;
    end
end    
if isstruct(dname)  %%plot/check probe data
    if checkprobes
        errs = CheckProbes(dname, verbose);
    end
    if plottype ==1
    PlotProbes(dname,res);
    end
    if readfiles
        if readsuffix
            if readnames
                probes = ReadProbeSuffixNames(dname,prefix,readsuffix);
            else
                probes = ReadProbeSuffix(dname,prefix,readsuffix);
            end
        else
            probes = ReadProbeFiles(dname,prefix,trange);
        end
    end
    return;
end
d = dir([dname '/*spkblk*']);
outfile = [dname '/probes.mat'];
if exist(outfile,'file') && rescan
    load(outfile);
end
oldsuff = [];
if isfield(probes,'suffix')
    oldsuff = unique([probes.suffix]);
elseif exist(outfile,'file') && ~recalc && isempty(probes)
    fprintf('Reading existing index %s from disk\n',outfile);
    load(outfile);
    if checkprobes
    errs = CheckProbes(probes, verbose);
    end
    return;
end



subfiles = {};
files = {};
for j = 1:length(d)
    sid = regexp(d(j).name,'[0-9,A]\.[0-9]*\.spkblk');
    bid = regexp(d(j).name,'\.spkblk');
    if length(sid)
        subfiles = {subfiles{:} {d(j).name}};
        filesuff(length(subfiles)) = sscanf(d(j).name(sid(1)+2:end),'%d');
        blkn(length(subfiles)) = sscanf(d(j).name(bid(1)+7:end),'%d');
        ids(length(subfiles)) = j;
    elseif length(bid)
        files = {files{:} {d(j).name}};
        blkn(length(files)) = sscanf(d(j).name(bid(1)+7:end),'%d');
        filesuff(length(files)) = 0;
        ids(length(files)) = j;
    end
end

if rescan == 2  %updating any new spkblks
    [sfiles, aid] = unique({probes.file});
    sfiledate = [probes(aid).filetime];
    needid = [];
    oldid = []; %track which elemetns will be replaced
    for j = 1:length(ids)
        a = find(strcmp(d(ids(j)).name,sfiles));
        if length(a) == 1 && sfiledate(a) < d(ids(j)).datenum
            needid = [needid j];
            oldid = [oldid find(strcmp(sfiles{a},{probes.file}))];
        end
    end
    if isempty(needid)
        return;
    else
        filesuff = filesuff(needid);
        ids = ids(needid);
        if ~isempty(subfiles)
            subfiles = subfiles(needid);
        end
        if ~isempty(files)
            files = files(needid);
        end
        blkn = blkn(needid);
    end
    probes = probes(setdiff(1:length(probes),oldid));
    nc = length(probes);
    oldsuff = [];
end

if sum(filesuff > 0) == 0
    filesuff(filesuff == 0) = 1;
end
if isempty(dosuffs)
    dosuffs = unique(filesuff);
end
pos = filesuff.*1000+blkn;
[a,idlist] = sort(filesuff);

ltimes = [];
ts = now;
for m = 1:length(ids)
    j = ids(idlist(m));
    sid = strfind(d(j).name,'spkblk');
    if length(sid) & ~ismember(filesuff(idlist(m)),oldsuff) & ismember(filesuff(idlist(m)),dosuffs);
        fname = [dname '/' d(j).name];
        tl = now;
        timefix = 0;
        a = load(fname);
        if isfield(a,'timefix')
            timefix = a.timefix;
            cprintf('blue','Adjusting times in %s by %.2f\n',fname,timefix);
        end
        ltimes(m) = mytoc(tl);
        f = fields(a);
        filenames{j} = fname;
        for k = 1:length(f)
            C = a.(f{k});
            if sum(isfield(C,{'start' 'length' 'interval'})) == 3
            nc = nc+1;
            probe = sscanf(C.title,'Spike %d');
            probes(nc).probe = probe;
            probes(nc).start = a.(f{k}).start+timefix;
            probes(nc).end = timefix+C.start+C.length .* C.interval;
            probes(nc).var = f{k};
            probes(nc).file = d(j).name;
            probes(nc).fileid = sscanf(d(j).name(sid+6:end),'%d');
            probes(nc).suffix = filesuff(idlist(m));
            probes(nc).filetime = d(j).datenum;
            end
        end
    end
end
if verbose > 1
fprintf('Processing took %.1f (%.1f loading)\n',mytoc(ts),sum(ltimes));
end
if nosave == 0
save(outfile,'probes');
end
errs = CheckProbes(probes, verbose);
if plottype
    ts = now;
    PlotProbes(probes);
    fprintf('Plotting Took %.1f\n',mytoc(ts));
end

function errs = CheckProbes(probes, verbose)

errs = [];

if isfield(probes,'suffix')
    ts = now;
    exs = unique([probes.suffix]);
    probelist = unique([probes.probe]);
    np = length(probelist);
    for j = exs
        errs(j) = 0;
        id = find([probes.suffix] == j);
        nep = length(unique([probes(id).probe]));
        if nep < np
            fprintf('Expt %d only %d/%d probes\n',j,nep,np);
            errs(j) = 1;
        end
        for k = 1:nep
            pid = find([probes(id).suffix] == j & [probes(id).probe] == k);
            pid = id(pid);
            [starts, sid] = sort([probes(pid).start]);
            ends = [probes(pid(sid)).end];
            bid = find(starts(2:end) + 0.002 < ends(1:end-1));
            if length(bid)
                fprintf('End/Start Errors Expt %d P%d at %s\n',j,k,sprintf('%.2f ',starts(1+bid)));
                errs(j) = bitor(errs(j),2);
            end
        end
    end
    if verbose > 1
    fprintf('Checking Took %.1f\n',mytoc(ts));
    end
end


function files = ReadProbeFiles(p, prefix,trange)
if isempty(trange)
    trange = [min([p.start]) max([p.end])];
end
k = 0;
allfiles = {};
for j = 1:length(p)
    if p(j).start >= trange(1) && p(j).end <= trange(2)
        if isempty(strmatch(p(j).file,allfiles))
        k = k+1;
        fname = [prefix '/' p(j).file];
        files{k} = load(fname);
        files{k}.name = fname;
        allfiles{k} = p(j).file;
        end
    end
end

function files = ReadProbeSuffix(p, prefix,suffix)
k = 0;
allfiles = {};
for j = 1:length(p)
    if p(j).suffix == suffix
        if isempty(strmatch(p(j).file,allfiles))
        k = k+1;
        fname = [prefix '/' p(j).file];
        files{k} = load(fname);
        files{k}.name = fname;
        allfiles{k} = p(j).file;
        end
    end
end

function files = ReadProbeSuffixNames(p, prefix,suffix)
k = 0;
allfiles = {};
for j = 1:length(p)
    if p(j).suffix == suffix
        if isempty(strmatch(p(j).file,allfiles))
            k = k+1;
            fname = [prefix '/' p(j).file];
            files{k} = fname;
            allfiles{k} = p(j).file;
        end
    end
end


function PlotProbes(p, varargin)
res = [];
j= 1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'blkstart')
        res = varargin{j};
    end
    j = j+1;
end

hold off;
for j = 1:length(p); 
    plot([p(j).start p(j).end],[p(j).probe p(j).probe],'-'); 
    hold on; 
    if p(j).probe == 1
        plot([p(j).start p(j).start],[1 16],'r-');
        text(p(j).start,0.5,sprintf('%d',p(j).fileid));
    elseif p(j).probe == 24
        plot([p(j).start p(j).start],[17 24],'r-');
        text(p(j).start,16.5,sprintf('%d',p(j).fileid));
    end
end
y = 15.5;
if ~isempty(res)
    for j = 1:length(res.blkstart)
        plot([res.blkstart(j) res.blkstart(j)+res.blklen(j) .* res.samper],[y y],'m');
    end
end