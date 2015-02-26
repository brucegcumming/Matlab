function res = BuildAllFullV(name, varargin)
%Builds FullV Matlab files from component files exported by spike 2. 
%BuildallFullV(dir, .....
%For all expts in dir. Also calls AllVPcs to cut clusters automatically.
%
% BuidlAllFullV(dir, 'scan')
% checks at the end to see if any more files have been made by spike2, and
% keeps going if new files exist. If you made RefClusters then the Typical use 
% for the end of a day is
% 
% BuidlAllFullV(dir, 'scan','refcutauto','matchcounts')
%
% BuidlAllFullV(dir, 'byname') Finds existing FullV files by name, and uses
% only these
% BuidlAllFullV(dir, 'expts', exlist)   %just does listed expts.
%
%BuildAllFullV(..., 'recalc')
%
%BuildAllFullV(..., 'reindex') forces rebuilding of index files (from .spkblk.)
%                     in case you have remade/modified spkblk files
%BuildAllFullV(..., 'recalcall') forces rebuilding of index files (from .spkblk.)
%
%BuildAllFullV(..., 'nocut') make FullV files but don't do auto
%cutting.(E.G. if Clustering done already)
%BuildAllFullV(..., 'refcut') appplies ref cut
%BuildAllFullV(..., 'refcutauto') appplies refcut AND does auto cut
%
%BuildAllFullV(..., 'remakeV')    forces rebuilding of FullV file
%
%BuildAllFullV(..., 'expts', explist, 'remakeV','recalc')    forces rebuilding of FullV file,
%and re-indexing of spkblk files, for just the named experiments
%
%...,'automode',mode)  passes the mode argument to AllVPcs to determine the
%                method used for automatic cutting
%to build a FullV file it must first build an index of what time samples
%are available for what probes in all the .mat file made by spike2.
%jan 2012  Changed so that default saved format uses int16 rather than
%double
%see also MakeProbeIndex, 


bysuffix = 0;
scanfiles = 0;
recalc = 0;
reindex = 0;
X.checkv = 0;
sizecheck = 0;  %set to max safe bytes for FullV to control if FullV is kept in memory
runautocut = 2; %run cut after build
memsize = CheckPhysicalMemory;
sizecheck = (memsize * 1e6) .* 0.8;
X.spkrate = 50;
template = [];
args = {'ndives' 0}; %don't allow chasing of trigger by default
version = 1.0;
parallel = 0;
autocutmode = 'mahal';
checkfiles = 1;

if isfield(name,'fullvdata') %previous result
    CheckErrors(name);
    return;
elseif isdir(name)
%    idxfile = [name / 'probes.mat'];
    idxfile = strrep(name,'.mat','idx.mat');
    bysuffix = 1;
    logfile = [name '/BuildLog.txt'];
else
    idxfile = strrep(name,'.mat','idx.mat');
    logfile = strrep(name,'.mat','log.mat');
end
expts = [];
byname = 0;
logfid = fopen(logfile,'a');
X.logfid = logfid;
X.nocut = 0;

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'V')
            FullV = varargin{j};
        end
    elseif strncmpi(varargin{j},'automode',7)
        j = j+1;
        autocutmode = varargin{j};
    elseif strncmpi(varargin{j},'byname',3)
        byname = 1;
    elseif strncmpi(varargin{j},'bysuffix',3)
        byname = 1;
    elseif strncmpi(varargin{j},'checkv',6)
        X.checkv = 1;
    elseif strncmpi(varargin{j},'expts',3)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'nocut',4)
        runautocut = 0;
        X.nocut = 1;
    elseif strncmpi(varargin{j},'quickautocut',4)
        runautocut = 2;
        args = {args{:} 'quickautocut'};
    elseif strncmpi(varargin{j},'refcut',4) %includes 'refcluster
        if strncmpi(varargin{j},'refcutauto',8) ||  strncmpi(varargin{j},'refclusterauto',12)
            runautocut = 2;
        else
            runautocut = 0;
        end
        X.nocut = 2;
        args = {args{:} 'reapply' 'refclusters'};
    elseif strncmpi(varargin{j},'parallel',6)
        parallel = 1;
    elseif strncmpi(varargin{j},'reindex',5)
        reindex = 1;
    elseif strncmpi(varargin{j},'remakeV',7)
        recalc = 4;
        runautocut = 0;
    elseif strncmpi(varargin{j},'recalcsuffix',8)
        recalc = 5;
    elseif strncmpi(varargin{j},'recalcall',8)
        recalc = 3;
    elseif strncmpi(varargin{j},'rebuild',5)
        recalc = 2;
    elseif strncmpi(varargin{j},'recalc',3)
        if recalc == 4
            recalc = 5;
            reindex = 1;
        elseif byname %just redo auto cuts
            recalc = 1;
        else
            recalc = 2;
        end
    elseif strncmpi(varargin{j},'scan',4)
        scanfiles = 1;
        runautocut = 2;
    elseif strncmpi(varargin{j},'template',6)
        args = {args{:} varargin{j:j+1}};
        j = j+1;
        template = varargin{j};
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
res.starts(1) = now;
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:MaxIterations');

if isdir(name)
    path = name;
else
[path, fname] = fileparts(name);
end
X.catcherrs = 0;
if exist('FullV','var')
    Expts(1)  = 1;
    probelist = 1:size(FullV.V,1);
    expts = 1;
else
    if ~exist(idxfile,'file') & ~isdir(name)
        fprintf('No Index File %s\n',idxfile);
        return;
    end
    if ~isdir(name)
        if ~exist(idxfile)
            fprintf('%s not a directory, and no index file. Check Path',name);
            return;
        end
        load(idxfile);
        if exist('Expt','var') && isfield(Expt,'DataType') && strmatch('GridData',Expt.DataType)
            ProcessGridFullV(name);
        end
    end
    

X.recalc = recalc;
X.reindex = reindex;
X.sizecheck = sizecheck;
X.runautocut = runautocut;

ns = 1;
res.starts(2) = now;
res.name = name;
if byname
    d = dir(name);
    for j = 1:24; probes(j).probe = j; end
    for j = 1:length(d)
        sid = regexp(d(j).name,'Expt[0-9]*FullV\.');
        if length(sid)
            ex = sscanf(d(j).name(sid(1)+4:end),'%d');
            if isempty(expts) || ismember(ex,expts)
                fprintf(logfid,'Processing %s %s\r\n',d(j).name,datestr(now));
               [res.cls{ns}, details] =  ProcessSuffix(path, probes, ex, 0, X, args);
               res.loadtimes(ns) = details.loadtime;
               ns = ns+1;
            end
        end
    end
    return;
elseif scanfiles
    d = dir(name);
    newfiles = 1;
    ts = now;
    [probes, errs] = MakeProbeIndex(path,'recalc','noplot');
    fprintf(logfid,'Indexing took %.2f at %s\r\n',mytoc(ts),datestr(now));
    
    donesuff = [];
    while newfiles
        suffixes = unique([probes.suffix]);
        %dont make the last one - may not be finished yet...
        suffixes = suffixes(1:end-1);
        newsuff = setdiff(suffixes,donesuff);
        if length(newsuff) > 0
            fprintf(logfid,'Suffixes %s\n',sprintf('%d ',suffixes));
            fprintf(logfid,'NewSuff %s\n',sprintf('%d ',newsuff));
            fprintf(logfid,'DoneSuff %s\n',sprintf('%d ',donesuff));
        end
        for j = 1:length(newsuff);
            fprintf(logfid,'%s Processing %s %d %s\r\n',datestr(now),path,newsuff(j),datestr(now));
            res.cls{newsuff(j)} = ProcessSuffix(path, probes, newsuff(j), 0, X, args);
            donesuff = [donesuff newsuff(j)];
        end
        ti = now;
        probes = MakeProbeIndex(path,'recalc','noplot','quiet');
        suffixes = unique([probes.suffix]);
        if length(newsuff) > 1 %did make some new ones this time
        fprintf('Suffixes done: %s     ',sprintf(' %d',donesuff));
        fprintf('Suffixes ready %s\n',sprintf(' %d',suffixes));
        end
        fprintf(logfid,'Indexing took %.2f at %s. New Suff %s\r\n',mytoc(ti),datestr(now),sprintf(' %d',newsuff));
        %if there is only on file left, and no new spkblks have been made,
        %finish
        if length(donesuff) && max(donesuff) >= max(suffixes)-1 && now-max([probes.filetime]) > 0.01
            fprintf('Finished\n');
            newfiles = 0;
        end
    end
    fprintf(logfid,'Processing %s %d %s (Last)\r\n',path,max(suffixes),datestr(now));
    res.cls{max(suffixes)} = ProcessSuffix(path, probes, max(suffixes), 0, X, args);
    if runautocut == 2 || X.nocut ==2
        X.runautocut = 1;
        X.recalc = 0;
        X.nocut = 0;
        res.startcuts = now;
        if parallel
            parfor (ex = donesuff)
                ProcessSuffix(path, probes, ex,0, X,args);
            end
        else
            for ex = donesuff
                ProcessSuffix(path, probes, ex,0, X,args);
            end
        end
    end
    fullv.Check(name);
return;
else

    if reindex
        probes = MakeProbeIndex(path,'recalc');
    elseif length(expts) && recalc == 4 %take probe ids from exist
        probes = MakeProbeIndex(path,'suffix',expts);
    elseif length(expts) && recalc == 5  %reread spkblk
        probes = MakeProbeIndex(path,'suffixrecalc',expts);
        res.probes = probes;
    else
        probes = MakeProbeIndex(path);
    end
    if ~isfield(probes,'probe')
        cprintf('red','Missing Probe data in Index\n');
        return;
    end
    probelist = unique([probes.probe]);
    if isempty(expts)
        if isfield(probes,'suffix')
            expts = unique([probes.suffix]);
        else
            expts = 1:length(expts);
        end
    end
end
end


if parallel
    X.logfid = -1;
    X.catcherrs = 1;
    cls = {};
    parfor (j = 1:length(expts))
        addpath('/b/bgc/matlab/dev');
        ex = expts(j);
        starts(j) = now;
        err = 0;
        id = find([probes.suffix] == ex);
        eprobes = unique([probes(id).probe]);
        if length(eprobes) < length(probelist)
            fprintf('Expt %d Has only %d probes\n',ex,length(eprobes));
            err = 1;
        end
        outname = [path '/Expt' num2str(ex) 'FullV.mat']
        if err > 0
        elseif bysuffix
            fprintf('%s from suffix\n',outname);
            cls{j} = ProcessSuffix(path, probes, ex,0, X,args);
        elseif exist(outname,'file') && recalc < 2
                fprintf('Loading %s\n',outname);
                load(outname);
        else
            fprintf('Making FullV for Expt %d\n',ex)
            ts = now;
            if bysuffix
                cls{j} =  ProcessSuffix(path, probes, ex,0, X, args);
            else
                
            end
            ProcessSuffix(path, probes, ex,0, X,args);
        end
    end
    res.fullvdata(expts)=cls;
    res.starts(expts) = starts;
else
    X.catcherrs = 0;        
    for ex = expts
        res.starts = [res.starts now];
        err = 0;
        id = find([probes.suffix] == ex);
        eprobes = unique([probes(id).probe]);
        outname = [path '/Expt' num2str(ex) 'FullV.mat'];
        if isempty(id) && exist(outname)
            fprintf('%s exists, but not in probes.mat\n',outname);
            
        elseif length(eprobes) < length(probelist)
            fprintf('Expt %d Has only %d probes\n',ex,length(eprobes));
            err = 1;
        end
        if ~exist('FullV','var')
            if err > 0
            elseif bysuffix
                res.fullvdata{ex} = ProcessSuffix(path, probes, ex,0, X,args);
            elseif exist(outname,'file') && recalc < 2
                fprintf('Loading %s\n',outname);
                load(outname);
            else
                fprintf('Making FullV for Expt %d\n',ex)
                ts = now;
                if bysuffix
                    res.fullvdata{ex} =  ProcessSuffix(path, probes, ex,0, X, args);
                else
                    t = [Expts(ex).start Expts(ex).end]./10000;
                    files = MakeProbeIndex(path,'readfiles',t);
                    FullV = PlotSpikeC(files,5,'probes',probelist,'sumv','submean','makev');
                    FullV.name = name;
                    FullV.exptno = ex;
                    outname = [path '/Expt' num2str(ex) 'FullV.mat'];
                    FullV.buildtime = mytoc(ts);
                    FullV.builddate = now;
                    res.cls{ex}.loadtime = FullV.buildtime;
                    save(outname,'FullV','-v7.3');
                    clear files;
                end
                if bysuffix == 0 && runautocut == 1
                    fprintf('Cutting Clusters for Expt %d\n',ex);
                    AllVPcs(FullV,'tchan',[1:size(FullV.V,1)],'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode',autocutmode,args,'savespikes','autocutall');
                    clear FullV;
                else
                    res.fullvdata{ex} = ProcessSuffix(path, probes, ex,0, X,args);
                    pack;
                end
            end
        end
    end
end


if runautocut == 2  || X.nocut == 2%do all teh cutting after building all the files
    X.runautocut = 1;
    res.startcuts = now;
    X.recalc = 0;
    X.nocut = 0;
    if parallel
        parfor ex = expts
            x = ProcessSuffix(path, probes, ex,0, X,args);
            if isempty(res.fullvdata{ex})
                fullvdata{ex} = x;
            end
        end
    else
        for ex = expts
            x = ProcessSuffix(path, probes, ex,0, X,args);
            if ex > length(res.fullvdata) || isempty(res.fullvdata{ex})
                res.fullvdata{ex} = x;
            end
        end
    end
end
CheckErrors(res);
res.end = now;
if logfid > 0
    try
        fclose(logfid);
    end
end

%CheckFullV elsewhere in here should really use fullv.check. But needs to
%be checked
if checkfiles %run a final check on all the files built
   fullv.Check(name);
end

function CheckErrors(res)

if ~isfield(res,'fullvdata')
    fprintf(' NO fulldata in Return result\n');
    return;
end
nerr = 0;
for j = 1:length(res.fullvdata)
    if isfield(res.fullvdata{j},'errmsg') || isfield(res.fullvdata{j},'error')
        fprintf('Errors in Expt%d\n',j);
        nerr = nerr+1;
    end
end

if nerr == 0
    fprintf('No errors\n')
end

function [res, details] = ProcessSuffix(path, probes, ex, lastblk, X, args)

details.loadtime = 0;
template = [];
mkint = 1;
readint = 0;
res = {};
exptargs = {};
refcut = 0;

j = 1;
checkexpts = 0;
while j <= length(args)
    if strncmpi(args{j},'check',5)
        checkexpts = 1;
        if strncmpi(args{j},'checkandbuild',8)
            checkexpts = 2;
        end
    elseif strncmpi(args{j},'mkdouble',5)
        mkint = 0;
    elseif strncmpi(args{j},'mkint',5)
        mkint = 1;
    elseif strncmpi(args{j},'noframes',5)
        exptargs = {exptargs{:} args{j}};
    elseif strncmpi(args{j},'refclusters',9)
        refcut = 1;
    elseif strncmpi(args{j},'refcut',6)
        refcut = 1;
    elseif strncmpi(args{j},'readint',5)
        readint = 1;
    elseif strncmpi(args{j},'template',6)
        j = j+1;
        template = args{j};
        cargn = j;
    end
    j = j+1;
end

if lastblk > 0
    outname = [path '/Expt' num2str(ex) 'aFullV.mat'];
else
    outname = [path '/Expt' num2str(ex) 'FullV.mat'];
end

if readint
    outname = strrep(outname,'FullV','FullVi');
end
probelist = unique([probes.probe]);
eid = find([probes.suffix] == ex);
d = dir(outname);
ClustersUptoDate = 0;
if isempty(eid) %can't make new FullV without this. But if it exists, its up do date
    FullvUptoDate = 1;    
elseif isfield(probes,'filetime') && (isempty(d) || d(1).datenum < max([probes(eid).filetime]))
    FullvUptoDate = 0;
else
    FullvUptoDate = 1;
end
if ~isempty(d)
    cname = regexprep(outname,'FullV','ClusterTimes');
    b = dir(cname);
    if ~isempty(b) && b(1).datenum > d(1).datenum
        fprintf('%s is up to date\n',cname);
        ClustersUptoDate = 1;
    end
end
Expt = LoadExpt(path, ex, X.recalc > 1, X.logfid);
if exist(outname,'file') && X.recalc < 2 && FullvUptoDate
    ts = now;
    if checkexpts == 0 && ClustersUptoDate || X.runautocut == 2 %Don't need to load FullV
        fprintf('Dont need to load %s\n',outname);
        return;
    end
    fprintf('Loading %s\n',outname);
    d = dir(outname);
    clear FullV;
    load(outname);
    if isfield(FullV, 'intscale')
        FullV.V = double(FullV.V) .* FullV.intscale(1)/FullV.intscale(2);
    end
    FullV.loadname = outname;
    if isfield(FullV,'matfile')
        [a,b] = fileparts(FullV.matfile);
        [c,d] = fileparts(outname);
        FullV.matfile = [c '/' b '.mat']; 
    end
    res.loadtime = mytoc(ts);
%    FullV.name = outname;
    if X.logfid > 0
        fprintf(X.logfid,'Load Took %.1f\r\n',res.loadtime);
    end
else
    ts = now;
%    Expt = LoadExpt(path, ex, X.recalc > 1, X.logfid);
    FullV = BuildFullVFile(path, probes, ex, lastblk, outname, Expt);
    res.loadtime = mytoc(ts);
    res = CopyFields(res, FullV, {'errmsg' 'error'},'-noempty');
    if X.logfid > 0
        if isfield(FullV,'error') && FullV.error > 0
            fprintf(X.logfid,'File %s not made: error %d\r\n',outname,FullV.error);
            return;
        else
            fprintf(X.logfid,'Build %s  Took %.1f\r\n',outname,res.loadtime);
        end
    end
end
x = whos('FullV');

fprintf('File Load/Build tood %.2f (%s) size %.2fGb\n',res.loadtime,outname,x.bytes./(1024 * 1024 * 1024));
details.loadtime = res.loadtime;
[~, FullV] = CheckFullV(FullV,Expt);
if mkint && ~isfield(FullV,'intscale')
    if ~isfield(FullV,'V')
        cprintf('red','%s Missing V!!\n',outname);
        return;
    end
    for j = 1:size(FullV.V,1)
        vm(j) = max(abs(FullV.V(j,:)));
    end
    vm = max(vm);
    FullV.V = int16(FullV.V .* 32700/vm);
    FullV.intscale = [vm 32700];
    if mkint == 2 %save to separate filename for comparisons
        outname = strrep(outname,'FullV','FullVi');
        save(outname,'FullV','-v7.3');
    return;
    end
    save(outname,'FullV','-v7.3');
end

if readint
    return;
end
if checkexpts
    if isempty(Expt)
        fprintf('Can''t read Expt %d\n',ex);
    else
        if checkexpts == 2
            [rebuild, FullV] = CheckFullV(FullV,Expt,'rebuild');
            if rebuild
                FullV = BuildFullVFile(path, probes, ex, lastblk, outname, []);
            end
        else
            [rebuild, FullV] = CheckFullV(FullV,Expt);
        end
    end
    res = CopyFields(res,FullV,{'errmsg' 'errdata'});
    return;
end
res = CopyFields(res,FullV,{'errmsg' 'errdata'});
if ismember(X.recalc, [4 5]) %just make FullV
    return;
end
if isfield(FullV,'matfile')
    Expt = LoadExpt(FullV.matfile, ex, X.recalc > 1, X.logfid, exptargs{:});
else
    Expt = LoadExpt(path, ex, X.recalc > 1, X.logfid);
%if Expt is empty, next step will fail.      
%    Expt = [];
end
if isfield(Expt,'Trials')
    ts = Expt.Trials(1).Start./10000;
    te = Expt.Trials(end).End./10000;
    if FullV.t(1) > ts || FullV.t(end) < te
        fprintf('ERROR!!!: FullV time range < Expt\n');
        if X.logfid > 0
            fprintf(X.logfid,'ERROR!!!: FullV time range (%.2f) < Expt (%.2f)\n',FullV.t(end),te);
        end
    end
    fprintf('Ex %d: %d Ch x %d samples %.1f (%.1f)- %.1f(%.1f)\n',ex,size(FullV.V,1), size(FullV.V,2),FullV.t(1),ts,FullV.t(end),te);
else
    fprintf('Expt %d has no trials\n',ex);
    return;
end
if X.nocut == 1
    return;
end
fprintf('Cutting Clusters for Expt %d\n',ex);
if isfield(FullV,'lastblk')
    lastblk = FullV.lastblk;
else
    lastblk = 0;
end
if X.sizecheck > 0
    v = whos('FullV');
    if v.bytes > X.sizecheck
        clear FullV;
        FullV = outname;
    end
end
if refcut && X.nocut == 0%make refcut if Clustertimes is out of date/non-existent
    if X.runautocut == 1 %asked for both
        res = AllVPcs(FullV,'tchan',probelist,'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal','ndives',0,'savespikes','noninteractive','autocutall');
    end
    res = {};
    for j = probelist
        if X.catcherrs
        try
        res{j} = AllVPcs(FullV,'tchan',j,'savespikes','noninteractive',args{:});
        catch ME
            res{j}.errstate = ME;
            if isfield(FullV,'loadname')
            res{j}.fullvfile = FullV.loadname;
            end
        end
        else
            res{j} = AllVPcs(FullV,'tchan',j,'savespikes','noninteractive',args{:});
        end
    end
    res{j} = CopyFields(res{j}, FullV, {'errmsg' 'error'},'-noempty');
    clear ms;



elseif ~isempty(template)
    for j = probelist
        res{j} = AllVPcs(FullV,'tchan',j,'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal','logfid',X.logfid,args{:},'noninteractive');
        res{j} = CopyFields(res{j}, FullV, {'errmsg' 'error'},'-noempty');
    end
    clear ms;
elseif X.runautocut == 1
    res = AllVPcs(FullV,'tchan',probelist,'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal',args{:},'savespikes','noninteractive','autocutall');
    res = CopyFields(res, FullV, {'errmsg' 'error'},'-noempty');
else 
    res = {};
end
if lastblk
    PrintMsg(X.logfid,sprintf('%s Too big for one file',outname));
    clear FullV
    ProcessSuffix(path, probes, ex, lastblk, X, args);
end
   
function FullV = BuildFullVFile(path, probes, ex, lastblk, outname, Expt)

version = 1.0;
probelist = unique([probes.probe]);
ts = now;
id = find([probes.suffix] == ex);
Expt = LoadExpt(path, ex, 1, 0);
matfile = [];
for j = 1:length(id)
    mfiles{j} = regexprep(probes(id(j)).file,'\.spkblk[0-9]*\.mat','.mat');
    mfiles{j} = regexprep(mfiles{j},'A.([0-9]*).mat','.$1.mat');
end
mfiles = unique(mfiles);
matfile = [path '/' mfiles{1}];
if isfield(Expt,'Trials')
    args = {'Expt' Expt};
else
    args = {};
end
    
FullV = PlotSpikeC(probes(id),5,'probes',probelist,'sumv','submean','makev','lastblk',lastblk,'prefix',path,args{:});
if isdir(path)
    FullV.name = path;
end
FullV.exptno = ex;
FullV.buildtime = mytoc(ts);
FullV.builddate = now;
FullV.buildversion = version;
FullV.matfile = matfile;

if isfield(FullV,'error') && FullV.error > 0
    return;
end
if ~isempty(Expt)
    [~,FullV] = CheckFullV(FullV, Expt);
end
        
fprintf('Writing %s at %s\n',outname,datestr(now));
save(outname,'FullV','-v7.3');

function [rebuild, FullV] = CheckFullV(FullV, Expt, varargin)
rebuild = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'rebuild',5)
        rebuild = 1;
    end
    j = j+1;
end

if ~isfield(Expt,'Trials')
    cprintf('red','No Completed Trials in Expt%d\n',GetExptNumber(FullV))
    return;
end
sds = std(FullV.V,[],2);
FullV.chstd = sds;

id = find(sds < 0.01); %Check for empty channels
if length(id)
    FullV = AddError(FullV, 'Expt %d Empty/Flat Channels %s\n', FullV.exptno, sprintf(' %d',id));
end
tid = [];
for j = 1:length(FullV.blkstart)
    t(1) = FullV.blkstart(j);
    t(2) = t(1) + FullV.blklen(j).*FullV.samper;
    t = t.*10000;
    id = find([Expt.Trials.Start] > t(1)& [Expt.Trials.Start] < t(2));
    tid = [tid id];
end
id = setdiff(1:length(Expt.Trials),tid);
if length(id)
    FullV = AddError(FullV,'Missing %d Trials for %s: %s\n',length(id),Expt.Header.name,sprintf(' %d',id));
    if rebuild
        name = strrep(Expt.Header.name,'idx.mat','.smr');
        cmd = ['C:\\Spike7\\sonview.exe /M ' name ' C:\\Spike7\\MakeSpikeWaves.s2s'] 
        system(cmd);
        name = regexprep(name,'\.[0-9]*\.smr','A$0');
        cmd = ['C:\\Spike7\\sonview.exe /M ' name ' C:\\Spike7\\MakeSpikeWaves.s2s'] 
        system(cmd);
    end
else
    rebuild = 0;
end

function Ex = LoadExpt(name, ei, rebuild, logfid, varargin)
        Ex = [];
        if isdir(name)
            [a,b,c,d] = GetMonkeyName(name);
            smrname = [name '/' a c];
            exfile = [smrname '.' num2str(ei) 'idx.mat'];
            matfile = [smrname '.' num2str(ei) '.mat'];
       elseif exist(name,'file')
                matfile = name;
                exfile = strrep(matfile,'.mat','idx.mat');
        else
            cprintf('red','%s does not exist\n',name);
            return;
        end
        if exist(matfile) && (~exist(exfile,'file') || rebuild)
            PrintMsg(logfid,'Building %s\n',exfile);
            APlaySpkFile(matfile,'bysuffix','noerrs', varargin{:});
        end
        if exist(exfile,'file')
            fprintf('Loading %s\n',exfile);
            load(exfile);
            id = find(Expt.Trials.Result == 1);
            for t = length(id):-1:1
                Ex.Trials(t).Start = Expt.Trials.Start(id(t));
                Ex.Trials(t).End = Expt.Trials.End(id(t));
                Ex.Trials(t).ed = Expt.Trials.ed(id(t));
                Ex.Trials(t).id = Expt.Trials.id(id(t));
                Ex.Trials(t).Trial = id(t);
            end
            if isempty(ExptList)
                Ex.Header.expname = 'None';
            else
                Ex.Header.expname = ExptList(end).expname;
            end
            Ex.Header.name = exfile;
            Ex.Stimvals.ed = mean(Expt.Trials.ed(id));
        else
            fprintf('Cant Find Expt %s\n',exfile);
        end