function res = BuildAllFullV(name, varargin)
%BuildallFullV(dir, .....
%Builds FullV Matlab files from component files exported by spike 2. For
%all expts in dir. Then calls AllVPcs to cut clusters automatically.
%
% BuidlAllFullV(dir, 'scan')
% checks at the end to see if any more files have been made by spike2, and
% keeps going if new files exist.
% BuidlAllFullV(dir, 'byname') Finds existing FullV files by name, and uses
% only these
% BuidlAllFullV(dir, 'expts', exlist)   %just does listed expts.
%
%BuildAllFullV(..., 'recalc')
%
%BuildAllFullV(..., 'recalcall') forces rebuilding of index files (from .spkblk.)
%
%BuildAllFullV(..., 'nocut') make FullV files but don't do auto cutting.(E.G. if Clustering done already)
%
%BuildAllFullV(..., 'remakeV')    forces rebuilding of FullV file
%
%BuildAllFullV(..., 'expts', explist, 'remakeV','recalc')    forces rebuilding of FullV file,
%and re-indexing of spkblk files, for just the named experiments
%
%to build a FullV file it must first build an index of what time samples
%are available for what probes in all the .mat file made by spike2.
%jan 2012  Changed so that default saved format uses int16 rather than
%double

bysuffix = 0;
scanfiles = 0;
recalc = 0;
X.checkv = 0;
sizecheck = 0;  %set to max safe bytes for FullV to control if FullV is kept in memory
runautocut = 1;
memsize = CheckPhysicalMemory;
sizecheck = (memsize * 1e6) .* 0.8;
X.spkrate = 50;
template = [];
args = {'ndives' 0}; %don't allow chasing of trigger by default
version = 1.0;
parallel = 0;


if isdir(name)
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

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'V')
            FullV = varargin{j};
        end
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
        args = {args{:} 'nocut'};
    elseif strncmpi(varargin{j},'refcut',4)
        runautocut = 0;
        args = {args{:} 'reapply' 'refclusters'};
    elseif strncmpi(varargin{j},'parallel',6)
        parallel = 1;
    elseif strncmpi(varargin{j},'remakeV',7)
        recalc = 4;
    elseif strncmpi(varargin{j},'recalcsuffix',8)
        recalc = 5;
    elseif strncmpi(varargin{j},'recalcall',8)
        recalc = 3;
    elseif strncmpi(varargin{j},'recalc',3)
        if recalc == 4
            recalc = 5;
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
X.sizecheck = sizecheck;

ns = 1;
res.starts(2) = now;
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
            if runautocut == 2
            res.cls{newsuff(j)} = ProcessSuffix(path, probes, newsuff(j), 0, X, {args{:} 'nocut'});
            else
            res.cls{newsuff(j)} = ProcessSuffix(path, probes, newsuff(j), 0, X, args);
            end
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
    if runautocut == 2
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
return;
else

    if recalc == 3
        probes = MakeProbeIndex(path,'recalc');
    elseif length(expts) && recalc == 4 %take probe ids from exist
        probes = MakeProbeIndex(path,'suffix',expts);
    elseif length(expts) && recalc == 5  %reread spkblk
        probes = MakeProbeIndex(path,'suffixrecalc',expts);
        res.probes = probes;
    else
    probes = MakeProbeIndex(path);
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
    cls = {};
    parfor (ex = expts)
        starts(ex) = now;
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
            cls{ex} = ProcessSuffix(path, probes, ex,0, X,args);
        elseif exist(outname,'file') && recalc < 2
                fprintf('Loading %s\n',outname);
                load(outname);
        else
            fprintf('Making FullV for Expt %d\n',ex)
            ts = now;
            if bysuffix
                cls{ex} =  ProcessSuffix(path, probes, ex,0, X, args);
            else
                
            end
            ProcessSuffix(path, probes, ex,0, X,args);
        end
    end
    res.cls=cls;
    res.starts = starts;
else
    
    
    for ex = expts
        res.starts = [res.starts now];
        err = 0;
        id = find([probes.suffix] == ex);
        eprobes = unique([probes(id).probe]);
        if length(eprobes) < length(probelist)
            fprintf('Expt %d Has only %d probes\n',ex,length(eprobes));
            err = 1;
        end
        if ~exist('FullV','var')
            outname = [path '/Expt' num2str(ex) 'FullV.mat'];
            if err > 0
            elseif bysuffix
                ProcessSuffix(path, probes, ex,0, X,args);
            elseif exist(outname,'file') && recalc < 2
                fprintf('Loading %s\n',outname);
                load(outname);
            else
                fprintf('Making FullV for Expt %d\n',ex)
                ts = now;
                if bysuffix
                    res.cls{ex} =  ProcessSuffix(path, probes, ex,0, X, args);
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
                    AllVPcs(FullV,'tchan',[1:size(FullV.V,1)],'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal',args,'savespikes','autocutall');
                    clear FullV;
                else
                    ProcessSuffix(path, probes, ex,0, X,args);
                    pack;
                end
            end
        end
    end
end
if runautocut == 2  %do all teh cutting after building all the files
    res.startcuts = now;
    for ex = expts
            ProcessSuffix(path, probes, ex,0, X,args);
    end
end
res.end = now;
if logfid > 0
fclose(logfid);
end

function [res, details] = ProcessSuffix(path, probes, ex, lastblk, X, args)

details.loadtime = 0;
template = [];
nocut = 0;
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
    elseif strncmpi(args{j},'nocut',5)
        nocut = 1;
    elseif strncmpi(args{j},'noframes',5)
        exptargs = {exptargs{:} args{j}};
    elseif strncmpi(args{j},'refclusters',9)
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
if exist(outname,'file') && X.recalc < 2
    ts = now;
    fprintf('Loading %s\n',outname);
    d = dir(outname);
    clear FullV;
    load(outname);
    if isfield(FullV, 'intscale')
        FullV.V = double(FullV.V) .* FullV.intscale(1)/FullV.intscale(2);
    end
    FullV.loadname = outname;
    res.loadtime = mytoc(ts);
%    FullV.name = outname;
    if X.logfid > 0
        fprintf(X.logfid,'Load Took %.1f\r\n',res.loadtime);
    end
else
    ts = now;
    FullV = BuildFullVFile(path, probes, ex, lastblk, outname);
    res.loadtime = mytoc(ts);
    if X.logfid > 0
        if isfield(FullV,'error')
            fprintf(X.logfid,'File %s not made: error %d\r\n',outname,FullV.error);
            return;
        else
            fprintf(X.logfid,'Build %s  Took %.1f\r\n',outname,res.loadtime);
        end
    end
end
fprintf('File Load/Build tood %.2f\n',res.loadtime);
details.loadtime = res.loadtime;
if mkint && ~isfield(FullV,'intscale')
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
    Expt = LoadExpt(path, ex, X.recalc > 1, X.logfid);
    if isempty(Expt)
        fprintf('Can''t read Expt %d\n',ex);
    else
        if checkexpts == 2
            rebuild = CheckFullV(FullV,Expt,'rebuild');
            if rebuild
                FullV = BuildFullVFile(path, probes, ex, lastblk, outname);
            end
        else
            CheckFullV(FullV,Expt);
        end
    end
    return;
end
if ismember(X.recalc, [4 5]) %just make FullV
    return;
end
Expt = LoadExpt(FullV.matfile, ex, X.recalc > 1, X.logfid, exptargs);
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
if nocut
    
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
if refcut
    res = {};
    for j = probelist
        res{j} = AllVPcs(FullV,'tchan',j,'logfid',X.logfid,'savespikes','noninteractive',args{:});
    end
    clear ms;



elseif ~isempty(template)
    for j = probelist
        res{j} = AllVPcs(FullV,'tchan',j,'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal','logfid',X.logfid,args{:},'noninteractive');
    end
    clear ms;
else
    res = AllVPcs(FullV,'tchan',probelist,'nprobepc',1,'tryall','spkrate',X.spkrate,'cutmode','mahal',args{:},'savespikes','noninteractive','autocutall');
end
if lastblk
    PrintMsg(X.logfid,sprintf('%s Too big for one file',outname));
    clear FullV
    ProcessSuffix(path, probes, ex, lastblk, X, args);
end
   
function FullV = BuildFullVFile(path, probes, ex, lastblk, outname)

version = 1.0;
probelist = unique([probes.probe]);
ts = now;
id = find([probes.suffix] == ex);
matfile = [];
for j = 1:length(id)
    mfiles{j} = regexprep(probes(id(j)).file,'\.spkblk[0-9]*\.mat','.mat');
    mfiles{j} = regexprep(mfiles{j},'A.([0-9]*).mat','.$1.mat');
end
mfiles = unique(mfiles);
matfile = [path '/' mfiles{1}];
FullV = PlotSpikeC(probes(id),5,'probes',probelist,'sumv','submean','makev','lastblk',lastblk,'prefix',path);
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
fprintf('Writing %s at %s\n',outname,datestr(now));
save(outname,'FullV','-v7.3');

function rebuild = CheckFullV(FullV, Expt, varargin)
rebuild = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'rebuild',5)
        rebuild = 1;
    end
    j = j+1;
end

sds = std(FullV.V,[],2);
id = find(sds < 0.01); %Check for empty channels
if length(id)
    fprintf('Expt %d Empty Flat Channels %s\n', FullV.exptno, sprintf(' %d',id));
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
    fprintf('Missing %d Trials for %s: %s\n',length(id),Expt.Header.name,sprintf(' %d',id));
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
            if regexp(name,'lem/M[0-9]*')
                smrname = regexprep(name,'lem/M([0-9]*)','$0/lemM$1');
                smrname = regexprep(smrname,'online/lemM([0-9]*)','$0/lemM$1');
                smrname = regexprep(smrname,'/$','');
            elseif regexp(name,'dae/M[0-9]*')
                smrname = regexprep(name,'dae/M([0-9]*)','$0/daeM$1');
                smrname = regexprep(smrname,'online/daeM([0-9]*)','$0/daeM$1');
                smrname = regexprep(smrname,'/$','');
            elseif regexp(name,'ica/M[0-9]*')
                smrname = regexprep(name,'ica/M([0-9]*)','$0/icaM$1');
                smrname = regexprep(smrname,'online/iacM([0-9]*)','$0/icaM$1');
                smrname = regexprep(smrname,'/$','');
            end
            exfile = [smrname '.' num2str(ei) 'idx.mat'];
            matfile = [smrname '.' num2str(ei) '.mat'];
       elseif exist(name,'file')
                matfile = name;
                exfile = strrep(matfile,'.mat','idx.mat');
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