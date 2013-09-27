function res = RunAllVPcs(name, varargin)
%res = RunAllVPcs(dir, .....
%Reads all FullV Matlab files in a directory and calls AllVPcs
%
%res = RunAllVPcs(dirname,'quantifyall','savespikes')
%Uses existing cluster boundary an recalculates statisitcs (finishining
%"quick" saves)
%
%res = RunAllVPcs(dirname,'reclassifyall','savespikes')
%rebuilds all files using clusters currently defined on disk.
%
%RunAllVPcs(dirname,'expts',exlist, 'tchan',tlist, 'savespikes', 'refineall',Clusters)
%Uses clusters defined in Clusters and reapplies these to probes named in
%tlist. If Clusters is not provided, it loads existing Cluster Files for
%each expt
%
% Clusters{tlist(j)} is applied to probe tlist(j)
%
%chk = RunAllVPcs(res,...) reads through a returned strucutre to find
%           errors/outcomes
%cmds = RunAllCPcs(chk,'run', args..)  builds a list of commands to run to
%           fix errors in chk, including unsafetosave, needing quantifying,
%           etc. args are added to arg list for AllVPcs
%
%RunallVpcs(cmds,'run')  runs the listed commands
%RunallVpcs(cmds,'runinteractive')  runs the listed commands
expts = [];
bysuffix = 0;
scanfiles = 0;
recalc = 0;
sizecheck = 0;  %set to max safe bytes for FullV to control if FullV is kept in memory
runautocut = 1;
%memsize = CheckPhysicalMemory;
%sizecheck = (memsize * 1e6) .* 0.8;
X.spkrate = 50;
findunsafe = 0;
template = [];
forcerun = 0;
parallel = 0;
userefcluster = 0;
condenselist = 0;
runinteractive = 0;
nskip = 0;

args = {};
addargs = {};
checktype = 'reclassify';

if ischar(name)
    logfile = [name '/' 'RunAll.log'];
    logfid = fopen(logfile,'a');
elseif isfield(name,'prefix')
    logfile = [name.prefix '/' 'RunAll.log'];
    logfid = -1;
else
    logfid = -1;
end
X.logfid = logfid;


j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        args = {args{:} varargin{j}};
    elseif strncmpi(varargin{j},'args',4)
        j = j+1;
        addargs = varargin{j};
    elseif strncmpi(varargin{j},'checktype',9)
        j = j+1;
        checktype = varargin{j};
    elseif strncmpi(varargin{j},'expts',3)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'finderrs',6)
        if length(varargin) > j
            errs = FindErrorMessages(name, varargin{j+1},0);
        else
            errs = FindErrorMessages(name, 0,0);
        end
        PlotTimes(name);
        return;
    elseif strncmpi(varargin{j},'findunsafe',6)
        j = j+1;
        findunsafe = varargin{j};
    elseif strncmpi(varargin{j},'nocut',4)
        runautocut = 0;
        args = {args{:} 'nocut'};
    elseif strncmpi(varargin{j},'paralell',4)
        parallel = 1;
    elseif strncmpi(varargin{j},'refcluster',9)
        userefcluster = 1;
        args = {args{:} varargin{j}};
    elseif strncmpi(varargin{j},'run',3)
        forcerun = j+1;
        if strncmpi(varargin{j},'runinteractive',4)
            runinteractive = 1;
        end
        
    elseif strncmpi(varargin{j},'skip',4)
        j = j+1;
        nskip = varargin{j};
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
res.starts(1) = now;
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:MaxIterations');

if iscellstr(name) && isdir(name{1}) %given a list of directories
   for j = 1:length(name)
       res.dres{j} = RunAllVPcs(name{j},varargin{:});
   end
   return;
end

if iscell(name) &&  length(name) > 1
     name = CellToStruct(name);
end

     
if isstruct(name)
    if isfield(name,'name') && isfield(name,'probes') && isfield(name,'args');
        if isfield(name,'overlapn') && forcerun == 0
            res = CheckReclassify(name,'reclassify');
        else
            res = RunCommandList(name, varargin{:});
        end
    elseif isfield(name,'cmds') && forcerun > 0 %actually run a command list
        res = {};
        for j = 1:length(name)
            res{j}.starts = now;
            name(j).interactive = runinteractive;
            res{j} = RunSelectedCommands(name(j), nskip, varargin{forcerun:end});
        end
    elseif isfield(name,'acts') && forcerun > 0
        res = {};
        for j = 1:length(name)
            res{j} = RunSelectedCommands(name(j), nskip, varargin{forcerun:end});
        end
        if length(res) > 1
            a.cmds = {};
            a.args = {};
            a.reason = [];
            a.expts = [];
            a.probes = [];
            a.files = {};
        for j = 1:length(res)
            if ~isempty(res{j})
            a.args = {a.args{:} res{j}.args{:}};
            a.cmds = {a.cmds{:} res{j}.cmds{:}};
            a.files = {a.files{:} res{j}.files{:}};
            a.reason = [a.reason res{j}.reason];
            a.expts = [a.expts res{j}.expts];
            a.probes = [a.probes res{j}.probes];
            end
        end
        res = a;
        end
        
    elseif isfield(name,'overlapn') && forcerun == 0
            res = CheckReclassify(name,'reclassify');
    elseif isfield(name,'acts')
        if ~isfield(name.acts,'name')
            for j = 1:length(name.acts)
                name.acts(j).name = sprintf('%s/%s',name.prefix,name.name{name.acts(j).exptid});
            end
        else
            for j = 1:length(name.acts)
                if strncmp(name.acts(j).name,'Expt',4)  %no directory prefix
                    name.acts(j).name = sprintf('%s/%s',name.prefix,name.acts(j).name);
                end
            end
        end
        res = RunCommandList(name.acts,'log',logfile,varargin{:});
    elseif isfield(name,'dres') && forcerun
        res ={};
        for j = 1:length(name.dres)
            res{j} = RunSelectedCommands(name.dres{j}, nskip, varargin{forcerun:end});
        end
        for j = 1:length(res)
                good(j) = length(res{j}.acts);
               for k = 1:length(res{j}.acts)
                   [a,b] = fileparts(res{j}.prefix);
                   fprintf('%s E%dP%d\n',b,res{j}.acts(k).exptno,res{j}.acts(k).probes(1));
               end
        end
        res = res(find(good));
    elseif isfield(name,'dres')
        res ={};
        for j = 1:length(name.dres)
            fprintf('Dir%d %s\n',j,name.dres{j}.prefix);
            res{j} = CheckReclassify(name.dres{j},checktype);
            if findunsafe > 0
                id = FineUnsafes(res{j}, findunsafe);
                res{j}.acts = res{j}.acts(id);
            else
                if condenselist
                    res{j} = rmfields(res{j},'acts');
                end
                if isfield(name.dres{j},'exptlist')
                    res{j}.exptlist = name.dres{j}.exptlist;
                end
            end
        end
        if findunsafe > 0
            for j = 1:length(res)
                good(j) = length(res{j}.acts);
               for k = 1:length(res{j}.acts)
                   [a,b] = fileparts(res{j}.prefix);
                   fprintf('%s E%dP%d\n',b,res{j}.acts(k).exptno,res{j}.acts(k).probes(1));
               end
            end
            res = res(find(good));
        end
    elseif isfield(name,'cls')
        res = CheckReclassify(name,checktype);
        if findunsafe > 0
            id = FineUnsafes(res, findunsafe);
            res.acts = res.acts(id); 
        elseif condenselist
            res = rmfields(res{j},'acts');
        end
    elseif isfield(name,'exptlist')
        if isfield(name,'args') &&length(addargs)
            name.args = {name.args{:} addargs{:}};
        end
        name = CheckExptList(name);
        res = CheckReclassify(name, checktype);
    end
    return;
 elseif isdir(name)
    path = name;
else
[path, fname] = fileparts(name);
end

d = dir(name);
ns = 1;
res.prefix = name;
for j = 1:length(d)
        sid = regexp(d(j).name,'Expt[0-9,a]*FullV.mat');
        X.logfid = logfid;
        if length(sid)
            ex = sscanf(d(j).name(sid(1)+4:end),'%d');
            if isempty(expts) || ismember(ex,expts)
                if logfid > 0
                fprintf(logfid,'Processing %s %s\n',d(j).name,datestr(now));
                end
                name = [path '/' d(j).name];
                if parallel == 0
                    res.cls{ns} =  ProcessSuffix(name, ex, args{:});
                else
                    res.cls{ns} = [];
                end
 %              fclose all;
               res.exptlist(ns) = GetExptNumber(d(j).name);
               res.name{ns} = d(j).name;
               filenames{ns} = name;
               ns = ns+1;
            end
        end
        
end

if ns <= 1
    fprintf('No FullV Files in %s\n',name);
end
if parallel && ns > 1
    args = {args{:} 'nowatch' 'nocheck' 'noninteractive' 'verbose'};
    workers = [];
    parfor  (j = 1:ns-1)
        t = mygetCurrentTask();
        fprintf('Expt%d workder id is %d\n',res.exptlist(j),t.ID);
        cls{j} =  ProcessSuffix(filenames{j}, res.exptlist(j), args{:});
        workers(j) = t.ID;
    end
    res.cls = cls;
    res.workerid = workers;
end

 res.end = now;
 res.args = args;
 if logfid > 0
     fclose(logfid);
 end

function need = FineUnsafes(res, findunsafe)
    need = [];
    if isempty(res)
        return;
    end
    need = zeros(size(res.acts));
    for j = 1:length(res.acts)
        if isfield(res.acts,'unsafetosave') && res.acts(j).unsafetosave > 0
        if bitand(res.acts(j).unsafetosave,findunsafe) > 0
            need(j) = 1;
        end
        end
    end
    need = find(need);
    
function res = ProcessSuffix(path, ex, varargin)
template = [];
nocut = 0;
res = {};
j = 1;
checkexpts = 0;
makesmall = 1;
quantify = 0;
checkspikes = 0;
userefcluster = 0;

args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'checkspikes',10)
        checkspikes = 1;
    elseif strncmpi(varargin{j},'check',5)
        checkexpts = 1;
        if strncmpi(varargin{j},'checkandbuild',8)
            checkexpts = 2;
        end
    elseif strncmpi(varargin{j},'nocut',5)
        nocut = 1;
    elseif strncmpi(varargin{j},'quantify',5)
        quantify = 1;
        args = {args{:} varargin{j}};
    elseif strncmpi(varargin{j},'refcluster',9)
        userefcluster = 1;
        args = {args{:} varargin{j}};
    elseif strncmpi(varargin{j},'template',6)
        j = j+1;
        template = varargin{j};
        cargn = j;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

nid = [];
res.filename = path;
outname = path;
if exist(outname,'file')
    ts = now;
    fprintf('%s Exists\n',outname);
    d = dir(outname);
else
    fprintf('Cant Read %s\n',outname);
    return;
end

if quantify || checkspikes
    cfile = strrep(outname,'FullV.mat','ClusterTimes.mat');
    if ~exist(cfile)
        if userefcluster
            fprintf('Using RefClussters for all Probes\n');
            nid = 1:24;
        else
            fprintf('No Cluster File %s\n',cfile);
            return;
        end
    else
        load(cfile);
        for j = 1:length(Clusters)
            if (isfield(Clusters{j},'quick') && Clusters{j}.quick > 0)
                needc(j) = 1;
            elseif userefcluster && (isfield(Clusters{j},'auto') && Clusters{j}.auto > 0)
                needc(j) = 5;
            elseif (isfield(Clusters{j},'manual') && Clusters{j}.manual ==2) 
                needc(j) = 4;
            elseif Clusters{j}.dropi(3) == 0
                needc(j) = 3;                
            elseif CheckSpikeFile(Clusters{j}, outname) == 0
                needc(j) = 2;
            else
                needc(j) =0;
            end
            for k = 1:length(Clusters{j}.next)
                if isfield(Clusters{j}.next{k},'manual') && Clusters{j}.next{k}.manual == 2
                    needc(j) = 4;
                elseif isfield(Clusters{j}.next{k},'quick') && Clusters{j}.next{k}.quick > 0
                    needc(j) = 1;
                end
            end
        end
        if sum(needc) == 0
            fprintf('No Clusters need quantifying in %s\n',cfile);
            return;
        else
            nid = find(needc > 0);
            fprintf('Need %s (%s) in %s\n',sprintf('%d ',nid),sprintf('%d ',needc(nid)),cfile);
        end
    end
end

if checkspikes
    id = find(needc);
    if isempty(args)
        fprintf('Need Probes %s\n',sprintf(' %d',id));
    else
    for j = 1:length(id)
        res.cls{j} = AllVPcs(outname,'tchan',id(j),args{:});
    end
    end
    return;
end


ts = now;
try
res = AllVPcs(outname,args{:});
for j = 1:length(nid)
    res{nid(j)}.needed = 1;
end

if ~iscell(res)
    res.runtime = mytoc(ts);
else

    toplevel = 0;
    for j = 1:length(res)
        if isfield(res{j},'toplevel')
            toplevel = res{j}.toplevel;
        end
        if makesmall
            res{j} = rmfields(res{j},'t','Evec','pcs');
            if isfield(res{j},'cluster')
            res{j}.cluster = rmfields(res{j}.cluster,'times');
            end
        end
    end
    if toplevel > 0 && isappdata(toplevel,'Vall')
        rmappdata(toplevel,'Vall');
    end
end
catch ME
    fprintf('!!!ALLVPCS FAILED %s for file %s (line %d, m-file %s)\n',ME.message,outname,ME.stack(1).line,ME.stack(1).name);
    res.exception = ME;
    res.time = now; 
    t = mygetCurrentTask();
    res.worker = t.ID;
end


function ok = CheckSpikeFile(C, fullvname)
    ok = 1;
    monkey = GetMonkeyName(fullvname);
    id = findstr(fullvname,monkey);
    if ~isfield(C,'spkfile')
        fprintf('No spike file in Cluster %d\n',C.probe(1));
        return;
    end
    bid = findstr(C.spkfile,monkey);
    spkfile = [fullvname(1:id(1)) C.spkfile(1+bid(1):end)];
    spkfile = strrep(spkfile,['/Spikes/' monkey],'/Spikes/');

    ddir = fileparts(fullvname);
    [a,fname] = fileparts(C.spkfile);
    spkfile = [ddir '/Spikes/' fname '.mat'];
    newspkfile = strrep(spkfile,'/Spikes/',['/Spikes/' monkey]);
    if ~exist(spkfile) && ~exist(newspkfile)
        fprintf('Missing Spk file %s and %s\n',spkfile,newspkfile);
        ok = 0;
    else
        if exist(newspkfile)
            spkfile = newspkfile;
        end
        d = dir(spkfile);
        if d.datenum < C.savetime(1)-0.01
            ok = 0;
            fprintf('Spk file %s is older (%s) than Cluster (%s)\n',spkfile,d.date,datestr(C.savetime(1)));
        elseif d.datenum > C.savetime(1)+1
            ok = 0;
            fprintf('Spk file %s is much newer (%s) than Cluster (%s)\n',spkfile,d.date,datestr(C.savetime(1)));
        end
    end
        
    
function Ex = LoadExpt(name, ei, rebuild, logfid)
        Ex = [];
        smrname = regexprep(name,'lem/M([0-9]*)','$0/lemM$1');
        exfile = [smrname '.' num2str(ei) 'idx.mat'];
        matfile = [smrname '.' num2str(ei) '.mat'];
        if exist(matfile) && (~exist(exfile,'file') || rebuild)
            PrintMsg(logfid,'Building %s\n',exfile);
            APlaySpkFile(matfile,'bysuffix','noerrs');
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
        end
        
function res = RunSelectedCommands(X, skip, varargin)
    res = [];
    nc = 0;
    docells = [0 1];
    parallel = 0;
    reasons{1} = 'Error';
    reasons{3} = '';
    reasons{8} = 'Unsafe <0 (auto cut)';
    reasons{9} = '';
    reasons{7} = 'Unsafe > 1';
    reasons{10} = 'Unsafe but >1 trigger chan';
    reasons{128} = 'Spk File Old';
    allvargs = {};
    j = 1;
    while j <= length(varargin)
        if strncmp(varargin{j},'cellsonly',4)
            docells = 1;
        elseif strncmp(varargin{j},'nocells',5)
            docells = 0;
        elseif strncmp(varargin{j},'parallel',5)
            parallel = 1;
        else
            allvargs = {allvargs{:} varargin{j}};
        end
        j = j+1;
    end
    if isfield(X,'acts') && isfield(X.acts,'exptno')
        cellim = X.cellimage > 0;
        for k = 1+skip:length(X.acts)
            vname = [X.prefix '/Expt' num2str(X.acts(k).exptno) 'FullV.mat'];
            p = X.acts(k).probes;
            e = X.acts(k).exptno;
%actimage = 1 + CheckResult()
% 1 = Error,  3 = auto = 2, 8 = unsafetosave < 0, 7 = unsafetosave > 0
% 10 = unsafe > 0 but trigger > 1 channel, so check
% 11 = Emtpy 'next' fields with non-empty ones later. Usuall can ignore
% 128 = needed becuase of dates/files
            if ismember(X.actimage(e,p),[1 3 7 8 10 128]) & ismember(cellim(e,p),docells)
                nc = nc+1;
                res.reason(nc) = X.actimage(e,p);
                if X.cellimage(e,p) > 0
                    cellchr = '* ';
                else
                    cellchr = ' ';
                end
                cmd = ['AllVPcs(' vname, ',''tchan'',' num2str(p) ',' sprintf('%s, ', allvargs{:}) ') %%' num2str(res.reason(nc)) cellchr reasons{res.reason(nc)} ];
                fprintf([cmd '\n']);
                res.cmds{nc} = cmd;
                res.args{nc} = {vname 'tchan' p allvargs{:}};
                res.expts(nc) = e;
                res.probes(nc) = p;
                [a,res.files{nc}] = fileparts(X.prefix);
            end
        end
        if ~isempty(res)
        [a,b] = sort(res.expts);
        res.cmds = res.cmds(b);
        res.args = res.args(b);
        res.reason = res.reason(b);
        res.expts = res.expts(b);
        res.probes = res.probes(b);
        res.files = res.files(b);
        end
    elseif isfield(X,'cmds')
        for j = 1+skip:length(X.files)
            files{j-skip} = [X.files{j} 'Ex' num2str(X.expts(j))];
        end
        ufiles = unique(files);
        res = {};
        if parallel && length(unique(ufiles)) > 1
            parfor j = (1:length(ufiles));
                id = skip+find(strcmp(ufiles{j},files));
                for k = 1:length(id)
                    t = mygetCurrentTask();
                    fprintf('Worker %d Running %s\n',t.ID,X.cmds{id(k)});
                    cls{j}{k} = AllVPcs(X.args{id(k)}{:},allvargs{:},'noninteractive');
%                    cls{j}{k} = X.args{id(k)};
                end
            end
            for j = (1:length(ufiles));
                id = find(strcmp(ufiles{j},files));
                for k = 1:length(id)
                    res{id(k)} = cls{j}{k};
                end
            end
            return;   
        end
        for j = 1+skip:length(X.args)
            res{j} = AllVPcs(X.args{j}{:},allvargs{:});
            res{j}.saved = 0;
            if X.interactive
            options.WindowStyle = 'normal';
            options.Interpreter = 'none';
            options.Resize = 'on';
            options.Default = 'No';
            if X.reason == 9
                needed = 1;
            else
                needed = 0;
            end
            if isfield(res{j},'unsafetosave')
                if isfield(res{j},'needed') && needed == 0;
                    needed = res{j}.needed;
                end
                msg = sprintf('Unsafe %d Trig %d/%d unmatch %d/%d Needed %d',...
                    res{j}.unsafetosave,sum(res{j}.trigmatch(1:3)),res{j}.trigmatch(4),...
                    sum(res{j}.matchcounts([1 3])),res{j}.matchcounts(4),needed);
            else
                msg = 'Hit OK to save. Modify this to cancel and record reason';
            end
                
            res{j}.saved = 0;
            yn =  mydialog('Save Spikes? (Change message to cancel+comment)','Call AllVPCs',...
                'text',msg,'parent',res{j}.toplevel,...
                'extrabutton','Stop');
            drawnow; %get rid of dialog
            if ~isempty(yn) %hit OK
                if strcmp(msg,yn) %message unchanged
                    AllVPcs(res{j}.toplevel,'savespikes');
                    res{j}.saved = 1;
                elseif strcmp(yn,'Stop')
                    return;
                else
                    res{j}.comment = yn;
                end
            end
            end
        end
    end

        
    function res = RunCommandList(X, varargin)
     
        checkmode = 0;
        logfile = [];
        addargs = {};
        j = 1;
        if ~isfield(X,'probes') && isfield(X,'acts')
        else
            probes = unique([X.probes]);
            expts = unique([X.exptid]);
        end
        
        while j <= length(varargin)
            if strncmpi(varargin{j},'check',5);
                checkmode = 1;
            elseif strncmpi(varargin{j},'args',4);
                j = j+1;
                addargs = varargin{j};
            elseif strncmpi(varargin{j},'expts',5);
                j = j+1;
                expts = varargin{j};
            elseif strncmpi(varargin{j},'probes',5);
                j = j+1;
                probes = varargin{j};
            elseif strncmpi(varargin{j},'log',3);
                j = j+1;
                logfile = varargin{j};
            end
            j = j+1;
        end

        rmid = zeros(size(X));
        for j = 1:length(X)
            allnames{j} = X(j).name;
            if isempty(X(j).args)
                rmid(j) = 1;
            end
        end
        id = setdiff(1:length(X),find(rmid));
        X = X(id);
        allnames = allnames(id);
        id = find(ismember([X.exptid],expts) & ismember([X.probes],probes));
        X = X(id);
        allnames = allnames(id);
        
    names = unique(allnames);
    n = 1;
    for k = 1:length(names)
        id = strmatch(names{k},allnames);
        name = names{k};
        ex = GetExptNumber(name);
        fprintf('Loading %s\n',name);
            if checkmode == 0 && sum([X(id).probes]) > 0
                load(name);
            end
        for j = id'
            ex = X(j).exptid;
            for k = 1:length(X(j).probes)
                cmd = sprintf(['AllVPcs(FullV, ''tchan'',' num2str(X(j).probes) ',' sprintf('%s, ', X(j).args{:}) ')']);
                fprintf([cmd '\n']);
                if checkmode == 0
                    if length(logfile)
                        logfid = fopen('logfile','a');
                        if logfid > 0
                            fprintf(logfid,'%s %s\n',datestr(now),cmd);
                        end
                    end
                    res{n} = AllVPcs(FullV, 'tchan', X(j).probes(k), X(j).args{:}, addargs{:});
                    fclose all;
                    res{n}.name = name;
                    res{n}.args = X(j).args;
                    res{n}.probes = X(j).probes(k);
                else
                    res.cmds{n} =  cmd;
                    if ex > 0
                    res.map(ex, X(j).probes(k)) = 1;
                    end
                end
                n = n+1;
            end
        end
        if checkmode == 0 & isfield(res{n-1},'toplevel') %AllVPCs was called
            close(res{n-1}.toplevel);
            clear FullV;
        end
    end


 function X = CheckExptList(X)
             %if args starts by listing all probes, strip this - want to call
        %one at at time next
        if strcmp(X.args{1},'tchan') && length(X.args{2}) > 2
            X.args = {X.args{3:end}};
        end
        if strcmp(X.args{1},'reclassifyall') 
            X.args{1} = 'reclassify';
        end

     for j = 1:size(X.cls)
         exptno = [];
         for k = 1:size(X.cls{j})
             exptno(k) = X.cls{j}{k}.cluster.exptno;
         end
         X.exptlist(j) = median(exptno);
     end
            
function result = CheckReclassify(X, checktype, varargin)
    
    exptno = 0;
    j = 1;
    result.allerrs = {};
    while j <= length(varargin)
        if strncmp(varargin{j},'exptno',6)
            j = j+1;
            exptno = varargin{j};
        end
        j = j+1;
    end
    if isfield(X,'cls')
        result.msg = {};
        result.msgp = [];
        result.msge = [];
        result.auto = [];
    [exps, exi] = sort(X.exptlist);
        acts = [];
        startid = 0;
        for j = 1:length(X.cls)
            e = X.exptlist(j);
            %e is the correct row for the results structre
            e = find(exps == X.exptlist(j));
            if iscell(X.cls{j})
                for k = 1:length(X.cls{j})
                    C = X.cls{j}{k};
                    if ~isfield(C,'spkrate')
                        C.spkrate = NaN;
                    end
                    if ~isfield(X.cls{j}{k},'probes')
                        X.cls{j}{k}.probes = k;
                    end
                    if ~isfield(X.cls{j}{k},'exptno')
                        X.cls{j}{k}.exptno = j;
                    end
                    if isfield(X.cls{j}{k},'errs')
                        a = length(result.allerrs)+1;
                        result.allerrs = {result.allerrs{:} X.cls{j}{k}.errs{:}};
                        b = length(result.allerrs);
                        result.errid(1,a:b) = X.cls{j}{k}.probes(1);
                        result.errid(2,a:b) = X.exptlist(j); %expt no, not row
                    end
                    X.cls{j}{k}.exptid = e;
                    if ~isfield(X.cls{j}{k},'args')
                        X.cls{j}{k}.args = X.args;
                    end                        
                end
                
                fprintf('%sExpt%d:\n',X.prefix,X.exptlist(j));
                res = CheckReclassify(X.cls{j},checktype);
                e = X.exptlist(j);
                if isfield(res.acts,'errstate')
                    for k = 1:length(res.acts)
                        iserr(k) = ~isempty(res.acts(k).errstate);
                    end
                    result.errstate = [res.acts.errstate];
                    result.errstateid = find(iserr);
                    result.errstatex = ones(size(result.errstateid)) .*e;
                end
                np = length(res.ctime);
                result.msgp = [result.msgp res.msgid];
                result.msge = [result.msge ones(size(res.msgid)).*X.exptlist(j)];
                result.actimage(e,:) = res.actimage(end,:);
                result.cellimage(e,:) = res.cellimage(end,:);
                if isfield(res.acts,'unsafetosave')
                    for u = 1:length(res.acts)
                        if isempty(res.acts(u).unsafetosave)
                            result.unsafe(e,u) = NaN;
                        else
                            result.unsafe(e,u) = res.acts(u).unsafetosave;
                        end
                    end
                end
                result.dateimage(e,:) = res.ctime;
                result.memimage(e,:) = res.memsz(:,1);
                result.auto(e,1:length(res.auto)) = res.auto;
                result.xclim(e,1:length(res.xcl)) = res.xcl;
                result.name{e} = X.name{j};
                result.msg = {result.msg{:} res.msg{:}};
                startid = startid + length(X.cls{j});
                if isfield(res,'readmethod')
                    result.readmethod(e,1:length(res.readmethod)) = res.readmethod;
                end
                try
                    if ~isempty(acts)
                        a = fields(acts);
                        b = fields(res.acts);
                        c = setdiff(a,b);
                        for k = 1:length(c)
                            [res.acts.(c{k})] = deal(NaN);
                        end
                        c = setdiff(b,a);
                        for k = 1:length(c)
                            [acts.(c{k})] = deal(NaN);
                        end
                    end
                    acts = [acts res.acts];
                catch
                    fprintf('Error catting acts E%d\n',e);
                end
            elseif isfield(X.cls{j},'exception')
                t = mygetCurrentTask();
                s = sprintf('!!!Error (Exception) %s in %s line %d at %s Worker %d',X.cls{j}.exception.message,X.cls{j}.exception.stack(1).file,X.cls{j}.exception.stack(1).line,datestr(X.cls{j}.time),t.ID);
                result.msg = {result.msg{:} s};
                result.msgp = [result.msgp 0];
                result.msge = [result.msge X.exptlist(j)];
                fprintf('%s\n',s);
            elseif isfield(X.cls{j},'filename')
                s = sprintf('%s empty',X.cls{j}.filename);                
            else
                s = sprintf('!!!Error Unrecognized');                
                fprintf('%s\n',s);
            end
        end
        result.acts = acts;
        if isfield(X,'prefix')
            result.prefix = X.prefix;
        end
        return;
    end
    if iscell(X)
        res = CellToStruct(X);
    else
        f = {'probes' 'exptno'};
        for j = 1:length(X)
            for k = 1:length(f);
                if isfield(X,f{k});
                    res(j).(f{k}) = X(j).(f{k});
                end
            end
        end
    end
    msgs = {};
    msgid = [];
    actimage = [];
    cellimage = [];
    for j = 1:length(X) 
        if iscell(X)
            res(j).cell = NaN;
            [a, msg] = CheckResult(X{j}, checktype);
            if isfield(a,'errstate')
                res(j).result = -1;
                res(j).errstate = a.errstate;
            else
                res(j).result = a;
            end
            if ~isempty(msg)
                msgs = {msgs{:} msg};
                msgid = [msgid j];
            end
            if isfield(X{j},'cluster')
                if isfield(X{j},'savetime')
                    res(j).ctime = X{j}.cluster.savetime(1);
                else
                    res(j).ctime = X{j}.cluster.ctime;
                end
                res(j).auto = X{j}.cluster.auto;
                if isfield(X{j}.cluster,'excludetrialids')  && length(X{j}.cluster.excludetrialids) > 1
                    res(j).xcl = 1;
                else
                    res(j).xcl = 0;
                end
                res(j).memsz = X{j}.memsz;
                if isfield(X{j}.cluster,'exptreadmethod')
                    res(j).readmethod = X{j}.cluster.exptreadmethod;
                end
                res(j).auto = X{j}.cluster.auto;
                if isfield(X{j},'cell')
                    res(j).cell = X{j}.cell;
                end
            else
                res(j).memsz = [0 0];
                res(j).ctime = 0;
                res(j).auto = 0;
                res(j).xcl = 0;
            end
        elseif isfield(X(j),'cluster')
            [res(j).result, msg] = CheckResult(X(j),checktype);
            if ~isempty(msg)
                msgs = {msgs{:} msg};
                msgid = [msgid j];
            end
        else
                res(j).memsz = [0 0];
            res(j).result = 0;
            res(j).auto = 0;
            res(j).xcl = 0;
        end
        
        if isfield(res,'args')
            if res(j).result == 1 %autocut
                res(j).args = {res(j).args{:} 'autocut' 'savespikes'};
            elseif res(j).result == 2
                res(j).args = {res(j).args{:} 'savespikes'};
            elseif res(j).result == 3
                res(j).args = {res(j).args{:} 'savespikes'};
            elseif res(j).result == 4
                res(j).args = {res(j).args{:} 'savespikes'};
            else
                res(j).args = {};
            end
        end
        if isfield(res(j),'name')
            res(j).exptno= GetExptNumber(res(j).name);
        end
        if isfield(res(j),'exptno')
            if ~isfield(res(j),'unsafetosave') || isempty(res(j).unsafetosave) %did nothing
                actimage(res(j).exptno,res(j).probes) = 0;
            else
                actimage(res(j).exptno,res(j).probes) = res(j).result+1;
            end
            cellimage(res(j).exptno,res(j).probes) = res(j).cell;
            if actimage(res(j).exptno,res(j).probes) == 0 && isfield(res(j),'needed') && sum(res(j).needed)
                actimage(res(j).exptno,res(j).probes) = 128;
            end
        end
    end
    result.acts = rmfields(res,'cluster');
    result.msg = msgs;
    result.msgid = msgid;
    result.actimage = actimage;
    result.cellimage = cellimage;
    if isfield(res,'ctime')
    result.ctime = [res.ctime];
    result.auto = [res.auto];
    result.xcl = [res.xcl];
    result.memsz = cat(1,res.memsz);
    end
    if isfield(res,'readmethod')
    result.readmethod = [res.readmethod];
    end
    
function [res, s] = CheckResult(a, type)
%[res, s] = CheckResult(a, type)  Checks outcome of call to AllVPcs
% res = 6 unsafetosave > 0
% res = 7 unsafetosave < 0 (-1000 = needs quantifying)
% res = 0 error recorded
% res = 4 recluster a good match
% res = 3 recluster a good match, auto cut
% res = 2 recluster a good match, auto = 2 = Plotclusters 
% res = 1  automatic cut, bad match or error
% res = 8 automatic cut, OK. 
%res = 9 = unsafe, but also used > triggerchan - check if works wiht just
%                    one
%res = 10 empty .next{j} in cluest, j < length(next)
%res = 5 no spikes either before or after

    s = '';
    if ~isfield(a,'auto') 
        if isfield(a,'err') && ischar(a.err) %terminated by exception
            s = ['Exception: ' a.err];
            res.errstate = a.errstate;
        else
            res = 0;
            s = sprintf('E%dP%d no auto',a.exptno,a.probes(1));
        end
        fprintf('%s\n',s);
        return;
    end
    if a.auto == 2 %old method for recording plotcluster cuts.
        a.auto = 2;
    end
      if strcmp(type,'reclassify')
          if a.err >0 && a.auto == 1
              res = 1;
          elseif a.unsafetosave ~= 0
              s = sprintf('E%dP%d Unsafe %d Trig %s space %s xc %.2f %.2f overlap %d(%s)',...
                  a.cluster.exptno,a.cluster.probe(1),...
                  a.unsafetosave,sprintf(' %d',a.trigmatch),sprintf(' %d',a.cluster.space),...
                  a.xcorr(1),a.xcorr(2),...
                  a.overlapn,sprintf(' %d',a.matchcounts));
              if a.unsafetosave < 0
                  res = 7;
              elseif a.unsafetosave == 128
                  res = 10;
              elseif bitand(a.unsafetosave,64) > 0
                  res = 9;
              else
                  res =6;
              end
          elseif isfield(a,'needed') && a.needed > 0
              s = sprintf('E%dP%d Needed %d',a.cluster.exptno,a.cluster.probe(1),a.needed);
              res = 8;
          elseif a.err > 0
              fprintf('E%dP%d  err, not auto\n',a.cluster.exptno,a.cluster.probe(1)) 
              res = 0;
          elseif abs(a.xcorr(1)) > 0.95 && abs(a.xcorr(2)) > 0.95 && a.matchcounts(1)./a.matchcounts(2) < 0.15 && a.overlapn./a.matchcounts(2) > 0.5
              res = 4-a.auto; %4 or 3 or 2. All = good result
          elseif a.auto ==1
              res = 8;
          elseif a.matchcounts(1) ==0 && a.matchcounts(3) == 0
              res = 5;
          elseif length(a.matchcounts) > 5 && sum(a.matchcounts(6:end)) == 0
              res = 4;
          else
              s = sprintf('E%dP%d Unsafe %d Trig %s space %s xc %.2f %.2f overlap %d(%s)',...
                  a.cluster.exptno,a.cluster.probe(1),...
                  a.unsafetosave,sprintf(' %d',a.trigmatch),sprintf(' %d',a.cluster.space),...
                  a.xcorr(1),a.xcorr(2),...
                  a.overlapn,sprintf(' %d',a.matchcounts));
              res = 0;
          end
      elseif strcmp(type,'baddim')
          if a.err ==1 && a.auto == 1
              res = 1;
          elseif a.err ==1 && a.auto == 0
              res = -1;
          else
              res = 0;
          end
      end
      if ~isempty(s)
          fprintf('%s\n',s);
      end
      
    function CheckErrorList(errs)
        unsafe = 0;
        missing = 0;
        gmfits = 0;
        DataClusters = 0;
        suspicious = 0;
        autos = 0;
        for j = 1:length(errs)
            if strfind(errs{j},'Missing')
                missing(j) = 1;
            elseif strfind(errs{j},'Unsafe')
                unsafe(j) = 1;
            elseif strfind(errs{j},'DataClusters')
                DataClusters(j) = 1;
            elseif strfind(errs{j},'suspicous')
                suspicious(j) = 1;
            elseif regexp(errs{j},'[Gg][Mm][Ff]it')
                gmfits(j) = 1;
            elseif regexp(errs{j},'Error Fitting')
                gmfits(j) = 1;
            elseif regexp(errs{j},'no auto')
                autos(j) = 1;
            else
                unknown(j) = 1;
            end
        end
        total = sum(missing)+sum(unsafe)+sum(DataClusters)+sum(gmfits)+sum(suspicious)+sum(autos);
        fprintf('%d/%d Errors: %d Missing, %d unsafe %d DataClusters %d gmfits %d suspicious %d auto errot\n',...
            length(errs),total,sum(missing),sum(unsafe),sum(DataClusters),....
            sum(gmfits),sum(suspicious),sum(autos));
    
 function PlotTimes(S)
 %summarize time used
 if iscell(S)
     starts = [];
     memsz = [];
     for j = 1:length(S)
         if isfield(S{j},'dateimage')
         id = find(S{j}.dateimage > 0);
         starts = [starts S{j}.dateimage(id)'];
         memsz = [memsz S{j}.memimage(id)'];
         end
     end
     starts = starts .* 60 * 24; %minutes
     plot(starts-min(starts),memsz,'o');
 end
        
    function [errs, allerrs] = FindErrorMessages(S, msg, offset)%S is a summmary made by calling RunAllVPcs with one of its own results
        
        errs = [];
        allerrs = {};
       if iscell(S)
           for j = 1:length(S)
               [a, newerrs] = FindErrorMessages(S{j},msg, j);
               errs = cat(2,errs,a);
               allerrs = {allerrs{:} newerrs{:}};
           end
           CheckErrorList(allerrs);
           return;
       end
       
       if isfield(S,'prefix')
       [a,prefix] = fileparts(S.prefix);
       end
       if strcmpi(msg,'unsafe+')
           msg = 'Unsafexx';
           if isfield(S,'unsafe')
               [a,b] = find(S.unsafe > 0);
               [d,name] = fileparts(S.prefix);
               for j = 1:length(a)
                   fprintf('%s (%d) E%dP%d Unsafe %d\n',name,offset,a(j),b(j),S.unsafe(a(j),b(j)));
               end
           end
       end
       id = find(strncmp(msg,S.allerrs,length(msg)));
       ne = 0;
       e = 0;
       p = 0;
       allerrs = {allerrs{:} S.allerrs{:}};
       allerrs = {allerrs{:} S.msg{:}};       
       for j = 1:length(S.allerrs)
           if ~isempty(regexp(S.allerrs{j},msg))
               p = S.errid(1,j);
               e = S.errid(2,j);
               ne = ne+1;
               errs(ne).s = S.allerrs{j};
               fprintf('%s(%d) E%dP%d %s\n',prefix,offset,e,p,S.allerrs{j});
           end
       end
       for j = 1:length(S.msg)
           e = S.msge(j);
           p = S.msgp(j);
           if ~isempty(regexp(S.msg{j},msg))
               ne = ne+1;
               errs(ne).s = S.msg{j};
               fprintf('%s(%d) E%dP%d %s\n',prefix,offset,e,p,S.msg{j});
           end
       end
       if isfield(S,'errstate')
           for j = 1:length(S.errstate)
               allerrs = {allerrs{:} S.errstate(j).message};
               e = S.errstatex(j);
               p = S.errstateid(j);
               fprintf('%s(%d) E%dP%d %s\n',prefix,offset,e,p,S.errstate(j).message);
           end
       end
        