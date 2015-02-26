function DATA = FixClusterFiles(name, varargin)
%FixClusterFile(dir, ......) Read Clusters, check, modify if necessary
%reads in Cluster files and fixes missing fields/lists errors
%FixClusterFiles(dir, 'checkexpt','Expts',Expts) Checks that numbering of Expts
%FixClusterFiles(dir, 'checkmark') Checks that clusers marked as bad probe
%                               matches ArrrayConfig
%FixClusterFiles(dir, 'errors') Builds Error Listing
%FixClusterFiles(dir, 'details') Checks for empty ClusterDetails Structs
%FixClusterFiles(dir, 'details','write','savedetails') sasves the fixed ClusterDetails Structs
%and replaces them with contents of the latest backup
%FixClusterFiles(dir, 'expts', eid) Just checks expts in eid
%FixClusterFiles(name, 'write')  writes out any modified Cluster structs,
%       provided the modification was caused by a check named on the
%       command line. 
%in the ClusterTimes Files matches Expts
%
cmpdrive = [];
forcewrite = 0;
TAGTOP = 'PlotClusters';
args = {};
Expts = [];
errs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        TOPTAG = varargin{j};
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
initcall =0 ;
if ishandle(name)  % callback
    gui = name;
    DATA = GetDataFromFig(name);
    name = varargin{2};
    if length(varargin) > 2
    varargin = varargin{3:end};
    else
        varargin = {};
    end
end


if iscellstr(name)
    for j = 1:length(name)
        DATA{j} = FixClusterFile(name{j}, varargin{:});
    end
    return;
elseif isstruct(name) || iscell(name)
    DATA.name = 'ClusterArg';
    DATA.strings{1} = 'ClusterArg';
    if isstruct(name) && isfield(name,'cls') %template fits
        Clusters = ReadTemplateClusters(name.cls);
        DATA.Templates = name.Templates;
        DATA.datatype = 2;
    elseif isstruct(name) && isfield(name,'Clusters')
        Clusters = name.Clusters;
        DATA = FixClusters(Clusters,args{:});
        return;
    elseif isstruct(name) %one cluster result;
        Clusters{1} = name.Clusters;
    elseif isfield(name{1},'mahal')
        Clusters{1} = name;
    elseif initcall == 1
        Clusters = name;
    end
    if initcall
    DATA = InitInterface(DATA); 
    [DATA, Clusters] = ReadClusterResults(DATA, Clusters);
    setappdata(DATA.toplevel,'Clusters',Clusters);
    end
elseif isdir(name)
    d = dir(name);
    strings = {};
    xid = [];
    for j  = 1:length(d)
        if strfind(d(j).name,'ClusterTimes.mat') 
            strings = {strings{:} [name '/' d(j).name]};
            if strfind(d(j).name,'OldClusterTimes.mat')
                xid = [xid j];
            elseif strfind(d(j).name,'Copy of')
                xid = [xid j];
            end
        end
    end
    id = setdiff(1:length(strings),xid);
    strings = strings(id);
    Array = GetArrayConfig(name);
    if ~isempty(Array)
        args = {args{:} 'Array' Array};
    end
elseif exist(name,'file')
    d = mydir([fileparts(name) '/*ClusterTimes.mat']);
    strings = {d.name name};
%    strings{1} = name;
end

DATA = FixClusters(strings,args{:});

if ~isempty(DATA.errorlist)
    errorfile = [name '/ClusterErrors.mat'];
    if exist(errorfile)
        E = load(errorfile);
    end
    errorlist = DATA.errorlist;
    save(errorfile,'errorlist');
end


function [DATA, errs] = FixClusters(strings, varargin)
    useman = 1;
    loadexpts = 1;
    writemode = 0; %0 = no write, 1 = write with confrim, 2 = write all
    timecheck = now;
    verbose = 1;
    j = 1;
    cmpdrive = [];
    cmp = [];
    checks = [0 1 1];
    checknext = 0;
    checkmean = 0;
    clstcheck = 0;
    copyallnewer = 0;
    checkauto = 0;
    checkexpt = 0;
    checkdetails = 0; 
    checktrigger = 0;
    
    fixlatestbug = 0;
    filldummy = 0;
    checkmark = 0;
    showolderrs = 1;
    
    detailsbackup = [];
    findspace = [];
    sublist = [];
    Expts = [];
    Array = [];
    eids = [];
    DATA.saves = {};
    savedetails = 0;
    saveclusters = 0;
    removetrigset = 0;

while j <= length(varargin)
    if strncmpi(varargin{j},'Array',6)
        j = j+1;
        Array = varargin{j};
    elseif strncmpi(varargin{j},'templatesrc',11)
        j = j+1;
        DATA.templatesrc = varargin{j};
    elseif strncmpi(varargin{j},'dummy',5)
        j = j+1;
        filldummy = 1;
        writemode = 2;
        saveclusters = 1;
    elseif strncmpi(varargin{j},'cmpdrive',6)
        j = j+1;
        cmpdrive =  varargin{j};
    elseif strncmpi(varargin{j},'clstcheck',6)
        clstcheck =  1;
    elseif strncmpi(varargin{j},'checkauto',9)
        checkauto = 1;
    elseif strncmpi(varargin{j},'checkexpt',9)
        checkexpt = 1;
    elseif strncmpi(varargin{j},'checknext',9)
        checknext = 1;
    elseif strncmpi(varargin{j},'checkmean',9)
        checkmean = 1;
    elseif strncmpi(varargin{j},'checktrigger',9)
        checktrigger = 1;
        showolderrs = 0;
    elseif strncmpi(varargin{j},'checkmark',9)
        checkmark = 1;
        saveclusters = 1;
    elseif strncmpi(varargin{j},'copynew',7)
        copyallnewer = 1;
    elseif strncmpi(varargin{j},'details',7)
        checkdetails = 1;
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        if isnumeric(varargin{j})
            sublist = varargin{j};
        else
            Expts = varargin{j};
        end
    elseif strncmpi(varargin{j},'latest',6)
        fixlatestbug = 1;
    elseif strncmpi(varargin{j},'sublist',7) % just checkcertain files
        j = j+1;
        sublist = varargin{j};
    elseif strncmpi(varargin{j},'findspace',6)
        j =j+1;
        findspace = varargin{j};
    elseif strncmpi(varargin{j},'quiet',5)
        verbose = 0;
    elseif strncmpi(varargin{j},'timecheck',6)
        j = j+1;
        timecheck =  varargin{j};
    elseif strncmpi(varargin{j},'savedetails',9)
        savedetails = 1;
    elseif strncmpi(varargin{j},'switchtrigger',9)
        removetrigset = 2;  %special case for fixing pauls trigset{1} -> 0
        checktrigger = 1;
        saveclusters = 1; %will need saving if found;
    elseif strncmpi(varargin{j},'write',5)
        writemode = 2;
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        TOPTAG = varargin{j};
    end
    j = j+1;
end

    errs = {};

    AllClusters = {};
    AllCClusters = {};
    errorlist = [];
   
   msgmode = 0;
   newlog = [];
   if length(cmpdrive) && exist([cmpdrive '/ClusterLog.mat'],'file')
       newlog = load([cmpdrive '/ClusterLog.mat']);
       newlog.probes = CellToMat(newlog.ClusterLog,'probe');
       newlog.expts = CellToMat(newlog.ClusterLog,'exptno');
       newlog.savetime = CellToMat(newlog.ClusterLog,'savetime');
   end
   
   if checkexpt && isempty(Expts)
       [a, Expts] = APlaySpkFile(strings{end});
   end
   if checkexpt && isempty(strfind(strings{end},'ClusterTimes'))
       strings = strings(1:end-1);
   end
   exptlist = [];
   if ~isempty(sublist)
       exids = GetExptNumber(strings);
       id = find(ismember(exids,sublist));
       strings = strings(id);
   end

    for j = 1:length(strings)
        nerr = length(errs);
        eid = GetExptNumber(strings{j});
        DATA.exptno(j) = eid;
        if ischar(strings{j})
        dname = strrep(strings{j},'.mat','Details.mat');
        if ~exist(dname,'file')
            dname = strrep(dname,'ClusterTimes','ClusterTimesDetails.mat');
        end
        if checkdetails
            bname = strrep(dname,'/Expt','/backup/Expt');
            bname = strrep(bname,'.mat','*');
            db = dir(bname);
            if ~isempty(db)
                [a,b] = sort([db.datenum]);
                s = [fileparts(bname) '/' db(b(end)).name];
                fprintf('Loading Backup %s\n',s);
                detailsbackup = load(s);
            end
        end

        go = 1;
        if isempty(strfind(strings{j},'ClusterTimes'))
            go = 0;
        end
        if go
            if length(cmpdrive)
                a = fileparts(strings{j});
                afile = strrep(strings{j},a,cmpdrive);
                bfile = strrep(dname,a,cmpdrive);
                if exist(afile,'file')
                    if verbose
                        fprintf('Loading %s\n',afile);
                    end
                    load(afile);
                    cmp.Clusters = Clusters;
                    clear Clusters;
                    load(bfile);
                    cmp.ClusterDetails = ClusterDetails;
                    clear ClusterDetails;
                else
                    cmp.Clusters = [];
                    cmp.ClusterDetails = [];
                end
            else
                cmp.Clusters = [];
            end
            if verbose
                fprintf('Loading %s\n',strings{j});
            end
            load(strings{j});
            autofile = length(strfind(strings{j},'AutoCluster'));
            if exist(dname,'file')
                load(dname);
            else
                ClusterDetails = {};
            end
            AllClusters{j} = Clusters;
            if isempty(exptlist)
            exptlist = 1:length(Clusters);
            end
            AllCClusters{j} = cmp.Clusters;
        else
            Clusters = {};
        end
        else
            Clusters = strings{j};
        end
        
        modified = 0;
        if diff(size(exptlist)) <0
            exptlist = exptlist';
        end
        nerr = 0;
        if checkexpt && ~isempty(Expts)
            if isfield(Clusters{1},'usealltrials') && Clusters{1}.usealltrials
                fprintf('Calling APlaySpkFile with usealltrials');
                if ~isfield(Expts{1}.Header,'usebadtrials') || Expts{1}.Header.usebadtrials == 0
                    [a, Expts] = APlaySpkFile(Expts{1}.Header.loadname,'usealltrials');
                end
            end
            t = minmax(Clusters{1}.times);
            for e = 1:length(Expts)
                matched(e) = 0;
                if t(1) > Expts{e}.Header.trange(1)./10000 && t(2) < Expts{e}.Header.trange(2)./10000 
                    fprintf('%s Matches Expt %d\n',strings{j},e);
                    matched(e) = 1;
                elseif t(1) > Expts{e}.Header.trange(1)./10000 && t(1) < Expts{e+1}.Header.trange(1)./10000
                    cprintf('red','%s Bigger than  Expt %d\n',strings{j},e);
                    matched(e) = 2;
                end
            end
            if sum(matched) == 0
                    cprintf('red','%s No matching Expt\n',strings{j});                    
            end
        end
        foundCD = 0;
        for k = 1:length(Clusters)
            cmodified = 0;
            errtype = 0;
            C = Clusters{k};
            if length(findspace) && isfield(C,'space') && length(C.space) >= length(findspace)
                a = sum(findspace == C.space(1:length(findspace)));
                if a == length(findspace)
                    fprintf('E%.0fP%d space is %s\n',C.exptno,k,sprintf(' %d',C.space));
                end

            end
            if isfield(C,'marked') && C.marked == 3
                fprintf('E%dP%d marked bad',eid,k);
                if isfield(Array,'badprobes') && ~ismember(k,Array.badprobes);
                    fprintf(' But not in ArrayConfig.badprobes\n');
                    if checkmark
                        Clusters{k}.marked = 0;
                        modified = 1;
                        saveclusters = 1;
                    end
                end
            end
            if length(ClusterDetails) >= k && isfield(ClusterDetails{k},'clst')
                D = ClusterDetails{k};
                if clstcheck && size(D.clst,1) < size(D.clst,2)
                fprintf('P%d Clst is %dx%d %s\n',k,size(D.clst,1),size(D.clst,2),datestr(C.savetime(1)));
                ClusterDetails{k}.clst = ClusterDetails{k}.clst';
                cmodified = 1;
                modified = 1;
                end
            end
            if fixlatestbug && Clusters{k}.trigdt == 4
                Clusters{k}.trigdt = 0;
                Clusters{k}.quick = 1;
                modified = 1;
                saveclusters = 1;
            end
            if checktrigger && isfield(C,'trigset')
                fprintf('E%dP%d extra Trigset\n',C.exptno, C.probe);
                [a,b] = GetTriggerSet(C,'warn');
                if ~isempty(b.duplicates)
                    s = sprintf('E%dP%d Duplicate TriggerCluster',C.exptno,C.probe(1));
                    errs = {errs{:} s};
                end
                if removetrigset == 2
                    modified = modified+1; 
                   C = CopyFields(C,C.trigset{1});
                   C.triggerset = 0;
                   for c = 1:length(C.next)
                       C.next{c}.triggerset = 0;
                   end
                   C = rmfield(C,'trigset');
                   Clusters{k} = C;
                end
            end
            if length(cmpdrive) && length(cmp.Clusters) >= k
                CC = cmp.Clusters{k};
                copynew = 0;
                if (CC.savetime(1) > C.savetime(1) || CC.savetime(1) > timecheck)  && CC.auto == 0
                    s =sprintf('Newer: %s P%d (%s) (%s)',afile,k,datestr(CC.savetime(1)),datestr(C.savetime(1)));
                    id = find(newlog.probes == C.probe & newlog.expts == C.exptno & newlog.savetime(:,1) > C.savetime(1)+0.1);
                    if length(id)
                        rectype = CellToMat(newlog.ClusterLog(id),'recluster');
                        L = newlog.ClusterLog{id(end)};
                        s = sprintf('%s re%d on %s',s,min(rectype(:,1)),L.hostname);
                        if sum(rectype(:,1) ==0) > 0  % a new manual cut was made
                            copynew = 1;
                        elseif copyallnewer
                            copynew = 1;
                        end
                    end
                    errtype = errtype+1;
                    
                    
                    if copynew
                        s = strrep(s,'Newer','Newer**');
                        Clusters{k} = cmp.Clusters{k};
                        ClusterDetails{k} = cmp.ClusterDetails{k};
                        C = Clusters{k};
                        cmodified = 1;
                        modified = 1;
                        savedetails = 1;
                        saveclusters = 1;
                        errtype = errtype+2;
                    end
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                end
                if length(ClusterDetails) >= k && isfield(ClusterDetails{k},'clst')
                    D = ClusterDetails{k};
                    if k > length(cmp.ClusterDetails) || isempty(cmp.ClusterDetails{k})
                        s=sprintf('%s:%d Missing ClusterDetails',afile,k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                        errtype = errtype+4;
                    else
                    DC = cmp.ClusterDetails{k};
                    a = length(D.clst);
                    b = length(DC.clst);
                    if a ~= b
                        s=sprintf('%s:%d  Clst length mistmatch %d vs %d',afile,k,a,b);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    end
                    [a,b] = Counts(D.clst);
                    [c,d] = Counts(DC.clst);
                    if length(c) ~= length(a)
                        s=sprintf('%s:%d  Nclusters mismatch',afile,k,a,b);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    end
                    end
                    a = length(D.xy);
                    b = length(DC.xy);
                    if a ~= b
                        s=sprintf('%s:%d  xy length mistmatch %d vs %d\n',afile,k,a,b);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};

                    else
                        xc = corrcoef(D.xy(:,1),DC.xy(:,1));
                        yc = corrcoef(D.xy(:,2),DC.xy(:,2));
                        if yc(1,2) < 0.99 || xc(1,2) < 0.99
                            s=sprintf('%s:%d  xy corr bad %.2f %.2f vs %d\n',afile,k,xc(1,2),yc(1,2));
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                        end
                    end

                    
                end
                if isfield(C,'next') && isfield(CC,'next')
                    if length(C.next) ~= length(CC.next)
                        s=sprintf('%s:%d  next length mistmacth)\n',afile,k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    end
                elseif isfield(CC,'next')
                        s=sprintf('Next Missing %s:%d',strings{j},k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                elseif isfield(C,'next')
                        s=sprintf('Next Missing %s:%d',afile,k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                end
                    
            end
            
            if isfield(C,'errs') && ~isempty(C.errs)              
                nerr = nerr+1;
                errorlist(nerr).errs = C.errs;
                errorlist(nerr).p = C.probe;
                errorlist(nerr).ex = C.exptno;     
                errorlist(nerr).auto = C.auto;     
                if showolderrs
                    for e = 1:length(C.errs)
                        cprintf('blue','E%dP%d %s',C.exptno,C.probe(1),deblank(C.errs{e}));
                    end
                    fprintf('\n');
                end
            end
            if isfield(C,'next') && ~iscell(C.next)
                fprintf('Converting C2 next to cell %s,%d\n',strings{j},k);
                next = Clusters{k}.next;
                Clusters{k} = rmfield(Clusters{k},'next');
                Clusters{k}.next{1} = next;
                modified = 2;
                cmodified = 1;
            end
            if checknext && isfield(C,'next') && iscell(C.next) && length(C.next) > 0
                if isfield(C.next{1},'next') 
                    if ~isfield(C.next{1},'MeanSpike') && isfield(C.next{1}.next{1},'MeanSpike')
                        C.next{1} = C.next{1}.next{1};
                        s= sprintf('Suspicious Nested Next Cells %s,%d',strings{j},k);
                    else
                        s = sprintf('Nested Next Cells %s,%d',strings{j},k);
                    end
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                    Clusters{k}.next{1} = rmfield(Clusters{k}.next{1},'next'); 
                    modified = 1;
                    saveclusters = 1;
                end
                isset = 0;
                missing = [];
                for c = length(C.next):-1:1
                    if isfield(C.next{c},'mahal')
                        isset = 1;
                    elseif isset
                        missing = [missing c];
                    end
                    if checkmean && isset && ~isfield(C.next{c},'MeanSpike')
                       s= sprintf('NextMean %d) missing Cells %s,%d',c,strings{j},k);
                       fprintf('%s\n',s);
                       errs = {errs{:} s};
                       if isfield(CC.next{c},'MeanSpike')
                           saveclusters = 1;
                           modified = 1;
                           Clusters{k}.next{c} = rmfields(CC.next{c},'next');
                       end
                          
                    end
                end
                if length(missing)
                    s = sprintf('Missing nexts (%s) %s,%d',sprintf('%d,',missing),strings{j},k);
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                end
            end

            if checknext && isfield(C,'xtimes')  %check if left over from old way
                if isfield(C,'next') && iscell(C.next) && ~isempty(C.next)
                % left over from old way
                    if isfield(C.next{1},'times')
                    s = sprintf('Redundant xtimes  %s,%d',strings{j},k);
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                    Clusters{k} = rmfields(Clusters{k},'xtimes');
                    modified = 1;
                    saveclusters = 1;
                    else
                        s = sprintf('Unmatched xtimes %s,%d',strings{j},k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    end
                else
                    s = sprintf('Need xtimes %s,%d',strings{j},k);
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                end
            end


            if checkauto
                if ~isfield(C,'auto')
                    Clusters{k}.auto = 1;
                elseif C.auto == 2
                    if autofile
                        s = sprintf('PlotClusters in Auto  %s,%d',strings{j},k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    end
                    if C.shape == 0
                        ctype = 'E';
                    else
                        ctype = 'L';
                    end
                    s = sprintf('%s P%d made by plotclusters %s',strings{j},k,ctype);
                    fprintf('%s\n',s);
                    errs = {errs{:} s};
                end
                if ~isfield(C,'manual')
                    modified = modified+1;
                    saveclusters = 1;
                    s = sprintf('%s P%d Missing field manual',strings{j},k);
                    fprintf('%s\n',s);
                    errs = {errs{:} s};                    
                    if C.auto
                        Clusters{k}.manual = 0;
                    else
                        Clusters{k}.manual = 1;
                    end
                end

            end


            if length(ClusterDetails) >= k && ~isempty(ClusterDetails{k})
                foundCD = foundCD+1;
            if ~isfield(Clusters{k},'clst') && ~isfield(ClusterDetails{k},'clst')
                if Clusters{k}.nspks == length(ClusterDetails{k}.t) && length(Clusters{k}.times) == Clusters{k}.ncut
                    s=sprintf('**Need to reconstruct clst for %s,%d\n',strings{j},k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                    ClusterDetails{k}.clst = ones(size(ClusterDetails{k}.t));
                    ClusterDetails{k}.buildby = 'fix';
                    modified = modified+1;
                    cmodified = 1;
                    id = find(ismember(ClusterDetails{k}.t,Clusters{k}.times));
                    if length(id) == Clusters{k}.ncut
                        ClusterDetails{k}.clst(id) = 2;
                    else
                        fprintf('Length Mismatch\n');
                    end
                else
                s=sprintf('Missing clst for %s,%d\n',strings{j},k);
                        fprintf('%s\n',s);
                        errs = {errs{:} s};
                end
            elseif ~isfield(ClusterDetails{k},'clst')
                fprintf('**clst in Clusters for %s,%d\n',strings{j},k);
                if length(Clusters{k}.clst) == length(ClusterDetails{k}.t)
                    ClusterDetails{k}.clst = Clusters{k}.clst;
                    ClusterDetails{k}.buildby = 'fix';
                    modified = modified+1;
                    cmodified = 1;
                end
            end

            if ~isfield(Clusters{k},'sign')
                Clusters{k}.sign= 0;
            end
            if checks(1)
            if isfield(C,'excludetrialids')
                fprintf('%s P%d %d excluded trials\n',strings{j},k,length(C.excludetrialids));
            end
            if isfield(C,'restricttimerange')
                if min(C.times) > C.restricttimerange(2) || max(C.times) < C.restricttimerange(1)
                fprintf('**%s P%d BAD Restrited Time %.2f - %.2f (%.2f - %.2f)\n',strings{j},k,C.restricttimerange(1),C.restricttimerange(2),min(C.times),max(C.times))
                Clusters{k} = rmfield(Clusters{k},'restricttimerange');
                modified = 2;
                cmodified = 1;
                else
                fprintf('%s P%d OK Restrited Time %.2f - %.2f (%.2f - %.2f)\n',strings{j},k,C.restricttimerange(1),C.restricttimerange(2),min(C.times),max(C.times))
                end
            end
            end
            elseif checkdetails
                fprintf('ClusterDetails %d is empty/missing\n',k);
                 if ~isempty(detailsbackup) && length(detailsbackup.ClusterDetails) >= k && ~isempty(detailsbackup.ClusterDetails{k})
                    ClusterDetails{k} = detailsbackup.ClusterDetails{k};
                    modified = modified+1;
                end
            end
            if filldummy
                force = 1; %temp
                if ~isfield(C,'shape')|| force
                    Clusters{k}.space = [1 3];
                    Clusters{k}.shape = 0;
                    Clusters{k}.dvdt = 0;
                    Clusters{k}.csd = 0;
                    Clusters{k}.auto = 0;
                    Clusters{k}.dropi = [0 0 0];
                    Clusters{k}.time = now;
                    Clusters{k}.ctime = now;
                    if ~ isfield(Clusters{k},'MeanSpike') || isempty(Clusters{k}.MeanSpike)
                        Clusters{k}.MeanSpike.ms = ones(1,32);
                    end
                    if ~isfield(Clusters{k}.MeanSpike,'mu')
                        Clusters{k}.MeanSpike.mu = ones(1,32);
                    end
                    DATA.modified(j,k) = 1;
                    modified = modified+1;
                end
                saveclusters = 1;
            end
            DATA.modified(j,k) = cmodified;
            DATA.nerr(j,k) = errtype;
        end
        if checkdetails 
            fprintf('E%d Found Details for %d probes\n',GetExptNumber(strings{j}),foundCD);
            DATA.ndetails(j) = foundCD;
        end
        if modified &&  writemode > 0
            if writemode > 1
                yn = 'y';
            else
                yn = input(sprintf('%s was modified - Save to Disk?\n',strings{j}),'s');
            end
            if length(yn) && yn(1) == 'y'
                if savedetails
                    s = sprintf('Saving %s',dname);
                    DATA.saves = {DATA.saves{:} s};
                    fprintf('%s\n',s);
                    BackupFile(dname);
                    if exist('FullVData')
                        save(dname,'ClusterDetails','FullVData');
                        clear FullVData;
                    else
                        save(dname,'ClusterDetails');
                    end
                end
                if saveclusters
                    s=sprintf('Saving %s',strings{j});
                    DATA.saves = {DATA.saves{:} s};
                    fprintf('%s\n',s);
                    BackupFile(strings{j});
                    if exist('FullVData','var')
                        save(strings{j},'Clusters','FullVData');
                    else
                        save(strings{j},'Clusters');
                    end
                end
            end
        end
    end
    fprintf('Done\n');

    DATA.errorlist = errorlist;
    DATA.errs = errs;
    DATA.Clusters = AllClusters;
    DATA.CC = AllCClusters;