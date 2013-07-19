function res = RunAllGridFiles(name, varargin)
%RunAllGridFiles(file, .....
%Reads all FullV Matlab files associated with file 'name' and calls AllVPcs, for Utah Array Expts
%N.B. This ususally uses the .mat filename, not the directory, as all expts are included in one .smr file
%
%To apply quantification of quick saves, use
%RunAllGridFiles(file,'quantify','savespikes')  
%


expts = [];
bysuffix = 0;
scanfiles = 0;
recalc = 0;
sizecheck = 0;  %set to max safe bytes for FullV to control if FullV is kept in memory
runautocut = 1;
%memsize = CheckPhysicalMemory;
%sizecheck = (memsize * 1e6) .* 0.8;
X.spkrate = 50;
template = [];
probes = [];
forcerun = 0;
args = {};
addargs = {};
checktype = 'reclassify';
classifyfromcluster = 1;
parallel = 0;
listtasks = 0;



parid = 0;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
    elseif strncmpi(varargin{j},'args',4)
        j = j+1;
        addargs = varargin{j};
    elseif strncmpi(varargin{j},'checktype',9)
        j = j+1;
        checktype = varargin{j};
    elseif strncmpi(varargin{j},'classify',6)
        j = j+1;
        classifyfromcluster = 1;
        Clusters = varargin{j};
    elseif strncmpi(varargin{j},'quantify',6)
        checktype = 'quantify';
    elseif strncmpi(varargin{j},'clusters',6)
        j = j+1;
        Clusters = varargin{j};
    elseif strncmpi(varargin{j},'expts',3)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'listtasks',8)
        listtasks = 1;
    elseif strncmpi(varargin{j},'nocut',4)
        runautocut = 0;
        args = {args{:} 'nocut'};
    elseif strncmpi(varargin{j},'parallel',4)
        parallel = 1;
        parid = j;
    elseif strncmpi(varargin{j},'probes',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'refcut',6)
        classifyfromcluster = 3;
    elseif strncmpi(varargin{j},'reclassify',6)
        classifyfromcluster = 0;
        args = {args{:} 'reclassify'};
    elseif strncmpi(varargin{j},'run',3)
        forcerun = 1;
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end
res.starts(1) = now;
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:MaxIterations');

if iscellstr(name) && isdir(name{1}) %given a list of directories
    if parallel
        varargin = varargin(setdiff(1:length(varargin),parid));
        jobs = {};
        groups = [];
        parfor (j = 1:length(name))
            t = mygetCurrentTask();
            fprintf('Worker %d scanning %s\n',t.ID,name{j})
            dres{j} = RunAllGridFiles(name{j},varargin{:},'listtasks');
            if ~isfield(dres{j},'names')
                dres{j}.names = {};
            end
        end
        tic;
        nf = 0;
        ng = 0;
        for j = 1:length(dres)
            lastex = 0;
            for k = 1:length(dres{j}.names)
                nf = nf+1;
                if dres{j}.exids(k) ~= lastex
                    lastex = dres{j}.exids(k);
                    ng = ng+1;
                end
                jobs{nf}{1} = dres{j}.names{k};
                jobs{nf}{2} = dres{j}.cfiles{k};
                cfromc(nf) = dres{j}.classifyfromcluster;
                probes(nf) = dres{j}.probes(k);
%                jobs = {jobs{:} {dres{j}.names{k} dres{j}.cfiles{k} dres{j}.classifyfromcluster}};
                groups(nf) = ng;
            end
        end
        toc
        clear dres;
        classifyfromcluster = 2;
        exgroups = unique(groups);
        probe = unique(probes);
        testing = 1;
        if testing
            for k = 1:length(probe)
                id = find(probes == probe(k));
                parfor (j = 1:length(id))
                    t = mygetCurrentTask();
                    job = jobs{id(j)};
                    fprintf('Worker %d working on %s\n',t.ID,job{1})
                    dres{j}{k} = DoFile(job{1},job{2},cfromc(id(j)),args{:});
                end
            end
        else
                    
        parfor (j = 1:length(exgroups))
            id = find(groups == exgroups(j));
            for k = 1:length(id)
                t = mygetCurrentTask();
                job = jobs{id(k)};
                fprintf('Worker %d working on %s\n',t.ID,job{1})
                dres{j}{k} = DoFile(job{1},job{2},cfromc(id(k)),args{:});
            end
%            dres{j} = RunAllGridFiles(name{j},varargin{:});
        end
        end
        res.dres = dres;
    else
        for j = 1:length(name)
            res.dres{j} = RunAllGridFiles(name{j},varargin{:});
        end
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
    elseif isfield(name,'exptlist')
        if isfield(name,'args') &&length(addargs)
            name.args = {name.args{:} addargs{:}};
        end
        D = CheckExptList(name);
    end
    return;
 elseif isdir(name)
    path = name;
else
[path, fname] = fileparts(name);
end

if listtasks
    logfid = -1;
elseif ischar(name) 
    logfile = [name '/' 'RunAll.log'];
    logfid = fopen(logfile,'a');
elseif isfield(name,'prefix')
    logfile = [name.prefix '/' 'RunAll.log'];
    logfid = -1;
else
    logfid = -1;
end
X.logfid = logfid;


if isdir(name)
res.prefix = name;
else
    res.prefix = fileparts(name);
end
if strcmp(checktype,'quantify')
    d = dir([res.prefix '/*ClusterTimes.mat']);
    ns = 1;
    needed = 0;
    vfiles = {};
    for j = 1:length(d)
        name = [path '/' d(j).name];
        load(name);
        ex = sscanf(d(j).name(5:end),'%d');
        for c = 1:length(Clusters)
            need = 0;
            if (isfield(Clusters{c},'quick') && Clusters{c}.quick > 0) || Clusters{c}.manual == 2
                need = 1;
                fprintf('%s need Probe %d C1\n',d(j).name,c);
            elseif length(Clusters{c}.next) > 0
                for k = 1:length(Clusters{c}.next)
                if isfield(Clusters{c}.next{k},'manual') && Clusters{c}.next{k}.manual == 2
                    need = 1;
                    fprintf('%s need Probe %d C%d\n',d(j).name,c,k+1);
                end
                end
            end
            if need
                needed = needed+1;
                vfile = [res.prefix '/Expt' num2str(ex) '.p' num2str(c) 'FullV.mat'];
                vfiles{needed} = vfile;
                cids(needed) = c;
                exids(needed) = ex;
            end
        end
        if needed == 0
            fprintf('%s is up to date\n',d(j).name);
        end
    end
    if parallel
        parfor (j = 1:needed)
            try
              cls{j} = AllVPcs(vfiles{j},'tchan',cids(j),'reclassify',args{:},'noninteractive');
            catch ME
                fprintf('Error!!!! %s:%d %s\n',vfiles{j},cids(j),ME.message)
            end
        end
    elseif listtasks
        for j = 1:needed
            cls{j} = { exids(j) vfiles{j} cids(j) 'reclassify' args{:} };
        end
    else
        for j = 1:needed
            cls{j} = AllVPcs(vfiles{j},'tchan',cids(j),'reclassify',args{:},'noninteractive');
        end
    end
    res = cls;
    return;
end

d = dir([name '/*.p*FullV.mat']);
ns = 1;
DistanceMatrices = {};
needfile = [];
for j = 1:length(d)
    sid = regexp(d(j).name,'Expt[0-9,a]*.p[0-9]*FullV');
    X.logfid = logfid;
    if length(sid)
        ex = sscanf(d(j).name(sid(1)+4:end),'%d');
        pid = regexp(d(j).name,'.p[0-9]*FullV');
        p = sscanf(d(j).name(pid(1)+2:end),'%d');
        if classifyfromcluster == 3 %ref cut
        cfile = [res.prefix '/RefClusters.mat'];
        else
        cfile = [res.prefix '/Expt' num2str(ex) 'ClusterTimes.mat'];
        end
        cfiles{j} = cfile;
        names{j} = [path '/' d(j).name];
        if (isempty(expts) || ismember(ex,expts)) && (isempty(probes) || ismember(p,probes))
            needfile(j) = 1;
        else
            needfile(j) = 0;
        end
        allexpts(j) = ex;        
        probes(j) = p;
    end
end

id = find(needfile);
if isempty(id)
     if logfid > 0
     fclose(logfid);
     end
    return;
end

if isempty(expts)
    expts = allexpts;
end
    x = zeros(size(needfile));
    for j = 1:length(id)
        x(id(j)) = j;
    end
    ex = unique(expts);
    exids = {};
    for j = 1:length(ex)
        exids{j} = id(allexpts(id) == ex(j));
    end
    cls = {};
if parallel  %have to do this by expt. If two probes in one expt done at the same
    %time, can both try to write file and corrupt it.
    fprintf('Calling parfor with %d expts\n',length(exids));
    parfor (e = 1:length(exids))
        efile = strrep(logfile,'RunAll','ExptRun');
        fid = fopen(efile,'a');
        t = getCurrentTask();
        for k = 1:length(exids{e})
            j = exids{e}(k);
            fprintf('Worker %d Processing %s %s\n',t.ID,d(j).name,datestr(now));
            jj = x(j);
            fprintf(fid,'Worker %d Processing %s %s\n',t.ID,d(j).name,datestr(now));
            cls{e}{k} = DoFile(names{j},cfiles{j},classifyfromcluster,args{:});
        end
        fclose(fid);
    end
    res.cls = cls;
elseif listtasks
    res.exids = allexpts;
    res.names = names;
    res.cfiles = cfiles;
    res.args = args;
    res.probes = probes;
    res.classifyfromcluster = classifyfromcluster;
    return;
else
    res.cls =  {};
    for (e = 1:length(exids))
        for k = 1:length(exids{e})
            j = exids{e}(k);
            fprintf('Processing %s %s\n',d(j).name,datestr(now));
            fprintf(logfid,'Processing %s %s\n',d(j).name,datestr(now));
            res.cls{x(j)} = DoFile(names{j},cfiles{j},classifyfromcluster,args{:});
        end
    end
end


for j = 1:length(res.cls)
if isfield(res.cls{j},'DistanceMatrix')
         DistanceMatrices{j} = rmfields(res.cls{j}.DistanceMatrix,'clid');
end
res.names{j} = names{id(j)};
res.exptlist(j) = GetExptNumber(d(id(j)).name);
end  

if ~isempty(DistanceMatrices)
    dname = [path '/DistanceMatrix.mat']; 
    fprintf('Saving %s\n',dname);
    save(dname, 'DistanceMatrices');
end
 res.end = now;
 res.args = args;
 if logfid > 0
     fclose(logfid);
 end
 
 
 function res = DoFile(name, cfile, ctype, varargin)
     res = [];
     if ctype ==2
         if exist(cfile)
             load(cfile);
             if (isfield(Clusters{p},'quick') && Clusters{p}.quick) || ...
                     (isfield(Clusters{p},'manual') && Clusters{p}.manual == 2)
                 res = ClassifyFromCluster(name,Clusters,'quantify',varargin{:});
             end
         end
     elseif ctype == 3
         try
             load(cfile);
             res = ClassifyFromCluster(name,Clusters,'reapply',varargin{:});
         catch ME
         end
     elseif ctype == 1 || ctype == 3
         load(cfile);
         res = ClassifyFromCluster(name,Clusters,varargin{:});
     elseif ~exist(cfile)
         fprintf('No Cluster File %s\n',cfile);
     else
         try
             load(cfile);
             try
                 res = AllVPcs(name, varargin{:},'noninteractive');
             catch ME
                 fprintf('!!!AllVPcs Error %s for file %s (line %d, m-file %s)\n',ME.message,name,ME.stack(1).line,ME.stack(1).name);
                 res = ME;
             end
         catch ME
             fprintf('!!!Error Reading %s\n',cfile);
             res = ME;
         end
     end
     if isfield(res,'logfid') && res.logfid > 0
         try
             fclose(res.logfid);
             fprintf('Successfully closed %d\n',res.logfid);
         catch
             fprintf('FID %d could not be closed\n',res.logfid);
         end
     else
         fprintf('No Log File in result for %s\n',name);
     end

 
 function res = ClassifyFromCluster(name, Clusters, varargin)
 
    requantify = 0;
    argon = {};
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'requantify',6)
            requantify = 1;
        elseif strncmpi(varargin{j},'reclassify',6)
            reclassify = 1;
        else
            argon = [argon varargin{j}];
        end
        j = j+1;
    end
    pid = regexp(name,'.p[0-9]*FullV');
    p = sscanf(name(pid(1)+2:end),'%d');
    args = {};
    if isfield(Clusters{p},'gmdprime') && Clusters{p}.gmdprime > 3
        args = {args{:} 'refine'};
    end
    res = AllVPcs(name,'tchan',p,'reapply',Clusters{p},args{:},argon{:},'noninteractive');

 function Summary = CheckExptList(X)
             %if args starts by listing all probes, strip this - want to call
        %one at at time next

        
        for j = 1:length(X.cls)
            for k = 1:length(X.cls{j})
                if isfield(X.cls{j}{k},'cluster')
                    Summary.gmdprime(j) = X.cls{j}{k}.cluster.gmdprime;
                    Summary.mahal(j,:) = X.cls{j}{k}.cluster.mahal;
                else
                    fprintf('%d.%d is empty\n',j,k);
                end
            end
        end
            
function result = CheckReclassify(X, checktype, varargin)
    
    exptno = 0;
    j = 1;
    while j <= length(varargin)
        if strncmp(varargin{j},'exptno',6)
            j = j+1;
            exptno = varargin{j};
        end
        j = j+1;
    end
    if isfield(X,'cls')
    [exps, exi] = sort(X.exptlist);
        acts = [];
        for j = 1:length(X.cls)
            e = X.exptlist(j);
            e = find(exps == X.exptlist(j));
            for k = 1:length(X.cls{j})
                if ~isfield(X.cls{j}{k},'probes')
                    X.cls{j}{k}.probes = k;
                end
                if ~isfield(X.cls{j}{k},'exptno')
                    X.cls{j}{k}.exptno = j;
                end
                X.cls{j}{k}.exptid = e;
                if ~isfield(X.cls{j}{k},'args')
                    X.cls{j}{k}.args = X.args;
                end
            end
            res = CheckReclassify(X.cls{j},checktype);
            result.actimage(e,:) = res.actimage(end,:);
            result.dateimage(e,:) = res.ctime;
            result.auto(e,:) = res.auto;
            result.xclim(e,:) = res.xcl;
            result.name{e} = X.name{j};
            if isfield(res,'readmethod')
                result.readmethod(e,:) = res.readmethod;
            end
            acts = [acts res.acts];
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
    for j = 1:length(X) 
        if iscell(X)
            res(j).result = CheckResult(X{j}, checktype);
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
            if isfield(X{j}.cluster,'exptreadmethod')
                res(j).readmethod = X{j}.cluster.exptreadmethod;
            end
            res(j).auto = X{j}.cluster.auto;
        else
            res(j).result = CheckResult(X(j),checktype);
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
        actimage(res(j).exptno,res(j).probes) = res(j).result+1;
    end
    result.acts = res;
    result.actimage = actimage;
    result.ctime = [res.ctime];
    result.auto = [res.auto];
    result.xcl = [res.xcl];
    if isfield(res,'readmethod')
    result.readmethod = [res.readmethod];
    end
    
function res = CheckResult(a, type)
    if a.auto == 2 %old method for recording plotcluster cuts.
        a.auto = 2;
    end
      if strcmp(type,'reclassify')
          if a.err >0 && a.auto == 1
              res = 1;
          elseif a.err > 0
              fprintf('E%dP%d  err, not auto\n',a.cluster.exptno,a.cluster.probe(1)) 
              res = 0;
          elseif abs(a.xcorr(1)) > 0.95 && abs(a.xcorr(2)) > 0.95 && a.matchcounts(1)./a.matchcounts(2) < 0.15 && a.overlapn./a.matchcounts(2) > 0.5
              res = 4-a.auto; %4 or 3 or 2
          elseif a.auto ==1
              res = 1;
          elseif a.matchcounts(1) ==0 && a.matchcounts(3) == 0
              res = 5;
          else
              fprintf('E%dP%d xc %.2f %.2f overlap %d(%s)\n',...
                  a.cluster.exptno,a.cluster.probe(1),a.xcorr(1),a.xcorr(2),...
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
        