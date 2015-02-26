function CD = CheckExptClusters(dname, varargin)
% CheckExptClusters(D, varargin) checks for misisng trials in Clusters
% if D is a folder, checks all Clusters and Expts.
%
% CheckExptClusters({D1...Dn}, varargin) checks for misisng trials in Clusters
% Does all folders named in a cell string array 
%CheckExptClusters(Clusters, Expts) Checks all elements of Clusters
%for which there is a matching Expt in Expts
%
%CheckExptClusters(..., 'fix')
%Adds any missing trials to the list in CellDetails 
%
%CellDetails = CheckExptClusters(...  returns a modified celldetails struct
%with the new list of excluded trials, regardless of whether its saved to
%disk
%

eid = [];
Expts = {};
CD  = [];
fixlist = 0;
check = '';
checkfullv = 0;
verbose = 0;
clearoldexclusions = 0;

j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        X = varargin{j};
        [a,b] = CellToMat(X,'Header');
        if sum(b.found) > 0
            Expts = X;
        end
    elseif isstruct(varargin{j})
        if isfield(varargin{j},'CellDetails')
            cells = varargin{j};
            CD = cells.CellDetails;
        elseif isfield(varargin{j},'trialids')
            CD = varargin{j};
        end
    elseif strncmpi(varargin{j},'checkfullv',10)
        checkfullv = 1;
    elseif strncmpi(varargin{j},'checkdetails',10)
        check = 'details';
    elseif strncmpi(varargin{j},'check',5)
        check = 'done';
    elseif strncmpi(varargin{j},'clearoldexclusions',7)
        clearoldexclusions = 1;
    elseif strncmpi(varargin{j},'trials',5)
        check = 'trials';
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        eid = varargin{j};
    elseif strncmpi(varargin{j},'fix',3)
        [j, fixlist] = GetOptionalArg(varargin,j,1);
    end
    j = j+1;
end

if ~isempty(Expts)
    exptid = GetExptNumber(Expts);
end
cdetails = [];

checktime = now;

if iscellstr(dname)
    for j = 1:length(dname)
        if strcmp(check,'done')
            CD{j} = CheckExptClusterDone(dname{j},varargin{:});
        elseif sum(strcmp(check,{'details' 'trials'}))
            CD{j} = CheckExptClusterDone(dname{j},varargin{:}, check);
        else
            CD{j} = CheckExptClusters(dname{j},varargin{:});
        end
    end
    return;
elseif iscell(dname) && iscell(dname{1}) %cell array of Cluster Cell arrays
    if isempty(Expts)
    end
    Clusters = dname;
    exptid = GetExptNumber(Expts);
    cexpts = prctile(CellToMat(Clusters,'exptno')',50);
    if ~isempty(eid)
        [id, cid] = ismember(exptid,cexpts);
        id = find(id);
        Clusters = Clusters(cid(id));
        cexpts = cexpts(cid(id));
        id = find(ismember(exptid,eid));
        Expts = Expts(exptid(id));
    else
        eid = exptid;
        id = find(ismember(cexpts,exptid));
        Clusters = Clusters(id);
        cexpts = cexpts(id);
    end
    id = find(ismember(exptid,eid));
    if fixlist
        cellname = [fileparts(Expts{1}.Header.loadname) '/CellList.mat']
        if exist(cellname)
            cells = load(cellname);
            CD = cells.CellDetails;
        else
            fixlist = 0;
        end
    end
elseif iscell(dname) %single set of clusters, for 1 expt
    while isempty(dname{1})
        dname = dname(2:end);
    end
    if isfield(dname{1},'')
        return;
    end
    if iscell(dname{1}) || isfield(dname{1},'minenergy')
        Clusters{1} = dname;
        cexpts = GetExptNumber(Clusters{1});
        for j = 1:length(Clusters{1})
            if ~isfield(Clusters{1}{j},'t')
                Clusters{1}{j}.t = Clusters{1}{j}.times;
            end
        end
    elseif isstruct(dname{1}) && isfield(dname{1},'trialids') %previos resuls
        for j = 1:length(dname)
            CD = ShowResult(dname{j},fixlist);
        end
        return;
    elseif isstruct(dname{1}) %combine AllClsuters struct
        cexpts = GetExptNumber(dname);        
        for j = 1:length(dname)
            for k = 1:length(dname{j})
            Clusters{j}{k} = dname{j}(k);
            if isfield(dname{j},'t')
                Clusters{j}{k}.t = dname{j}(k).t;
            elseif isfield(dname{j},'times')
                Clusters{j}{k}.t = dname{j}(k).times./10000;
            end
            Clusters{j}{k}.probe = k;
            end
        end
    end
    id = find(ismember(exptid,cexpts));
    Expts = Expts(id);
    if fixlist
        cellname = [fileparts(Expts{1}.Header.loadname) '/CellList.mat'];
        if exist(cellname)
            cells = load(cellname);
            CD = cells.CellDetails;
        else
            fixlist = 0;
        end
    end    
    dname = 'Passed Struct';
elseif isdir(dname)
    if isempty(Expts)
        Expts = ReadExptDir(dname);
    end
    exptid = GetExptNumber(Expts);
    if ~isempty(eid)
        id = find(ismember(exptid,eid));
        
        Expts = Expts(id);
    else
        id= [];
    end
    [Clusters, FullVData, cdetails] = LoadCluster(dname, GetExptNumber(Expts),'array');
    if length(id) == 1 && 0 %should be done by 'array' arg
        AllClusters{1} = Clusters;
        Clusters = AllClusters;
        a{1} = cdetails;
        cdetails = a;
    end
    cexpts = prctile(CellToMat(Clusters,'exptno')',50);
    cellname = [dname '/CellList.mat'];
    if exist(cellname)
        cells = load(cellname);
        CD = cells.CellDetails;
    end
    for j = 1:min([length(cdetails) length(Clusters)])
        needed = sum(CellToMat(Clusters{j},'needed'));
        if needed && isfield(cdetails{j},'errs')
            CD = AddError(CD, cdetails{j}.errs);
        end
    end
end

exptid = GetExptNumber(Expts); %selected and available
modified = 0;
excludetrials = {};
nprobes = 1;
if ~exist('cells')
    if ischar(dname)
        fprintf('No CellList file in %s\n',dname);
    end
    return;
end

CD.checktime = checktime;
CD.program = 'CheckExptClusters';
if isfield(CD,'excludetrials')
    oldxcl = CD.excludetrials;
else
    oldxcl = {};
end
for j = 1:length(Expts);
    cid = find(cexpts == exptid(j));
   sumerrs(j,1:4) = 0;
   if ~isempty(Expts{j}) && ~isempty(cid)
       C = Clusters{cid};
       E = Expts{j};
       clusterfile = '';
       if length(cdetails) >= cid && isfield(cdetails{cid},'loadname')
           clusterfile = cdetails{cid}.loadname;
       end
       clustermodified = 0;
       if isfield(CD,'exptids')
           cxid = find(CD.exptids == exptid(j));
       else
           cxid = [];
       end

       if clearoldexclusions
           oldx = CD.excludetrials{cxid,1,1};
           for p = 1:size(CD.excludetrials,2)
               for c = 1:size(CD.excludetrials,3)
                   oldx = intersect(oldx,CD.excludetrials{cxid,p,c});
               end
           end
           if ~isempty(oldx)
               modified = modified+1;
               for p = 1:size(CD.excludetrials,2)
                   for c = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,c} = setdiff(CD.excludetrials{cxid,p,c},oldx);
                   end
               end
           end
       end
       clear n;
       nt = length(E.Trials);
       clear starts;
       clear ends;
       for t = 1:length(E.Trials)
           starts(t) = E.Trials(t).Start(1);
           ends(t) = E.Trials(t).End(end);
           ids(t) = E.Trials(t).id;
       end
       missed = {};
       for c = 1:length(C)
           if isfield(C{c},'missingtrials')               
               missed{c} = C{c}.missingtrials(:);
           end
       end
       missingids = unique(cat(1,missed{:}));
       if ~isempty(missingids)
               fprintf('Expt%d Known Missing Trials %s\n',exptid(j),sprintf('%d ',missingids));
       end
       misserr = 0;
       for c = 1:length(missed)
           if sum(~ismember(missed{c},missingids))
               misserr(c) = 1;
           end
       end
       if sum(misserr)
           bid = find(misserr);
           fprintf('Probes %s incomplete missing Trials (%s)\n',sprintf(' %d',bid),sprintf(' %d',missingids))
       end
       expname = Expt2Name(E);
       CD.exptnames{exptid(j)} = expname;
       probes = CellToMat(C,'probe');
       nprobes = max(probes);
       for c = 1:length(C)
           if isfield(C{c},'errs')
               for k = length(C{c}.errs):-1:1
                   id = strfind(C{c}.errs{k},'Missing Trials');
                   if ~isempty(id)
                       xids = sscanf(C{c}.errs{k}(id(1)+14:end),'%d');
                       if sum(~ismember(xids,ids))
                           fprintf('Id %d not in Expt (%d), but has error message %s\n',xids(1),exptid(j),deblank(C{c}.errs{k}));
                           C{c}.errs(k) = [];
                           clustermodified(c) = 1;
                       end
                   end
               end
           end
           p = C{c}.probe(1);
           if p > nprobes
               nprobes = p;
           end
           if isfield(C{c},'t')
           spktimes = C{c}.t(:).*10000;
           ctimes(c) = C{c}.savetime(1);
           for t = 1:length(starts)
               n(t) = sum(spktimes > starts(t) & spktimes < ends(t));
           end
           if ~isfield(C{c},'chspk')
               if verbose
                   fprintf('E%dP%d No chspk\n',exptid(j),p);
               end
           elseif length(C{c}.chspk) ==1 && nprobes > 24
           elseif ~ismember(C{c}.probe(1),C{c}.chspk) 
               fprintf('E%dP%d Probe not in chspk %s\n',exptid(j),p,sprintf(' %d',C{c}.chspk));
           elseif max(abs(C{c}.chspk-p)) > 2
               fprintf('E%dP%d suspicious chspk %s\n',exptid(j),p,sprintf(' %d',C{c}.chspk));               
           end
           if isfield(C{c},'missingtrials')               
               missed = find(ismember(ids,C{c}.missingtrials));
           else
               missed = [];
           end
           missing = sum(n==0);
           if missing 
               mid = ids(n==0);
               mid = setdiff(mid, ids(missed)); %trials missing, not listed in Cluster.missingtrials
               if isfield(C{c},'excludetrialids') && ~isempty(C{c}.excludetrialids)
                   if iscell(C{c}.excludetrialids) && ~isempty(C{c}.excludetrialids{1})
                       fprintf('E%dP%d %d excluded trials in Cluster\n',exptid(j),p,length(C{c}.excludetrialids{1}));
                       mid = setdiff(mid, C{c}.excludetrialids{1});
                   elseif isnumeric(C{c}.excludetrialids)
                       fprintf('E%dP%d %d excluded trials in Cluster\n',exptid(j),p,length(C{c}.excludetrialids));
                       mid = setdiff(mid, C{c}.excludetrialids);
                   end
               end
               miserr(p) = sum(ismember(mid,missingids));
               miserrs{p} = mid(ismember(mid,missingids));               
               probehascell = AllV.isacell(cells, exptid(j), p);
               if probehascell
                   cstr = 'Cell';
               else
                   cstr = '';
               end
               if ~isempty(mid)
                   if isempty(oldxcl)
                       fprintf('E%dP%d%s Missing %d/%d trials (%s)\n',exptid(j),p,cstr,missing,nt,sprintf('%d,',mid));
                   else
                       if ndims(oldxcl) == 2
                           mid = setdiff(mid, oldxcl{cxid,p});
                       else
                           mid = setdiff(mid, unique(cat(2,oldxcl{cxid,p,:})));
                       end
                       if ~isempty(mid)  && ~strcmp(expname,'square.co')
                           fprintf('E%dP%d%s Missing %d/%d trials (%s)\n',exptid(j),c,cstr,missing,nt,sprintf('%d,',mid));
                       end
                   end
                   if isempty(mid) %miserrs already in excludetrials
                       miserr(p) = 0;
                   end
                   excludetrials{cid,p} = mid;
               end
           else
               miserr = [];
               excludetrials{cid,p} = [];
           end
           else
               CD = AddError(CD,'E%dP%d %s Missing .t in Details',exptid(j),p,expname);
           end
       end
       if sum(clustermodified) && exist(clusterfile)
           CX = load(clusterfile);
           cmod = find(clustermodified);
           for k = 1:length(cmod)
               CX.Clusters{cmod(k)}.errs = C{cmod(k)}.errs;
           end
           fprintf('Saving %s with new error list\n',clusterfile);
           BackupFile(clusterfile);
           save(clusterfile,'-struct','CX');
       end
       
       if size(excludetrials,1) >= cid && ~isempty(cxid) %something to check
           nconfirm = [0 0];
           if ~isfield(CD,'excludetrials')
               CD.excludetrials{length(CD.exptids),nprobes,1} = [];
           end
           missing = cat(2,excludetrials{cid,:});
           if sum(miserr > 0) > 7  && ~strcmp(expname,'square.co')             
               CD = AddError(CD, 'E%d %s %d probes lack %d/%d missing trials (%s)\n',exptid(j),expname,sum(miserr>0),max(miserr),nt,sprintf('%d,',miserrs{p}));
               modified = modified+1; 
               missed = find(miserr);
               for p = missed(:)'
                   for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = unique([CD.excludetrials{cxid,p,k} miserrs{p} excludetrials{cid,p}]);
                   end                  
               end
               modified = modified+1;
           end
           
           
           if nprobes ==96 
               [a,b] = Counts(missing);
               id = find(a > 50); %missing on 50 probes
               if ~isempty(id)
                   pid = find(CellToMat(excludetrials(cid,:)));
                   str = '';
                   FullV = FullVData{cid};
                   tid = b(id(1));
                   if checkfullv
                       for p = pid(:)'
%                           fullvname = BuildFileName(E,'fullv','probe',p);
%                           FullV = LoadFullV(fullvname);
                           if ~isempty(FullV)
                               vcheck = CheckFullV(FullV,E);
                               if isempty(vcheck.missingid)
                               end
                               if isfield(vcheck,'missingids')
                                   nconfirm(p) = sum(ismember(b(id), vcheck.missingids{1}));
                               end
                           end
                           if sum(nconfirm > 0) > 6
                               str = sprintf('Confirmed');                               
                               break;
                           end
                       end
                   else
                       vcheck = CheckFullV(FullV,E);
                       nconfirm(1) = sum(ismember(b(id), vcheck.missingids{1}));
                       if nconfirm(1)
                           str = sprintf('%d confirmed',nconfirm(1));
                       else
                           str = sprintf('FullVData OK');
                       end
                   end
                   CD = AddError(CD, 'E%d (%s) %d probes Missing %d/%d trials (%s) %s\n',exptid(j),expname,length(pid),length(id),nt,sprintf('%d,',b(id)),str);
                   %list of missing probes
                   for p = 1:size(excludetrials,2);
                       mp(p) = ismember(tid,excludetrials{cid,p});
                   end
                   mp = find(mp>0);
                   CD.errdata(end).nprobes = a(id(1));
                   CD.errdata(end).probes = mp;
                   CD.errdata(end).eid = exptid(j);
                   CD.errdata(end).ctimes = ctimes(mp);
                   modified = modified+1;

                   for p = pid(:)'
                       for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = unique([CD.excludetrials{cxid,p,k} b(id)]);
                       end
                   end
               end
           end
           if nprobes == 24
               FullV = [];
               missing = cat(2,excludetrials{cid,1:16});
               [a,b] = Counts(missing);
               id = find(a ==16); %missing on all probes
               if ~isempty(id) && ~strcmp(expname,'square.co')
                   tid = b(id(1));
               sumerrs(j,1) = length(id);
                   if checkfullv
                       fullvname = BuildFileName(E,'fullv');
                       FullV = LoadFullV(fullvname);
                       if ~isempty(FullV)
                           vcheck = CheckFullV(FullV,E);
                           if isempty(vcheck.missingid)
                           end
                           if isfield(vcheck,'missingids')
                               nconfirm(1) = sum(ismember(b(id), vcheck.missingids{1}));
                           end
                       end
                   end
                   fprintf('%s:%s ',dname,expname);
                   CD = AddError(CD, 'E%dP1-16 (%s) Missing %d/%d trials (%s) Confimed %d\n',exptid(j),expname,length(id),nt,sprintf('%d,',b(id)),nconfirm(1));
                   %list of missing probes
                   for p = 1:size(excludetrials,2);
                       mp(p) = ismember(tid,excludetrials{cid,p});
                   end
                   mp = find(mp>0);
                   CD.errdata(end).nprobes = a(id(1));
                   CD.errdata(end).probes = mp;
                   CD.errdata(end).eid = exptid(j);
                   CD.errdata(end).ctimes = ctimes(mp);
                   if checkfullv
                       fprintf('Confirmed %d missing ids\n',nconfirm(1));
                   end
                   errs{j,1} = CD.errs{end};
                   sumerrs(j,4) = nconfirm(1);
                   if checkfullv == 0 || nconfirm(1) > 0
                       for p = 1:16
                           for k = 1:size(CD.excludetrials,3)
                               CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} b(id)];
                           end
                       end
                       modified = modified+1;
                   end
               end
               if size(excludetrials,2) >= 24
               missing = cat(2,excludetrials{cid,17:24});
               [a,b] = Counts(missing);
               id = find(a == 8); %missing on all probes
               if ~isempty(id) && ~strcmp(expname,'square.co')
                   sumerrs(j,2) = length(id);
                   tid = b(id(1));
                   if checkfullv
                       fullvname = BuildFileName(E,'fullv');
                       if isempty(FullV)
                       FullV = LoadFullV(fullvname);
                       end
                       if ~isempty(FullV)
                       vcheck = CheckFullV(FullV,E);
                       if isempty(vcheck.missingid)
                       end
                       if isfield(vcheck,'missingids') && length(vcheck.missingids) > 1
                           nconfirm(2) = sum(ismember(b(id), vcheck.missingids{2}));
                       end
                       end
                   end
                   %list of missing probes
                   for p = 1:size(excludetrials,2);
                       mp(p) = ismember(tid,excludetrials{cid,p});
                   end
                   mp = find(mp>0);
                   fprintf('%s:%s ',dname,expname);
                   CD = AddError(CD, 'E%dP17-24 (%s) Missing %d/%d trials (%s)confirmed %d\n',exptid(j),expname,length(id),nt,sprintf('%d,',b(id)),nconfirm(2));
                   CD.errdata(end).nprobes = a(id(1));
                   CD.errdata(end).probes = mp;
                   CD.errdata(end).eid = exptid(j);
                   CD.errdata(end).ctimes = ctimes(mp);
                   if checkfullv
                       fprintf('Confirmed %d missing ids\n',nconfirm(2));
                   end
                   errs{j,2} = CD.errs{end};
                   if checkfullv == 0 || nconfirm(1) > 0
                       for p = 17:24
                           for k = 1:size(CD.excludetrials,3)
                               CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} b(id)];
                           end
                       end
                       modified = modified+1;
                   end
                   sumerrs(j,4) = nconfirm(2);
               end
               end
               if sum(sumerrs(j,:)) == 0
                   missing = cat(2,excludetrials{cid,:});
                   [a,b] = Counts(missing);
                   id = find(a >= nprobes/3);
                   if ~isempty(id) && ~strcmp(expname,'square.co')
                       tid = b(id(1));
                       %list of missing probes
                       for p = 1:size(excludetrials,2);
                           mp(p) = ismember(tid,excludetrials{cid,p});
                       end
                       mp = find(mp>0);
                       np = a(id(1));
                       badids = b(id);
                       trials = find(ismember(ids,b(id)));
                       if isfield(E.Trials,'st') && ~isempty(trials) && E.Header.rc ==0
%find ids that are not blanks. N.B. need this to work even if E.Trials.id has
%duplicates
                           blanks = find([E.Trials(trials).st] ==0);
                           [~,id] = setdiff(b(id),[E.Trials(trials(blanks)).id]);
                       end
                       if isempty(id)
                           fprintf('E%d %d probes empty but stim is blank (%s)\n',exptid(j),np,sprintf('%d,',badids));
                       else
                           CD = AddError(CD, 'E%d (%s) %d probes (%d-%d) Missing %d/%d trials (%s)\n',exptid(j),expname,a(id(1)),mp(1),mp(end),length(id),nt,sprintf('%d,',b(id)));
                           CD.errdata(end).nprobes = a(id(1));
                           CD.errdata(end).probes = mp;
                           CD.errdata(end).eid = exptid(j);
                           CD.errdata(end).ctimes = ctimes(mp);
                       end
                   end                   
               end

           end
       end
   end
end

if checkfullv
    ida = find(sumerrs(:,1) > 2 .* sumerrs(:,3));
    idb = find(sumerrs(:,2) > 2 .* sumerrs(:,4));
    if ~isempty(ida) || ~isempty(idb)
        fprintf('Missing in ClusterDetails but not in FullV:\n');
        for j = 1:length(ida)
            fprintf('%d/%d: %s',sumerrs(ida(j),[1 3]),errs{ida(j),1});
        end
        for j = 1:length(idb)
            fprintf('%d/%d: %s',sumerrs(idb(j),[2 4]),errs{idb(j),2});
        end
    end
end

if modified == 0 && fixlist
    fprintf('No errors affected entire probe groups');
    [good, cells.CellDetails] = CheckErrList(cells.CellDetails);
    if sum(good ==0) %some need to be removed
        BackupFile(cellname);
        save(cellname,'-struct','cells');
    end
end

if fixlist
    if modified || ~isfield(CD,'checkedtimes')
        CD.checkedtimes(exptid) = 1;
        for j = 1:size(CD.excludetrials,1)
        for k = 1:size(CD.excludetrials,2)
        for m = 1:size(CD.excludetrials,3)
            CD.excludetrials{j,k,m} = unique(CD.excludetrials{j,k,m}); 
        end
        end
        end
        cells.CellDetails = CD;
        BackupFile(cellname);
        save(cellname,'-struct','cells');
    end
end
CD.name = dname;



function result = CheckExptClusterDone(dname, varargin)

undo = 0;
fix = 0;
check = '';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'details',4)
        check = 'details';
    elseif strncmpi(varargin{j},'fix',3)
        fix =1;
    elseif strncmpi(varargin{j},'trials',4)
        check = 'trials';
    elseif strncmpi(varargin{j},'undo',4)
        undo = 1;
    end
    j = j+1;
end
cellfile = [dname '/CellList.mat'];
if ~exist(cellfile)
    result.done = 0;
    result.ok = 0;
    result.err = 'no CellList';
    return;
end
d = dir(cellfile);
latest = d.datenum;
backups = mydir([dname '/backup/CellList*.mat']);
if isempty(backups)
    result.done = 0;
    result.ok = 1;
    return;
else
    [a, b] = sort([backups.datenum]);
    backups = backups(b);
    oldfile = backups(end).name;
end


gotbackup = 0;
 for j = 1:length(backups)
     O = load(backups(j).name);
     result.version = O.CellDetails.version;
     if ~isfield(O.CellDetails,'checkedtimes')
         fprintf('%s is safe backup\n',backups(j).name);
         oldfile = backups(j).name;
         gotbackup = 1;
     end
 end
 if ~gotbackup
     backups = mydir([dname '/CellList*.mat']);
    [a, b] = sort([backups.datenum]);
    backups = backups(b);
    for j = 1:length(backups)
        O = load(backups(j).name);
        if isfield(O.CellDetails,'version')
            result.version = O.CellDetails.version;
        else
            result.version = 0;
        end
        if ~isfield(O.CellDetails,'checkedtimes')
            fprintf('%s is safe backup\n',backups(j).name);
            oldfile = backups(j).name;
            gotbackup = 1;
        end
    end
 end
 O = load(oldfile);
   


if strcmp(check,'trials')
    C = load(cellfile);
    for j = 1:length(C.CellDetails.excludetrials(:))
        C.CellDetails.excludetrials{j} = C.CellDetails.excludetrials{j}(:)'; 
    end

    if isfield(C.CellDetails,'checkedtimes')
        fprintf('%s has been fixed\n',dname);
    end
    Expts = ReadExptDir(dname);
    for j = 1:length(Expts)
        Tn{j} = [Expts{j}.Trials.id];
    end
    expid = GetExptNumber(Expts);
    result.idmatch = [];
    fixerr = 0;
    for j = 1:length(C.CellDetails.exptids)
        e = C.CellDetails.exptids(j);
        k = find(expid == e);
        if ~isempty(k)
            id = find(~ismember(C.CellDetails.trialids{j},[Expts{k}.Trials.id]));
            if ~isempty(id)
                if e <= length(Expts)
                    id = find(~ismember(C.CellDetails.trialids{j},[Expts{e}.Trials.id]));
                    if isempty(id)
                        result.idmatch(end+1,:) = [e j];
                        fprintf('%s Trial ids %d(%d) is for Expt %d\n',dname,e,j,e);
                    else
                        fprintf('%s Trial ids wrong in Expt %d\n',dname,e);
                        fixerr = fixerr+1;
                    end
                else
                    fprintf('%s Trial ids wrong in Expt %d\n',dname,e);
                        fixerr = fixerr+1;
                end
            end
            result.iderr(j) = length(id);
            result.idlist(j,:) = [length(C.CellDetails.trialids{j}) length([Expts{k}.Trials.id])];
        else
            fprintf('Expt %d missing from Exts in %s\n',e,dname);
        end
    end
    for j = 1:size(C.CellDetails.excludetrials,1)
        exclid = cat(2,C.CellDetails.excludetrials{j,:,:});
        if ~isempty(exclid)
            for k = 1:length(Tn)
                match(k) = sum(ismember(exclid,Tn{k}));
            end       
            id = find(match > 0);
            if sum(match) == 0
                fprintf('No match for row %d',j);
            elseif length(id) > 1
                cprintf('red','Exclusions row %d matches %d expts%s\n',j,length(id),sprintf(' %d',expid(id)));                
            elseif j <= length(C.CellDetails.exptids)
                fprintf('Exclusions row %d is Expt %d (CD) %d (Expts) \n',j,C.CellDetails.exptids(j),expid(id));
                if expid(id) ~= C.CellDetails.exptids(j)
                    e = find(C.CellDetails.exptids == expid(id));
                    exclidb = cat(2,C.CellDetails.excludetrials{e,:,:});
                    if isempty(setdiff(exclid,exclidb))
                        cprintf('red','%s: ExcludeTrialList line %d duplicates %d\n',dname,j,e)
                    else
                        cprintf('red','%s: ExcludeTrialList wrong line\n',dname)
                    end
                end
            else
                cprintf('red','Exclusions row %d is Expt %d longer than exptid\n',j,k);
            end
        elseif j > length(C.CellDetails.exptids)
            cprintf('red','Exclusions row %d/%d Has %d exclusions and is past end (%s)\n',j,size(C.CellDetails.excludetrials,1),length(exclid),dname);
            fixerr = fixerr+1;
        end
    end
    if fix && fixerr
        if size(C.CellDetails.excludetrials,1) > length(C.CellDetails.exptids)
            j = length(C.CellDetails.exptids)+1;
            exclid = cat(2,C.CellDetails.excludetrials{j:end,:,:});
            if isempty(exclid)
                C.CellDetails.excludetrials = C.CellDetails.excludetrials(1:j-1,:,:);
            end
        end
        for j =1:length(C.CellDetails.exptids)
            e = C.CellDetails.exptids(j);
            k = find(expid == e);
            if isempty(k)
                C.CellDetails.trialids{j} = [];
            else
                C.CellDetails.trialids{j} = Tn{k};
            end
        end
        BackupFile(cellfile,'show');
        save(cellfile,'-struct','C');
    end
    return;
elseif strcmp(check,'details')
    result.version = 0;
    C = load(cellfile);
    if isfield(C.CellDetails,'checkedtimes')
        fprintf('%s has been fixed\n',dname);
    end
    for j = 1:length(backups)
        O = load(backups(j).name);
        result.version = O.CellDetails.version;
        if ~isfield(O.CellDetails,'checkedtime')
            fprintf('%s is safe backup\n',backups(j).name);
        CD = O.CellDetails;
        result.size(1) = size(CD.excludetrials,1);
        result.size(2) = length(CD.exptids);
        result.size(3) = max(CD.exptids);
        missing = setdiff(1:result.size(1),CD.exptids);
        if diff(result.size(1:2)) < 0            
            fprintf('Size mismatch. Version %.2f',CD.version);
            exids = find(sum(sum(CellToMat(CD.excludetrials),3),2));
            if max(exids) > missing(1)
                fprintf('Exclusions in rows %s, missing expts %s\n',sprintf('%d ',exids),sprintf('%d ',missing));
                for e = exids(:)'
                    ids = cat(2,CD.excludetrials{e,:,:});
                    ids = unique(ids);
                end
            end
        elseif length(CD.exptids) > result.size(1)
            fprintf('%s Exptids is missing some exclude lists\n',dname);          
        elseif max(CD.exptids) > result.size(1)
            fprintf('%s Exptids skips some expts\n',dname);          
        end
        result.version = CD.version;
        result.exptids = CD.exptids;
        result.trialids = CD.trialids;
        result.excludetrials = CD.excludetrials;
        break;
        end
    end
    return;
end
if exist(cellfile) && exist(oldfile)


    fprintf('Checking %s\n',cellfile);
    C = load(cellfile);
    O = load(oldfile);
    if isfield(C.CellDetails,'checkedtimes')
        result.done = 1;
    else
        result.done = 0;
    end
    xcl = C.CellDetails.excludetrials;
    nerr = 0;
    if size(xcl,1) > length(C.CellDetails.exptids)
        fprintf('exclude list longer that exptids\n');
    end
    expids = C.CellDetails.exptids;
    for j = 1:size(xcl,1)
        ids = cat(2,xcl{j,:,:});
        ids = unique(ids);
        oldid = unique(cat(2,O.CellDetails.excludetrials{j,:,:}));
        e = find(expids == j);
        if ~isempty(e)
        errs =setdiff(ids,C.CellDetails.trialids{j});
        if ~isempty(errs)
            fprintf('%s not in id list for expt %d\n',sprintf('%d ',errs),O.CellDetails.exptids(e));
            if undo
                for k = 1:size(xcl,2)
                    for m = 1:size(xcl,3)
                        C.CellDetails.excludetrials{j,k,m} = setdiff(xcl{j,k,m},errs);
                    end
                end
            end
            nerr = nerr + length(errs);
        end
end
    end
    if undo
        if nerr && undo == 2
            changecount = zeros(size(xcl,1),size(xcl,2));
            for j = 1:size(xcl,1)
                for k = 1:size(xcl,2)
                    for m = 1:size(xcl,3)
                        changed{j,k,m} = setdiff(O.CellDetails.excludetrials{j,k,m}, C.CellDetails.excludetrials{j,k,m});
                        changed{j,k,m} = changed{j,k,m}(:)';
                        changecount(j,k) = changecount(j,k) + length(changed{j,k,m});
                        goodoldcount(j,k) = sum(ismember(oldid,C.CellDetails.trialids{j}));
                    end
                    oldid = unique(cat(2,changed{j,k,:}));
                    oldcount(j,k) = sum(ismember(oldid,C.CellDetails.trialids{j}));  %trials really gond missing
                end
            end
            C.CellDetails.fixdate = now;
            nc = length(cat(2,changed{:}));
            if nc == 0
                BackupFile(cellfile);
                save(cellfile,'-struct','C');
            else
                fprintf('%d changed elements\n',nc);
                if sum(changecount(:)) == nerr || sum(oldcount(:) ==0) %removede only bad id listings
                    fprintf('Saving %s\n',cellfile);
                    BackupFile(cellfile);
                    save(cellfile,'-struct','C');
                else
                    cprintf('red','Change count does not match\n');
                end
            end
        elseif undo == 1 && isfield(C.CellDetails,'checkedtimes') && gotbackup
            if sum(C.CellDetails.exptids ~= O.CellDetails.exptids) ==0
                C.CellDetails = O.CellDetails;
                [good, C.CellDetails] = CheckErrList(O.CellDetails);
                if sum(good ==0) %some old errors removed
                    BackupFile(cellfile);
                    save(cellfile,'-struct','C');
                end
            else
                fprintf('exptids mismatch\n');
            end
        elseif ~gotbackup
            fprintf('No backup for %s\n',cellfile);
        end
    end
    
end


function fixed = ShowResult(CD, fix)

fixed = {};
fixexpt = [];
if isempty(CD)
    return;
end
eids = [];
name = CD.name;
if isfield(CD,'errs')
for j = 1:length(CD.errs)
    estr = '';
    if CD.errdata(j).time > CD.checktime-0.1
            if isfield(CD.errdata,'eid') && ~isempty(CD.errdata(j).eid)
                eid = CD.errdata(j).eid;
            else
                [eid, p] = GetExptNumber(CD.errs{j});
            end
                if isfield(CD,'exptnames') && eid > 0
                    estr = [' ' CD.exptnames{eid}];
                end
        fprintf('%s%s:%s\n',name,estr,deblank(CD.errs{j}));
        eids(j) = eid;
        if fix && ~isempty(CD.errdata(j).nprobes) && CD.errdata(j).nprobes >= fix
            fixexpt(eid) = 1;
        end
    end
end
end

eids = eids(eids>0);
[a,b] = Counts(eids);
if max(a) > 1
fprintf('Expts With problems on > 1 probe: %s\n',sprintf(' %d',b(a>1)));
end

if fix
    expts = find(fixexpt);
    for j = 1:length(expts)
        fullv = sprintf('%s/Expt%dFullV.mat',CD.name,expts(j));
       fixed{j} = AllVPcs(fullv,'reclassifyall','checklast','savespikesifsafe');
    end
end

function [good, fixed] = CheckErrList(CD, fix)
good = [];
if ~isfield(CD,'errs')
    fixed = CD;
    return;
end
for j = 1:length(CD.errs)
    if regexp(CD.errs{j},'E[0-9]+P[0-9].* Missing [0-9]+/[0-9+]')
        good(j) = 0;
    else
        good(j) = 1;
    end
end


fixed = CD;
if sum(good==0)
    fixed.errs =fixed.errs(find(good));
    fixed.errdata =fixed.errdata(find(good));
end


