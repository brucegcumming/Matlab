function fixed = RemakeClusters(dname, varargin)
% fixed - RemakeClusters(D, varargin) Calls AllVPcs for Clusters that need
% remaking because of detecable errors in the ClusterTimes/Details files
% call on its own output to summarize results
% RemakeClusters(fixed)  displays results where events marked as a cell
% would be affected.
% RemakeClusters(fixed,'showall')  displays results of all calls to AllVPcs

eid = [];
Expts = {};
CD  = [];
fixlist = 0;
check = '';
checkfullv = 0;
checkonly = 0;
verbose = 0;
rebuildemptytrials = 0;
fixed = {};
saveargs = {};
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        X = varargin{j};
        [a,b] = CellToMat(X,'Header');
        if sum(b.found) > 0
            Expts = X;
        end
    elseif isstruct(varargin{j})
        if isfield(varargin{j},'trialids')
            CD = varargin{j};
        end
    elseif strncmpi(varargin{j},'checkfullv',10)
        checkfullv = 1;
    elseif strncmpi(varargin{j},'checkdetails',10)
        check = 'details';
    elseif strncmpi(varargin{j},'check',5)
        checkonly = 1;
    elseif strncmpi(varargin{j},'trials',5)
        check = 'trials';
    elseif strncmpi(varargin{j},'savespike',5)
        saveargs = {saveargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        eid = varargin{j};
    elseif strncmpi(varargin{j},'fix',3)
        fixlist = 1;
    end
    j = j+1;
end

if ~isempty(Expts)
    exptid = GetExptNumber(Expts);
end

if iscellstr(dname)
    for j = 1:length(dname)
        if strcmp(check,'done')
        elseif sum(strcmp(check,{'details' 'trials'}))
        else
            fixed{j} = RemakeClusters(dname{j},varargin{:});
        end
    end
    return;
elseif iscell(dname) %set of results
   ShowResults(dname, varargin{:});
elseif isdir(dname)
    if isempty(Expts)
        Expts = ReadExptDir(dname);
    end
    exptid = GetExptNumber(Expts);
    if ~isempty(eid)
        id = find(ismember(exptid,eid));
        
        Expts = Expts(exptid(id));
    else
        id= [];
    end
    [Clusters, FullVData] = LoadCluster(dname, GetExptNumber(Expts),'array');
    cexpts = prctile(CellToMat(Clusters,'exptno')',50);
    cellname = [dname '/CellList.mat'];
    if exist(cellname)
        cells = load(cellname);
        CD = cells.CellDetails;
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
if isfield(CD,'excludetrials')
    oldxcl = CD.excludetrials;
else
    oldxcl = {};
end
nfix = 0;
fixed = {};
for j = 1:length(Expts);
    cid = find(cexpts == exptid(j));
    sumerrs(j,1:4) = 0;
    if ~isempty(Expts{j}) && ~isempty(cid)
        
        C = Clusters{cid};
        E = Expts{j};
        fullvname = BuildFileName(E,'fullv');
        if isfield(CD,'exptids')
            cxid = find(CD.exptids == exptid(j));
        else
            cxid = [];
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
        for c = 1:length(C)
            spktimes = C{c}.t(:).*10000;
            p = C{c}.probe(1);
            if p > nprobes
                nprobes = p;
            end
            for t = 1:length(starts)
                n(t) = sum(spktimes > starts(t) & spktimes < ends(t));
            end
            if ~isfield(C{c},'chspk')
                if verbose
                    fprintf('E%dP%d No chspk\n',exptid(j),p);
                end
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
            if isfield(C{c},'needed') && C{c}.needed == 1
                if checkonly
                    fixed{end+1}.name = dname;
                    fixed{end}.probe = p;
                    fixed{end}.expt = exptid(j);
                    fixed{end}.missing = [];
                    fixed{end}.needed = 1;
                    cprintf('blue', 'Will need %s E%dP%d\n',dname,exptid(j),p);
                else
                    fixed{end+1} = AllVPcs(fullvname,'tchan',p,'reclassify',saveargs{:});
                    fixed{end}.needed = 1;
                end
            elseif missing
                mid = ids(n==0);
                mid = setdiff(mid, ids(missed)); %trials missing, not listed in Cluster.missingtrials
                if isfield(C{c},'excludetrialids') && ~isempty(C{c}.excludetrialids)
                    fprintf('E%dP%d %d excluded trials in Cluster\n',exptid(j),p,length(C{c}.excludetrialids));
                    mid = setdiff(mid, C{c}.excludetrialids);
                    if ~isempty(mid) && checkfullv
                        fmiss = CheckFullV(FullVData{j},Expts{j});
                        if ~isempty(fmiss.missingids{1})
                            mid = setdiff(mid,fmiss.missingids{1});
                        end
                    end
                    if ~isempty(mid) %trials with no evetnts, but are in FullV.  Do we need to run AllVPcs here?
                        if checkonly || rebuildemptytrials ==0
                            fixed{end+1}.name = dname;
                            fixed{end}.probe = p;
                            fixed{end}.expt = exptid(j);
                            fixed{end}.missing = mid;
                        else
                        fixed{end+1} = AllVPcs(fullvname,'tchan',p,'reclassify',saveargs{:});
                        fixed{end}.missing = mid;
                        end
                    end
                end
            else
            end
        end
    end
end

function ShowResults(fix, varargin)
showall = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'showall',6)
        showall = 1;
    end
    j = j+1;
end
for j = 1:length(fix)
    if iscell(fix{j})
        R = fix{j};
        for k = 1:length(R);
            if isfield(R{k},'cluster')
            C = R{k}.cluster;
            name = regexprep(R{k}.cluster.spkfile,'Spikes/.*','');
            eid = GetExptNumber(C);
            s = sprintf('%s E%dP%d: %d %d',name,eid,C.probe(1),R{k}.unsafetosave,R{k}.cell);
            if R{k}.cell > 0 && R{k}.unsafetosave ~= 0
                cprintf('red','%s Spks %d %d %d %d\n',s,R{k}.matchcounts(1:4));
            elseif showall
                fprintf('%s\n',s);
            end
            elseif isfield(R{k},'name') && showall
                if isfield(R{k},'needed')
                    cprintf('blue','%s needed %d\n',R{k}.name,R{k}.needed);
                else
                    fprintf('%s not needed E%dP%d %d empty trials\n',R{k}.name,R{k}.expt,R{k}.probe,length(R{k}.missing));
                end
            end
        end
    end
end