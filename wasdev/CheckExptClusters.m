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
        CD{j} = CheckExptClusters(dname{j},varargin{:});
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
    if iscell(dname{1})
        Clusters{1} = dname;
        cexpts = GetExptNumber(Clusters{1});
    elseif isstruct(dname{1}) %combine AllClsuters struct
        cexpts = GetExptNumber(dname);        
        for j = 1:length(dname)
            for k = 1:length(dname{j})
            Clusters{j}{k} = dname{j}(k);
            if isfield(dname{j},'t')
                Clusters{j}{k}.t = dname{j}(k).t;
            elseif isfield(dname{j},'times')
                Clusters{j}{k}.t = dname{j}(k).times;
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
    Clusters = LoadCluster(dname, GetExptNumber(Expts));
    if length(id) == 1
        AllClusters{1} = Clusters;
        Clusters = AllClusters;
    end
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
if isfield(CD,'excludetrials')
    oldxcl = CD.excludetrials;
else
    oldxcl = {};
end
for j = 1:length(Expts);
    cid = find(cexpts == exptid(j));
   if ~isempty(Expts{j}) && ~isempty(cid)
       C = Clusters{cid};
       E = Expts{j};
       if isfield(CD,'exptids')
           cxid = find(CD.exptids == exptid(j));
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
               missed{c} = C{c}.missingtrials;
           end
       end
       missingids = unique(cat(1,missed{:}));
       if ~isempty(missingids)
               fprintf('Expt%d Missing Trials %s\n',exptid(j),sprintf('%d ',missingids));
       end
       misserr = 0;
       for c = 1:length(missed)
           if length(missed{c}) ~= length(missingids)
               misserr(c) = 1;
           end
       end
       if sum(misserr)
           bid = find(misserr);
           fprintf('Probes %s incomplete missing Trials (%s)\n',sprintf(' %d',bid),sprintf(' %d',missingids))
       end
       for c = 1:length(C)
           spktimes = C{c}.t(:).*10000;
           p = C{c}.probe(1);
           if p > nprobes
               nprobes = p;
           end
           for t = 1:length(starts)
               n(t) = sum(spktimes > starts(t) & spktimes < ends(t));
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
               miserr(p) = sum(ismember(mid,missingids));
               miserrs{p} = mid(ismember(mid,missingids));
               if ~isempty(mid)
                   if isempty(oldxcl)
                       fprintf('E%dP%d Missing %d/%d trials (%s)\n',exptid(j),p,missing,nt,sprintf('%d,',mid));
                   else
                       if ndims(oldxcl) == 2
                           mid = setdiff(mid, oldxcl{cxid,p});
                       else
                           mid = setdiff(mid, unique(cat(2,oldxcl{cxid,p,:})));
                       end
                       if ~isempty(mid)
                           fprintf('E%dP%d Missing %d/%d trials (%s)\n',exptid(j),c,missing,nt,sprintf('%d,',mid));
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
       end
       if size(excludetrials,1) >= cid %something to check
           missing = cat(2,excludetrials{cid,:});
           if sum(miserr)
               CD = AddError(CD, 'E%d %d probes lack %d/%d missing trials (%s)\n',exptid(j),sum(miserr>0),max(miserr),nt,sprintf('%d,',id));
               modified = modified+1; 
               missed = find(miserr);
               for p = missed(:)'
                   for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} miserrs{p} excludetrials{cid,p}];
                   end                  
               end
               modified = modified+1;
           end
           if nprobes ==96 
               [a,b] = Counts(missing);
               id = find(a > 50); %missing on 50 probes
               if ~isempty(id)
                   pid = find(CellToMat(excludetrials(cid,:)));
                   CD = AddError(CD, 'E%d %d probes Missing %d/%d trials (%s)\n',exptid(j),length(pid),length(id),nt,sprintf('%d,',id));
                   modified = modified+1;
                   for p = pid(:)'
                       for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} b(id)];
                       end
                   end
               end
           end
           if nprobes == 24
               missing = cat(2,excludetrials{cid,1:16});
               [a,b] = Counts(missing);
               id = find(a ==16); %missing on all probes
               if ~isempty(id)
                   CD = AddError(CD, 'E%dP1-16 Missing %d/%d trials (%s)\n',exptid(j),length(id),nt,sprintf('%d,',id));
                   modified = modified+1;
                   for p = 1:16
                       for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} id];
                       end
                   end
               end
               missing = cat(2,excludetrials{cid,17:24});
               [a,b] = Counts(missing);
               id = find(a == 8); %missing on all probes
               if ~isempty(id)
                   CD = AddError(CD, 'E%dP1-16 Missing %d/%d trials (%s)\n',exptid(j),length(id),nt,sprintf('%d,',id));
                   for p = 17:24
                       for k = 1:size(CD.excludetrials,3)
                       CD.excludetrials{cxid,p,k} = [CD.excludetrials{cxid,p,k} id];
                       end
                   end
                   modified = modified+1;
               end
           end
       end
   end
end

if modified == 0 && fixlist
    fprintf('No errors affected entire probe groups');
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