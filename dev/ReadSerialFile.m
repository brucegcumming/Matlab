function [result, Expt, txt] = ReadSerialFile(name, varargin)
%[result, Expt, txt] = ReadSerialFile(name, 'readexpt')
% read a serial output file made by binoc
% Expt is and Expt structure with psychophysical results.
% If the file contains more than one psychophysical trial type, they ill
% all be included by default, which will break many scripts. To limit Expt
% to only one experiment type use one of the following:
%
%[result, Expt, txt] = ReadSerialFile(name, 'readexpt','exptype','ORBW')
%[result, Expt, txt] = ReadSerialFile(name, 'readexpt','exptype','DCOR')
%
% result lists expt types, and the values of 'id' for each trial 
% txt returns the contents of the file as a matlab character array
%result = ReadSerialFile(name,'chkFr') adds the values for Fr
%result = ReadSerialFile(name,'chkrpts') checks seed repeats for ob=130
%
%
maxlen = 12000;
checkfx =0;
checkFr = 0;
checkseq = 0;
checkrpts = 0;
forcelen = 0;
checkTimes = 1;
chectrials = 0;
readexpts = 0;
checkepos=0;
addfields = {'rw' 'mD' 'Dc'};
Expt = [];
exptlist = [];
mkufl = 0;
exptype = ''; %default
verbose = 1;
checkexpts = 0;

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j}, 'Trials')
        Expt = varargin{j};
    elseif strncmpi(varargin{j},'chkFr',5)
        checkFr= 1;
    elseif strncmpi(varargin{j},'chkexpts',4)
        checkexpts = 1;
    elseif strncmpi(varargin{j},'chkepos',4)
        checkepos = 1;
    elseif strncmpi(varargin{j},'chkfix',4)
        checkfx = 1;
    elseif strncmpi(varargin{j},'chkseq',6)
        checkseq = 1;
    elseif strncmpi(varargin{j},'chkrpt',6)
        checkrpts = 1;
    elseif strncmpi(varargin{j},'exptype',4)
        j = j+1;
        exptype = varargin{j};
    elseif strncmpi(varargin{j},'field',4)
        j = j+1;
        addfields = {addfields{:} varargin{j}};
    elseif strncmpi(varargin{j},'mkufl',5)
        mkufl = 1;
    elseif strncmpi(varargin{j},'forcelen',4)
        j = j+1;
        forcelen = varargin{j};
    elseif strncmpi(varargin{j},'maxl',4)
        j = j+1;
        maxlen = varargin{j};
    elseif strncmpi(varargin{j},'quiet',5)
        verbose = 0;
    elseif strncmpi(varargin{j},'readexpt',6)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            exptlist = varargin{j};
            checkrpts = 1;
        end
        readexpts = 1;
    elseif strncmpi(varargin{j},'trials',5)
        checktrials = 1;
    end
    j = j+1;
end

if iscell(name)
    if exist(name{1},'file')
        Expt = {};
        for j = 1:length(name)
            fprintf('Reading %s\n',name{j});
            [result{j}, Expt{j}] = ReadSerialFile(name{j},varargin{:});
        end
    else
        if checkexpts
            result = CheckExptLoops(name, checkexpts);
            return;
        end
        
    end
    return;
elseif ischar(name)
fid = fopen(name,'r');
Header.name = name;
if fid < 1
    result = [];
    fprintf('Can''t read %s\n',name);
    return;
end
tic;
result.filename = name;




a = textscan(fid,'%s','delimiter','\n','bufsize',maxlen);
fclose(fid);


if checkexpts
    result = CheckExptLoops(a{1});
    txt = a{1};
    return;
end

txt = char(a{1});
if forcelen
    txt = txt(:,1:forcelen);
end
elseif iscellstr(name)
    txt = name;
end

id = strmatch('#Not in',txt);
%gid = setdiff(1:size(txt,1),id);
%txt = txt(gid,:);
%get rid of lines saying Can't read - they mess up line conting
xid = strmatch('Can''t',txt);
gid = setdiff(1:size(txt,1),xid);
txt = txt(gid,:);

gid = strmatch('RG',txt);
wid = strmatch('RW',txt);
bid = strmatch('RB',txt);
lid = strmatch('RL',txt);
fid = strmatch('RF',txt);
rwid = strmatch('rw',txt);
idid = strmatch('id',txt);
dxid = strmatch('dx',txt);
stid = strmatch('O 5',txt);
psyid = strmatch('psyv=',txt);
endid = strmatch('O 3',txt);
exid = strmatch('Stimulus',txt);
types(gid) = 1;
types(wid) = 2;
types(bid) = 3;
types = types(find(types > 0));
pid = union(gid,wid);
pid = union(pid,bid);
pid = union(pid,lid);
pid = union(pid,fid);
scores = zeros(size(pid));
scores(ismember(pid, gid)) =1;
scores(ismember(pid, wid)) = -1;
scores(ismember(pid, lid)) = -2;
scores(ismember(pid, fid)) = 0;
for j = 1:length(pid)
    if strfind(txt(pid(j),:),'Corloop')
        scores(j) = scores(j)-10;
    end
end



if mkufl
    result.rfstr = MkUfl(name,txt);
end

if readexpts
    for j = 1:length(idid)
        ids(j) = sscanf(txt(idid(j),:),'id%d');
    end
    nfid = strmatch('nf',txt);
    for j = length(nfid):-1:1
        nfs(j) = sscanf(txt(nfid(j),:),'nf%d');
    end
    Stimvals.e3 = 'e0';
    Header.rc = 0;
    sfid = strmatch('sf',txt);
    wiid = strmatch('wi',txt);
    seid = strmatch('se',txt);
    orid = strmatch('or',txt);
    bhid = strmatch('bh',txt);
    bdid = strmatch('bd',txt);
    jxid = strmatch('jx',txt);
    soid = strmatch('so',txt);
    siid = strmatch('#id',txt);
    sqid = strmatch('mtse',txt);
    eiid = strmatch('mtei',txt);
    mtrSid = strmatch('mtrS',txt);
    trid = strmatch('tr',txt);
    etid = strmatch('et',txt);
    e2id = strmatch('e2',txt);
    mdid = strmatch('mD',txt);
    exid = strmatch('exvals',txt);
    for j = 1:length(addfields)
        xfid{j} = strmatch(addfields{j},txt);
    end
    lasteiid = 0;
        otherrcs = {};
        trialmode = 0;
    if isempty(sqid)
        sqid = strmatch('dx:',txt);
        if ~isempty(sqid)
            trialmode = 1;
            dxid = setdiff(dxid,sqid);
            ceid = strmatch('ce:',txt);
            otherrcs(1).id = ceid;
            otherrcs(1).type = 'ceseq';
        end
    end

    if length(exptlist)
        remid = strmatch('Remaining',txt);
%        exid = strmatch('Expt',txt);
        exendid = strmatch('Run end',txt);
        allsq = [];
        allrem = [];
        for j = 1:length(exptlist)
            k = exptlist(j);
            id = find(sqid > exid(k) & sqid < exendid(k));
            allsq = [allsq id'];
            id = find(remid > exid(k) & remid < exendid(k));
            allrem = [allrem id'];
            if length(exptlist) == 1
                id = find(seid > exid(k) & seid < exendid(k));
                seid = seid(id);
                id = find(remid > exid(k) & remid < exendid(k));
                remid = remid(id);
            end
            if checkrpts
                for k = 1:length(remid)
                    id = find(pid < remid(k));
%                    fprintf('%s %s\n',txt{pid(id(end))}(1:2),txt{remid(k)});
                end
            end
        end
        sqid = sqid(allsq);
    end
imid = strmatch('imve',txt);
badimid = strmatch('imve 0.00',txt);
imid = setdiff(imid,badimid);
%RB lines come before imse
%R[WG] com after imse.
% so first find the line that comes
x = 0;

addid = [];
for j = 1:length(mtrSid)
    aid = find(stid < mtrSid(j));
    if ~isempty(aid) && aid(end) < length(stid)
    id = find(sqid > stid(aid(end)) & sqid < stid(aid(end)+1));
    if isempty(id) %
        addid = [addid mtrSid(j)];
    end
    end
end
sqid = union(sqid,addid);

if isempty(sqid)
    sqid = pid+1;
end

for j = 1:length(sqid)
    id = find(nfid < sqid(j));
    if isempty(id)
        id = 1;
    end
    nframes(j) = nfs(id(end))-1;
    if trialmode == 1 && ~isempty(psyid)%sqid lines are before O 5 for current trial
        aid = find(pid < sqid(j));
        bid = find(psyid > sqid(j));
        if isempty(aid)
            startid = sqid(j) - 10;
        else
            startid = pid(aid(end));
        end
        if isempty(bid)
            nextid  = size(txt,1)
        else
            nextid = psyid(bid(1))-1;
        end
    else
        aid = find(stid < sqid(j));
        bid = find(endid < sqid(j));
        if ~isempty(aid) && aid(end) < length(stid)
            startid = stid(aid(end));
            nextid = stid(aid(end)+1);
        else
            startid = sqid(j)-50;
            nextid = sqid(j)+12;
        end
    end
    line = txt(sqid(j),:);
    if strncmp(line,'ce:',3)
        nc=5;
        rcvar = 'ceseq';
    elseif strncmp(line,'dx:',3)
        nc = 5;
        rcvar = 'dxseq';
    else
        rcvar = 'Seedseq';
        nc = 6;
    end
    a = sscanf(line(nc:end),'%f ');
    nrc(j) = length(a);
    Trials(j).nFr(1) = length(a);
    Trials(j).nFr(2) = length(a);
    
    
    id = find(idid < sqid(j));
    seqid(j) = ids(id(end));
%Response string is after seedseq, by 10 lines on completed trials
% but before seqid on bad trials
%can be much later on Late resp trials
    id = find(pid < nextid & pid > startid);
    if isempty(id)
        id = find(idid < sqid(j));
        if verbose
            fprintf('Missing Response String for line %d %s\n',sqid(j),txt(idid(id(end)),:));
        end
        id = find(pid < sqid(j));
        if isempty(id)
            d = 100;
            t = 1;
        else
            d = 100;
            t = id(end);
        end
    else
        if pid(id(1)) > sqid(j) && id(1) > 1
            nprelin(j) = sqid(j)-pid(id(1)-1);
            xid = find(imid < sqid(j));
            if ~isempty(xid) && imid(xid(end)) < pid(id(1)-1)
 %               id(1) = id(1)-1;
 %               nprelin(j) = -nprelin(j);
            end
        end
        t = id(1);
        d = pid(t) - sqid(j);
    end
    if d > 20 & pid(t) ~= -2 && length(id) > 2%only late trials should be this late
        if verbose
            fprintf('%s Suspicious Trial (no sequence) id %d\n',name,seqid(j));
        end
%        t = id(end-1);
    end    
        
    rcvars{j} = rcvar;
    Trials(j).(rcvar) = a;
    Trials(j).nrc = nrc(j);
    seedseq{j} = a;
    for k = 1:length(otherrcs)
        oid = otherrcs(k).id;
        id = find(oid > startid & oid < nextid);
        line = txt(oid(id(end)),:);
        nc = 6;
        a = sscanf(line(nc:end),'%f ');
        Trials(j).(otherrcs(k).type) = a;
    end
    
    a = textscan(txt(pid(t),:),'%2s %[^=]=%f %[^=]=%f %[^=]=%f');
    %only set stimvals.et, e2 if its a psych trial
    if sum(strcmp(a{1}{1},{'RG' 'RW'}))
        Stimvals.et = a{2}{1};
        Stimvals.e2 = a{4}{1};
    end
    if strcmp(a{2},'pR')
        Trials(j).pR = a{3}(1);
    end
    if strcmp(a{4},'sn')
        snfield = 5;
    else
        snfield = 7;
    end
    Trials(j).(a{2}{1}) = a{3}(1);
    Trials(j).(a{4}{1}) = a{5}(1);
    id = find(exid < sqid(j));
    if ~isempty(id)
        if length(id) > 1
        mx = sscanf(txt(exid(id(end-1)),:),'exvals %f %f %f %f');
        else
            mx = [];
        end
        x = sscanf(txt(exid(id(end)),:),'exvals%f %f %f %f');
        if length(mx) > length(x)
            x(length(mx)) = mx(end);
        end
        Trials(j).exvals = x;
    end
%    id = find(dxid < sqid(j));
%    if ~isempty(id)
%        Trials(j).dx = sscanf(txt(dxid(id(end)),:),'dx%f');
%    end
    
    if scores(t) == 0 || scores(t) == -2 || length(seedseq{j}) > 250
        Trials(j).result = -1;
        Trials(j).RespDir = 0;
    elseif scores(t) < -5
        Trials(j).RespDir = 0;
        Trials(j).result = scores(t);
    elseif ~isempty(a{snfield}) %sign is in the line
        Trials(j).result = 1;
        if a{snfield}(1) == 0
            Trials(j).correct = scores(t);
            id = find(psyid < sqid(j));
            if ~isempty(psyid)
                psyv = sscanf(txt(psyid(id(end)),:),'psyv=%f');
                Trials(j).RespDir = scores(t) .* sign(psyv);
                Trials(j).rwdir = sign(psyv);
                Trials(j).psyv = psyv;
            end
        else
            Trials(j).RespDir = scores(t) .* sign(a{snfield}(1)-0.5);
            Trials(j).rwdir = sign(a{snfield});
        end
    else
        Trials(j).result = -2;  %%for now, if rwsign not known, cant get RespDir;
        Trials(j).RespDir = 0;
    end
    

    
    id = find(mtrSid <= nextid);
    if ~isempty(id) && mtrSid(id(end)) > startid
        Trials(j).seqoffset = sqid(j)-mtrSid(id(end));
        Trials(j).seline = sqid(j);
        line = txt(mtrSid(id(end)),:);
        a = sscanf(line(6:end),'%f ');
        Trials(j).nFr(2) = length(a);
        if length(a) > nframes(j)
            a = a(1:nframes(j));
        end
        if Trials(j).nFr(2) > Trials(j).nFr(1)
            Trials(j).nrc = Trials(j).nFr(2);
        end
        id = find(eiid < sqid(j));
        if ~isempty(id)
            stvals = sscanf(txt(eiid(id(end)),6:end),'%f ');
            stvals(length(stvals+1):max(a)+1) = NaN;
            Trials(j).Stimseq = stvals(a+1);
            if eiid(id(end)) > lasteiid && Trials(j).result > 0 %first good trial of new expt
                Trials(j).Dmvals = stvals;
                lasteiid = eiid(id(end));              
            end
        else
            Trials(j).Stimseq = a;
        end
    elseif ~isempty(mtrSid)
        id = mtrSid(end);
    end
    

    id = find(bhid < sqid(j));
    if length(id)
        Trials(j).bh = sscanf(txt(bhid(id(end)),:),'bh%f');
    end
    id = find(jxid < sqid(j));
    if length(id)
        Trials(j).jx = sscanf(txt(jxid(id(end)),:),'jx%f');
    end
    id = find(soid < sqid(j));
    if length(id)
        Trials(j).so = sscanf(txt(soid(id(end)),:),'so%f %f %f %f');
    end
    id = find(orid < sqid(j));
    if length(id)
        Trials(j).or = sscanf(txt(orid(id(end)),:),'or%d');
    end
        id = find(siid < sqid(j));
    if length(id)
        Trials(j).stimi = sscanf(txt(siid(id(end)),:),'#id %d');
    end
    id = find(trid < sqid(j));
    if length(id)
        Trials(j).tr = sscanf(txt(trid(id(end)),:),'tr%f');
    end
    for k = 1:length(addfields)
        id = find(xfid{k} < sqid(j));
        if length(id)
            Trials(j).(addfields{k}) = sscanf(txt(xfid{k}(id(end)),:),[addfields{k} '%f']);
            x = Trials(j).(addfields{k});
            xidlist(j) = id(end);
        end
    end
    
    
    

    
    
    Trials(j).id = seqid(j);
    Trials(j).Start = j .* 3;
    id = find(imid < sqid(j));
    if ~isempty(id)
        a = sscanf(txt(imid(id(end)),:),'%*4s %f,%f %f %f');
        Trials(j).imseed = a(2);
        if length(a) == 4
            Stimvals.impx = a(3);
        end
    end
    id = find(sfid < sqid(j));
    a = sscanf(txt(sfid(id(end)),:),'%*2s%f');
    Trials(j).sf = a(1);
    id = find(wiid < sqid(j));
    a = sscanf(txt(wiid(id(end)),:),'%*2s%f');
    Trials(j).wi = a(1);
    id = find(seid < sqid(j));
    a = sscanf(txt(seid(id(end)),:),'%*2s%f');
    Trials(j).se = a(1);
    Trials(j).Trial = j;
end

if isfield(Trials,'nFr')
    for j = 1:length(Trials)
        nfa(j) = Trials(j).nFr(1);
        nfb(j) = Trials(j).nFr(2);
    end
    nrcframes = max([min(nfa) min(nfb)]);
    
    
    for j = 1:length(sqid)
        a = Trials(j).(rcvars{j});
        if sum(a) > 0
            if length(a) < nrcframes
                a(end+1:nrcframes) = a(end);
            else
                a = a(1:nrcframes);
            end
        end
    end
end

k = 1;
if length(exptlist)
    exid = exid(exptlist);
end
for j = 1:length(exid)
    id = find(idid > exid(j));
    tid = [];
    if length(id)
        n = find(ids(id) > 0);
        if ~isempty(n)
            tid = find([Trials.id] == ids(id(n(1))));
        end
        if isempty(tid)
            tid = find([Trials.id] > ids(id((1))));
            if length(tid)
                tid = tid(1);
                if verbose
                    fprintf('Id missing for line %d. Using %d instead of %d\n',exid(j),Trials(tid).id,ids(id((1))));
                end
            else
                fprintf('No Trials for Id > %d\n',ids(id((1))));
            end
        end
    else
        fprintf('Ids missing\n');
    end
    if length(tid) == 1
        Header.BlockStart(k) = Trials(tid).Trial;
        k = k+1;
    end
    result.exptlist(j).stim = txt(exid(j),:);
    id = find(etid < exid(j));
    if length(id)
        result.exptlist(j).et = txt(etid(id(end)),3:end);
    end
    id = find(e2id < exid(j));
    if length(id)
        result.exptlist(j).e2 = txt(e2id(id(end)),3:end);
    end
    result.exptlist(j).name = [result.exptlist(j).et 'X' result.exptlist(j).e2];
end

if ~exist('seedseq') && strcmp(exptype,'ORBW')
    return;
end
result.seedseq = seedseq;
result.seqid= seqid;
result.exptype = exptype;

    if ~isempty(Expt)
        for j = 1:length(Expt.Trials)
            id = find(seqid == Expt.Trials(j).id)
            if id
                Expt.Trials(j).Seedseq = seedseq{id(1)};
            end
        end
    else
        Expt.Trials = Trials;
        Stimvals.sf = median([Trials.sf]);
        Stimvals.wi = median([Trials.wi]);
        if isfield(Trials,'or')
            Stimvals.or = median([Trials.or]);
        end
        Expt.Stimvals = Stimvals;
        Expt.Header = Header;
    end
    Expt = FixSerialExpt(Expt, result);
if checkrpts 
    if isfield(Expt.Trials,'ob')
        id = find([Expt.Trials.ob] > 120 & abs([Expt.Trials.RespDir]) ==1);
        [result.secounts(:,1),result.secounts(:,2)] = Counts([Expt.Trials(id).se]);
        id = find(ismember(result.secounts(:,1), [1 3]));
        if length(id)
            fprintf('Uneven repeats for seeds %s\n',sprintf('%d ',result.secounts(id,2)));
        end
    end
end
return;
end

if checktrials
    fxid = strmatch('fx',txt);
    fyid = strmatch('fy',txt);
    for j = 1:length(fid)
        a = sscanf(txt(fid(j),:),'R%*c %*2s=%f %*2s=%f');
        s = regexprep(txt(fid(j),:),'.*st=none ','');
        d = sscanf(s,'%f %f');
        Trials(j).Fa = a(1);
        Trials(j).Fs = a(2);
        Trials(j).id = fid(j);
        Trials(j).Result = 1;
        Trials(j).dur = d(2);    
        Trials(j).End = d(1);
        xid = find(fyid < fid(j));
        if ~isempty(xid)
            d = sscanf(txt(fyid(xid(end)),:),'fy%f');
            Trials(j).fy = d(1);
        end
        xid = find(fxid < fid(j));
        if ~isempty(xid)
            d = sscanf(txt(fxid(xid(end)),:),'fx%f');
            Trials(j).fx = d(1);
        end
    end
    nt = j;
    for j = 1:length(bid)
        if ~isempty(strfind(txt(bid(j),:),'st=none'))
            a = sscanf(txt(bid(j),:),'R%*c %*2s=%f %*2s=%f');
            s = regexprep(txt(bid(j),:),'.*st=none ','');
            d = sscanf(s,'%f %f');
            xid = find(fid > bid(j));
            if ~isempty(xid)
                s = regexprep(txt(fid(xid(1)),:),'.*st=none ','');
                gd = sscanf(s,'%f %f');
                gd = gd(1)-d(1); %
            else
                gd = 1;
            end
            if gd > 0.9
                Trials(j+nt).Fa = a(1);
                Trials(j+nt).Fs = a(2);
                Trials(j+nt).id = bid(j);
                Trials(j+nt).Result = 0;
                Trials(j+nt).dur = d(2);
            else
                Trials(j+nt).Result = -2;
            end
            Trials(j+nt).End = d(1);
            xid = find(fyid < bid(j));
            if ~isempty(xid)
                d = sscanf(txt(fyid(xid(end)),:),'fy%f');
                Trials(j+nt).fy = d(1);
            end
            xid = find(fxid < bid(j));
            if ~isempty(xid)
                d = sscanf(txt(fxid(xid(end)),:),'fx%f');
                Trials(j+nt).fx = d(1);
            end
        else
            Trials(j+nt).Result = -1;
        end
        result.Trials = Trials;
    end
     return;
end

if checkfx
fxid = strmatch('#fx',txt);
for j = 1:length(fxid)
    line = txt(fxid(j),:);
    a = sscanf(line,'%*4s %f,%f');
    fx(j,:) = a;
end
plot(fx(:,1),'o');
hold on;
plot(fx(:,2),'ro');
result.fx = fx;
result.txt = txt;
fxid = strmatch('fy',txt); %NB has intermediates
for j = 1:length(fxid)
    line = txt(fxid(j),:);
    a = sscanf(line,'%*2s%f');
    fy(j) = a;
end
plot(fy,'g.');
return;
end
    

if checkepos
fxid = strmatch('epos',txt);
for j = 1:length(fxid)
    line = txt(fxid(j),:);
    a = sscanf(line,'%*5s %f %f %f %f %f');
    epos(j,:) = a;
end
epos(:,1:4) = epos(:,1:4)./204.8;
id = find(ismember(epos(:,5),[8 20 11]));
hold off;
plot(epos(id,1),'o');
hold on;
plot(epos(id,2),'ro');
plot(epos(id,3),'go');
plot(epos(id,4),'mo');
legend('RH','LH','RV','LV');
Expt.epos = epos;
return;
end


if checkFr
    fxid = strmatch('Fr',txt);
    for j = 1:length(fxid)
        id = find(idid < fxid(j));
        pre(j) = fxid(j)-idid(id(end));
        post(j) = idid(id(end)+1)-fxid(j);
       fid(j) = idid(id(end)+1); %id for this trial send AFTER Fr
       Fr(j) = sscanf(txt(fxid(j),3:end),'%f');
       
       stimid(j) = sscanf(txt(fid(j),3:end),'%f');
    end
    result.Fr = Fr;
    result.id = stimid;
    result.txt = txt;
    return;
end


if checkseq
    needid = strmatch('Repneed',txt);
    remid = strmatch('Remaining',txt);
    startid = strmatch('O 5',txt);
    mdid = strmatch('MD',txt);
    seid = strmatch('se',txt);
    good = abs(scores) ==1;
    id = find(pid < mdid(1));
    gs = id(end); 
    for j =1:length(mdid)
        md = sscanf(txt(mdid(j),:),'MD %d %d %f %f %d/%d %d (%d,%d,%d)');
%md(5) is n correct, md(6) is n completed        
        id = find(remid > mdid(j));
        a = find(pid < mdid(j));
        result.scores(j) = scores(a(end));
        ngood = sum(good(gs:a(end)));
        if length(id) && length(md) > 6
            sti = md(7);
            stdone(sti+1,1) = md(6);
            b = find(idid < mdid(j));
            fprintf('%d(%d):%s done(%d) = %d, N%d %s\n',sti,scores(a(end)),txt(remid(id(1)),:),sti,md(5),ngood,txt(idid(b(end)),:));
            [stimid, n, err, pos] = sscanf(txt(remid(id(1)),:),'Remaining(%d):');
            rem = sscanf(txt(remid(id(1)),pos:end),'%d');
            n = length(stdone);
            totals(j,1:n) = stdone+rem(1:n);
        else
            
        end
        %WIth Late resps, RL comes after the stimuli  are swapped, chaning
        %id. So go back before last Enstim
        id = find(startid < mdid(j));
        id = find(idid < startid(id(end)));
        if length(id)
            result.ids(j) = sscanf(txt(idid(id(end)),:),'%*2s%f');
        end
        id = find(seid < mdid(j));
        if length(id)
            result.ses(j) = sscanf(txt(seid(id(end)),:),'%*2s%f');
        end
    end
    plot(totals);
    result.totals = totals;
    return;
end


if checkTimes
    esid = strmatch('O 3',txt);
    if isempty(esid)
    esid = strmatch('O 0',txt); %Binoclean did this for a while
    end
for j = 1:length(esid)
    line = txt(esid(j),:);
    a = sscanf(line,'O %d %f %f');
    estimes(j) = a(2);
    id = find(idid < esid(j));
    ids(j) = sscanf(txt(idid(id(end)),:),'%*2s%f'); 
end
    result.estimes = estimes;
    result.ids = ids;
    return;
end
for j = 1:length(pid)
    line = txt(pid(j),:);
    a = sscanf(line,'%*2s %*3s%f%*3s%f');
    x(j) = a(1);
    id = strfind(line,'st=');
    if length(id)
    in = sscanf(line(id(1):end),'%*s %f %f %f');
    times(j) = in(1);
    durs(j) = in(2);
    end
    id = find(rwid < pid(j));
    rwsz(j) = sscanf(txt(rwid(id(end)),:),'%*2s%f'); 
    id = find(idid < pid(j));
    ids(j) = sscanf(txt(idid(id(end)),:),'%*2s%f'); 
end
find(diff(times) < 0.2);
plot(types,'o');
hold on; plot(rwsz.*4,'ro-');

gid = find(ismember(types,[1 2]));
bgid = gid(1+find(types(gid(2:end)-1) == 3)); %badfix preceding good
dt = times(bgid)-times(bgid-1);
result.bgseq = find(dt < 1);
result.types = types;
result.ids = ids;
result.rwsz = rwsz;
result.times = times;



function result = CheckExptLoops(txt, checkmode)
    qeid = find(strncmp('#qe',txt,3));
    stid = find(strncmp('Stimulus',txt,6));
    rseqid = find(strncmp('#From Run',txt,8));
    idid = find(strncmp('id',txt,2));
    rfid = find(strncmp('RF',txt,2));
    lastt = 0;
    if checkmode == 1
    for j = 1:length(stid)
        t = stid(j);
        id = find(rseqid < t & rseqid > lastt);
        lastt = t;
        if ~isempty(id)
            fprintf('*');
        end
        id = find(rfid < t);
        nt = length(id);
        id = find(idid < t);
        idl = idid(id(end));
        fprintf('%d %s %s\n',nt,txt{idl},txt{t}(1:50));
    end
    elseif checkmode == 2
    for j = 1:length(rseqid)
        t = rseqid(j);
        id = find(qeid < t);
        qe = txt{qeid(id(end))};
        ql = qeid(id(end));
        id = find(rfid < t);
        nt = length(id);
        id = find(stid < t);
        stiml = txt{stid(id(end))};
        fprintf('%d %s(%d) %s\n',nt,qe,ql,stiml);
        result.qe{j} = qe;
        result.nt(j) = nt;
    end
    end


function [rfstr, rf] = MkUfl(name, Text, varargin)
%first make .ufl file with rf boxes, so that can build pen maps the
%old way
overwrite = 0;
j = 1;
rfstr = [];
rf = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'overwrite',5)
        overwrite = 1;
    end
    j = j+1;
end
ufl = strrep(name,'.online','.ufl');
rid = strmatch('cm=rf',Text);

AddTxtFile = strrep(name,'.online','Add.txt');
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%d %s','delimiter','\n');
    fclose(fid);
    id = find(a{1} == -1);
    rfstrs = a{2}(id);
else
    rfstrs = {};
end

if isempty(rid) && isempty(rfstrs)
    return;
end
%a = textscan(Text.text(id,:),'cm=rf%f,%f:%fx%f,%fdeg pe%d %f %f%*s');
% trailing spaces seem to mess this up. text(id,1:65) works for most line
% but still barfs if a line is the wrong length
if isempty(rfstrs)
for j = 1:length(rid)
    a = sscanf(Text(rid(j),:)','cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    rfs(j,1:length(a)) = a;
end
else
    for j = 1:length(rfstrs)
        a = sscanf(rfstrs{j},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(j,1:length(a)) = a;
    end
end
% find lines suggesting RF was changed after a quantitative measure
oid = strmatch('RO',Text);
pid = strmatch('RP',Text);
sid = [oid pid strmatch('RO',Text)];
if length(sid)
    id = find(rid > max(sid))
end
for j = 1:size(rfs,2)
    rf(j) = mode(rfs(:,j));
end
if size(rfs,2) < 10
  rfstr = 'Missing RF data';
else
    rfstr = sprintf('cm=rf%.2f,%.2f:%.2fx%.2f,%.0fdeg pe%.0f %.1f,%.1f fx=%.2f,fy=%.2f\n',...
        rf(1),rf(2),rf(3),rf(4),rf(5),mode(rfs(:,6)),...
        mode(rfs(:,7)),mode(rfs(:,8)),mode(rfs(:,9)),mode(rfs(:,10)),rf);
end




if exist(ufl,'file') & ~overwrite
    return;
end


d = CreationDate(Text);
ds = [];
if d > 0
    ds = ['Created: ' datestr(d,'mm/dd/yyyy')];
end
of = fopen(ufl,'w');
if of > 0 
    fprintf(of,'%s\n',rfstr);
    for j = 1:length(rid)
%        fprintf(of,'%s\n',Text.text(id(j),:));
    end
    if ~isempty(ds)
        fprintf(of,'%s\n',ds);
    end
    fclose(of);
else
    questdlg(sprintf('Can''t Write %s',ufl),'test','OK','OK');
end


function d = CreationDate(Text)

did = strmatch('uf',Text);
if isempty(did) %online file
    did = strmatch('bt',Text);
end
if length(did)
    ds = Text(did(1),:);
    did = strfind(ds,'Creat');
    if length(did)
        d = datenum(ds(did(1)+8:end));
    else
        d = 0;
    end
else 
    d = 0;
end
