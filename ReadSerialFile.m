function [result, Expt, atxt] = ReadSerialFile(name, varargin)
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
%[result, Expt, txt] = ReadSerialFile(name, 'trials') returns a list of
%trials, not made into an expt
%
%result.totalrw sums the rewards that should have been given
%
% result lists expt types, and the values of 'id' for each trial 
% txt returns the contents of the file as a matlab character array
%result = ReadSerialFile(name,'chkFr') adds the values for Fr
%result = ReadSerialFile(name,'chkrpts') checks seed repeats for ob=130
%
%ReadSerialFile(...,'quiet') suppresses non-fatal warning messages
%
maxlen = 12000;
checkfx =0;
checkFr = 0;
checkseq = 0;
checkrpts = 0;
forcelen = 0;
checkTimes = 1;
checktrials = 0;
readexpts = 0;
needpsych = 1;
checkepos=0;
addfields = {'rw' 'mD' 'Dc'};
checkfields = {'ti' 'to'}; %include if unique >1
Expt = [];
exptlist = [];
mkufl = 0;
exptype = ''; %default
verbose = 1;
adderrflag = '-show';
datatype = 'monkey';
checkexpts = 0;
profiling = 1;
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
        if iscell(varargin{j})
            addfields = {addfields{:} varargin{j}{:}};
        else
            addfields = {addfields{:} varargin{j}};
        end
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
        adderrflag = '-silent';
    elseif strncmpi(varargin{j},'readexpt',6)
        readexpts = 1;
        if strcmp(varargin{j},'readexpts')
            readexpts = 2;
            needpsych = 0;
        end
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            exptlist = varargin{j};
            checkrpts = 1;
        end
    elseif strncmpi(varargin{j},'trials',5)
        checktrials = 1;
    end
    j = j+1;
end

result.verbose = verbose;
tt = [];

if iscell(name)
    if strncmp(name{1},'Reopened',8)
        atxt = name;
        name = 'Cell Array Input';
    else
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
    end
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
tt = TimeMark(tt, 'Start',0);

ts = now;


a = textscan(fid,'%s','delimiter','\n','bufsize',maxlen);
fclose(fid);
tt = TimeMark(tt, 'Read File', profiling);


if checkexpts
    result = CheckExptLoops(a{1});
    txt = a{1};
    return;
end

atxt = a{1};
if forcelen
    txt = txt(:,1:forcelen);
end
elseif iscellstr(name)
    txt = name;
    atxt = name;
end

id = find(strncmp('#Not in',atxt,7));
%gid = setdiff(1:size(txt,1),id);
%txt = txt(gid,:);
%get rid of lines saying Can't read - they mess up line conting
xid = find(strncmp('Can''t',atxt,5));
gid = setdiff(1:length(atxt),xid);
atxt = atxt(gid);
%txt = char(atxt);

gid = find(strncmp('RG',atxt,2));
wid = find(strncmp('RW',atxt,2));
bid = find(strncmp('RB',atxt,2));
lid = find(strncmp('RL',atxt,2));
fid = find(strncmp('RF',atxt,2));
rwid = find(strncmp('rw',atxt,2));
idid = find(strncmp('id',atxt,2));
dxid = find(strncmp('dx',atxt,2));
stid = find(strncmp('O 5',atxt,3));
psyid = find(strncmp('psyv=',atxt,5));
endid = find(strncmp('O 3',atxt,3));
exid = find(strncmp('Stimulus',atxt,8));
postid = find(strncmp('#PostPoststim',atxt,12));
prepid = find(strncmp('#Prep',atxt,5));
exendid = find(strncmp('Run end',atxt,7));
stopid = find(strncmp('#Stopped',atxt,7));

types(gid) = 1;
types(wid) = 2;
types(bid) = 3;
types = types(find(types > 0));
pid = union(gid,wid);
pid = union(pid,bid);
pid = union(pid,lid);
pid = union(pid,fid);
scores = zeros(size(pid));
scores(ismember(pid, gid)) =1; %correct choice
scores(ismember(pid, wid)) = -1; %wrong choice
scores(ismember(pid, lid)) = -2;
scores(ismember(pid, fid)) = 0; %fixation
scores(ismember(pid, bid)) = -3; %bad fixation
for j = 1:length(pid)
    if strfind(atxt{pid(j)},'Corloop')
        scores(j) = scores(j)-10;
    end
end



if mkufl
    result.rfstr = MkUfl(name,txt);
end

if readexpts
    dxseqid = [];
    stimno = [];
    for j = 1:length(idid)
        ids(j) = sscanf(atxt{idid(j)},'id%d');
    end
%only look for lines nfxxx not nf=xx
%nfxx is completed frames,written at End Stim
% COuld even used Nf fo this in more recent files?
    nfid = find(strncmp('nf=',atxt,3));
%    if ~isempty(nfid)
%        atxt(nfid) = strrep(atxt(nfid),'nf=','nf');
%    end
    xid = nfid;
    nfid = find(strncmp('nf',atxt,2));
    nfid = setdiff(nfid,xid);
    for j = length(nfid):-1:1
        nfs(j) = sscanf(atxt{nfid(j)},'nf%d');
    end
    Stimvals.e3 = 'e0';
    Stimvals.nf = prctile(nfs,95);
    Header.rc = 0;
   
    sfid = find(strncmp('sf',atxt,2));
    wiid = find(strncmp('wi',atxt,2));
    seid = find(strncmp('se',atxt,2));
    xid = find(strncmp('ser',atxt,3));
    seid = setdiff(seid,xid);
    xid = find(strncmp('seof',atxt,4));
    seid = setdiff(seid,xid);
    orid = find(strncmp('or',atxt,2));
    bhid = find(strncmp('bh',atxt,2));
    bdid = find(strncmp('bd',atxt,2));
    jxid = find(strncmp('jx',atxt,2));
    soid = find(strncmp('so',atxt,2));
    siid = find(strncmp('#id',atxt,3));
    sqid = find(strncmp('mtse',atxt,4));
    eiid = find(strncmp('mtei',atxt,4));
    mtrSid = find(strncmp('mtrS',atxt,4));
    trid = find(strncmp('tr',atxt,2));
    etid = find(strncmp('et',atxt,2));
    e2id = find(strncmp('e2',atxt,2));
    e3id = find(strncmp('e3',atxt,2));
    incid = find(strncmp('ei',atxt,2));
    inc2id = find(strncmp('i2',atxt,2));
    mdid = find(strncmp('mD',atxt,2));
    opid = find(strncmp('op=0',atxt,4));
%get a list of all expt types to track variables
    exptypes = {};
    expts = unique(atxt(etid));
    for j = 1:length(expts)
        ex = strrep(expts{j}(3:end),'*','');
        if ~ismember(ex,exptypes)
            exptypes{end+1} = ex;
        end
    end
    expts = unique(atxt(e2id));
    for j = 1:length(expts)
        ex = strrep(expts{j}(3:end),'*','');
        if ~ismember(ex,exptypes)
            exptypes{end+1} = ex;
        end
    end
    expts = unique(atxt(e3id));
    for j = 1:length(expts)
        ex = strrep(expts{j}(3:end),'*','');
        if ~ismember(ex,exptypes)
            exptypes{end+1} = ex;
        end
    end
    
    exvid = find(strncmp('exvals',atxt,2));
    for j = 1:length(addfields)
        xfid{j} = find(strncmp(addfields{j},atxt,length(addfields{j})));
    end
    
    if sum(strcmp('dx',exptypes)) && sum(strcmp('dx',exptypes))
        checkfields = {checkfields{:} 'cb' 'bd'};
    end
    for j = 1:length(checkfields)
        xid = find(strncmp(checkfields{j},atxt,length(checkfields{j})));
        if length(unique(atxt(xid))) > 1
            xfid{end+1} = xid;
            addfields{end+1} = checkfields{j};
        end
    end
    lasteiid = 0;
        otherrcs = {};
        trialmode = 0;
    if isempty(sqid)
        sqid = find(strncmp('O 3',atxt,3));
       dxseqid = find(strncmp('dx:',atxt,3));
       seseqid = find(strncmp('se:',atxt,3));
%       sqid = FindNext(dxseqid,sqid);
        if ~isempty(dxseqid)
            Header.rcvars = {'dxseq' 'ceseq'};
                 
       for j = 1:length(dxseqid)
           id = find(sqid < dxseqid(j));
           if isempty(id)
               good(j) = 0;
           else
               good(j) = id(end);
           end
       end
       sqid = sqid(good(good > 0));
            trialmode = 1;
            dxid = setdiff(dxid,sqid);
            ceid = find(strncmp('ce:',atxt,3));
            otherrcs(1).id = ceid;
            otherrcs(1).type = 'ceseq';
            otherrcs(1).pre = 4;
            id = find(strncmp('mycodes:',atxt,8));
            if ~isempty(id)
                otherrcs(2).id = id;
                otherrcs(2).type = 'mycodes';
                otherrcs(2).pre = 9;
                Header.rcvars = {Header.rcvars{:} 'mycodes'};
            end
            
        elseif ~isempty(seseqid)
            dxseqid = seseqid;
            Header.rcvars = {'seseq'};
            trialmode = 1;
            dxid = setdiff(dxid,sqid);
        end
    end

    if length(exptlist)
        remid = find(strncmp('Remaining',atxt,9));
%        exid = find(strncmp('Expt',txt);
        exendid = find(strncmp('Run end',atxt,7));
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
imid = find(strncmp('imve',atxt,4));
badimid = find(strncmp('imve 0.00',atxt,7));
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

%in human psych files, nf is sometimes not written for final trial
if length(nfs) < length(sqid)
    nfs(end+1:length(sqid)) = mean(nfs);
end

tt = TimeMark(tt, 'Found Strings', profiling);
rcvar = '';
Trials(length(sqid)).Start = 0;
if ~isempty(opid)
    x = CellToMat(strfind(atxt(opid),'+ff'));
    if sum(x ==0) < 10
        Header.flipdir = 1;
        fliptrials = [];
    elseif sum(x == 1) < 10  %all not flipped
        Header.flipdir = 0;
        fliptrials = [];
    else
        fliptrials = x;
    end
end
for j = 1:length(sqid)
    if trialmode == 1 %nf line is after O 3 .... This does not need 1 subtracting
        id = find(nfid > sqid(j));
        if isempty(id)
            id = 1;
        end
        nframes(j) = nfs(id(1));
    else
        id = find(nfid < sqid(j));
        if isempty(id)
            id = 1;
        end
        nframes(j) = nfs(id(end))-1;
    end
%O 5 lines are stid    
    if trialmode == 1 && ~isempty(psyid)%sqid lines are before O 5 for current trial
        aid = find(pid < sqid(j));
        bid = find(psyid > sqid(j));
        if isempty(aid)
            startid = sqid(j) - 10;
        else
            startid = pid(aid(end));
        end
        if isempty(bid)
            nextid  = length(atxt);
        else
            nextid = psyid(bid(1))-1;
        end
        if ~isempty(prepid)
            %#prep is in serial file just before setting stim. So anything before this
            %belongs to previous stim
            xid = find(prepid > startid & prepid < sqid(j));
            if ~isempty(xid)
                startid = prepid(xid(end));
                if length(prepid) > xid(end)+1
                    zid = prepid(xid(end)+1);
                    a = sscanf(atxt{startid}(7:end),'%f');
                    b = sscanf(atxt{zid}(7:end),'%f');
                    if length(a) > 3
                        zid = prepid(xid(end)+1);
                    end
                    if a(1) == b(1)
                        nextid = prepid(xid(end)+2);
                    end
                end
            elseif strncmp(atxt{startid},'RB',2)
                xid = find(prepid < startid);
                if ~isempty(xid)
%                    startid = prepid(xid(end));
                end                  
            end
        elseif strcmp(datatype,'Obsolete')

        xid = find(postid > startid & postid < nextid);
%this check for postpost stim is wrong for monkey files ? correct for human? 
%PostPostStim Not even in human files. So this seems obsolete
        if ~isempty(xid) %Don't go back before a line saying PostPosStim
            startid = postid(xid(1));
        end
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
    line = atxt{sqid(j)};

    if exist('ceid','var')
        if length(ceid) == length(dxseqid)
            celine = atxt{ceid(j)};
        else
            celine = [];
        end
    end
    if strncmp(line,'ce:',3)
        nc=5;
        rcvar = 'ceseq';
    elseif strncmp(line,'dx:',3)
        nc = 5;
        rcvar = 'dxseq';
    elseif strncmp(line,'se:',3)
        nc = 5;
        rcvar = 'Seedseq';
    elseif ~isempty(dxseqid)        
        a = find(dxseqid < nextid & dxseqid > startid);
        if ~isempty(a)
            dxline = atxt{dxseqid(a(end))};
            a = sscanf(dxline(5:end),'%f ');
            nc = 0;
            nrc(j) = length(a);
            Trials(j).nFr(1) = length(a);
            Trials(j).nFr(2) = nframes(j);
            if strncmp(dxline,'se:',3)
                rcvar = 'Seedseq';
            else
            rcvar = 'dxseq';
            end
            nc = -1; %done scanning
        else
            nc=0;
            nrc(j) = 0;
        end
    else
        rcvar = 'Seedseq';
        nc = 6;
    end
    if nc > 0
        a = sscanf(line(nc:end),'%f ');
        nrc(j) = length(a);
        Trials(j).nFr(1) = length(a);
        Trials(j).nFr(2) = nframes(j);
    elseif nc == 0
        Trials(j).nFr = [0 0];
    end    
    
    id = find(idid < sqid(j));
    seqid(j) = ids(id(end));
    idline(j) = id(end);
%Response string is after seedseq, by 10 lines on completed trials
% but before seqid on bad trials
%can be much later on Late resp trials
% need <= >= to handle cases where sqid is dx:
    id = find(pid <= nextid & pid >= startid);
    sid = find(stid <= nextid & stid >= startid);
    Trials(j).line(1) = sqid(j);
    if isempty(id)
        rid = find(pid > nextid);
        if ~isempty(rid) && strncmp(atxt{pid(rid(1))},'RL',2)
            fprintf('Late Resp Recorded after Prep\n');
            t = rid(1);
            d = pid(t) - sqid(j);
        else
            if ~isempty(rid) && ~isempty(sid)
                if ~isempty(prepid)
                    str = atxt{startid};
                else
                    str  = '';
                end
                fprintf('%s Next Resp is %s\n',str,atxt{pid(rid(1))});
            end
        id = find(idid < sqid(j));
        if j == length(sqid)
            xid = [];
        else
            xid = find(exendid > sqid(j) & exendid < sqid(j+1));
        end
        if isempty(sid) % no stim start either. A Prep call but never shown. Possibly end expt
            if isempty(xid)
                if j < length(sqid)
                    xid = find(stopid > sqid(j) & stopid < sqid(j+1));
                end
                if isempty(xid)
                    result = AddError(result,'Stim Prep but no response or stim start for line %d Trial %d %s\n',sqid(j),j,atxt{idid(id(end))});
                end
            else
               fprintf('Stim Prep at End Expt but unused line %d Trial %d %s\n',sqid(j),j,atxt{idid(id(end))});
            end
        else
            if isempty(xid) && needpsych %This error can happen at end expt
                result = AddError(result,'Missing Response String for line %d Trial %d %s\n',sqid(j),j,atxt{idid(id(end))});
            end
        end
        id = find(pid < sqid(j));
        if isempty(id)
            d = 100;
            t = 1;
        else
            d = 100;
            t = id(end);
        end
        end
        badstate = ~isempty(sid);
    else
        Trials(j).line(2) = stid(1);
        Trials(j).line(3) = nextid;

        if pid(id(1)) > sqid(j) && id(1) > 1
            nprelin(j) = sqid(j)-pid(id(1)-1);
            xid = find(imid < sqid(j));
            if ~isempty(xid) && imid(xid(end)) < pid(id(1)-1)
 %               id(1) = id(1)-1;
 %               nprelin(j) = -nprelin(j);
            end
        end
        if strncmp(atxt{startid},'#Prep',5)
            stimno = sscanf(atxt{startid}(7:end),'%f');
        end
        t = id(1);
        d = pid(t) - sqid(j);
        badstate = 1;
    end
    if d > 20 & pid(t) ~= -2 && length(id) > 2%only late trials should be this late
        if badstate && needpsych
            result = AddError(result,'%s Suspicious Trial (no sequence) id %d\n',name,seqid(j));
        end
%        t = id(end-1);
    end    
        
    Trials(j).stimno = stimno;
    if ~isempty(rcvar)
    rcvars{j} = rcvar;
    Trials(j).(rcvar) = a;
    Trials(j).nrc = nrc(j);
    if mean(nrc) > nfs(j)+2 %kludge to fix manual expts where nframes does not match nf
        nfs(j) = mean(nrc);
    end
    end
    seedseq{j} = a;
    for k = 1:length(otherrcs)
        oid = otherrcs(k).id;
        id = find(oid > startid & oid < nextid);
        if ~isempty(id)
        line = atxt{oid(id(end))};
        nc = otherrcs(k).pre;
        a = sscanf(line(nc:end),'%f ');
        Trials(j).(otherrcs(k).type) = a;
        else
            Trials(j).(otherrcs(k).type) = NaN;
        end
    end
    
    a = textscan(atxt{pid(t)},'%2s %[^=]=%f %[^=]=%f %[^=]=%f');
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
    id = find(exvid < sqid(j));
    if ~isempty(id)
        if length(id) > 1
        mx = sscanf(atxt{exvid(id(end-1))},'exvals %f %f %f %f');
        else
            mx = [];
        end
        x = sscanf(atxt{exvid(id(end))},'exvals%f %f %f %f');
        if length(mx) > length(x)
            x(length(mx)) = mx(end);
        end
        Trials(j).exvals = x;
    end
%    id = find(dxid < sqid(j));
%    if ~isempty(id)
%        Trials(j).dx = sscanf(atxt{dxid(id(end))},'dx%f');
%    end
    
    if scores(t) == 0 && trialmode == 1  %fixation trial for dx:  Take it
        Trials(j).result = 1;
        Trials(j).RespDir = 0;
    elseif scores(t) == -3 || scores(t) == -2 || (length(seedseq{j}) > 250 && nframes(j) < 250)
        Trials(j).result = -1;
        Trials(j).RespDir = 0;
        if scores(t) > -2
            fprintf('Trial %d ignored becuase of short RC seq %d values, but only %d frames\n',j,length(seedseq{j}),nfs(j));
        end
    elseif scores(t) == 0 %just fixation
        Trials(j).result = 1;
        Trials(j).RespDir = 0;
    elseif scores(t) < -5
        Trials(j).RespDir = 0;
        Trials(j).result = scores(t);
    elseif ~isempty(a{snfield}) %sign is in the line
        Trials(j).result = 1;
        if a{snfield}(1) == 0
            Trials(j).correct = scores(t);
            id = find(psyid < sqid(j));
            f = split(atxt{pid(t)});
            Trials(j).button = f{4}(end);
                
            if ~isempty(psyid) && ~isempty(id)               
                psyv = sscanf(atxt{psyid(id(end))},'psyv=%f');
                if psyv == 0
                    if strcmp(f{4},'sn=0R') %R button correct for - signal
                        Trials(j).RespDir = -1;
                        Trials(j).button = 'R';
                    else
                        Trials(j).RespDir = 1;
                        Trials(j).button = 'L';
                    end
                    Trials(j).rwdir = Trials(j).RespDir .* scores(t);                        
                else
                    Trials(j).RespDir = scores(t) .* sign(psyv);
                    Trials(j).rwdir = sign(psyv);
                end
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
        line = atxt{mtrSid(id(end))};
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
            stvals = sscanf(atxt{eiid(id(end))}(6:end),'%f ');
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
        Trials(j).bh = sscanf(atxt{bhid(id(end))},'bh%f');
    end
    id = find(jxid < sqid(j));
    if length(id)
        Trials(j).jx = sscanf(atxt{jxid(id(end))},'jx%f');
    end
    id = find(soid < sqid(j));
    if length(id)
        Trials(j).so = sscanf(atxt{soid(id(end))},'so%f %f %f %f');
    end
    id = find(orid < sqid(j));
    if length(id)
        Trials(j).or = sscanf(atxt{orid(id(end))},'or%d');
    end
        id = find(siid < sqid(j));
    if length(id)
        Trials(j).stimi = sscanf(atxt{siid(id(end))},'#id %d');
    end
    id = find(trid < sqid(j));
    if length(id)
        Trials(j).tr = sscanf(atxt{trid(id(end))},'tr%f');
    end
    id = find(opid < nextid);
    if length(id)
        Trials(j).OptionCode = atxt{opid(id(end))}(5:end);
    end
    for k = 1:length(addfields)
        id = find(xfid{k} < sqid(j));
        if length(id)
            if atxt{xfid{k}(id(end))}(length(addfields{k})+1) == '='
                Trials(j).(addfields{k}) = sscanf(atxt{xfid{k}(id(end))},[addfields{k} '=%f']);
            else
                Trials(j).(addfields{k}) = sscanf(atxt{xfid{k}(id(end))},[addfields{k} '%f']);
            end
            x = Trials(j).(addfields{k});
            xidlist(j) = id(end);
        end
    end
    
    
    

    if ~isempty(fliptrials)
        id = find(opid < nextid);
        if ~isempty(id)
            Trials(j).flipdir = fliptrials(id(end));
        end
    end
    
    Trials(j).id = seqid(j);
    Trials(j).tid = idline(j); %index in idid
    if ~isempty(sid)
        Trials(j).Start = sscanf(atxt{stid(sid(end))},'O 5 %f');
    else
        Trials(j).Start = 0;
    end
    id = find(imid < sqid(j));
    if ~isempty(id)
        a = sscanf(atxt{imid(id(end))},'%*4s %f,%f %f %f');
        Trials(j).imseed = a(2);
        if length(a) == 4
            Stimvals.impx = a(3);
        end
    end
    id = find(sfid < sqid(j));
    a = sscanf(atxt{sfid(id(end))},'%*2s%f');
    Trials(j).sf = a(1);
    id = find(wiid < sqid(j));
    a = sscanf(atxt{wiid(id(end))},'%*2s%f');
    Trials(j).wi = a(1);
    id = find(seid < sqid(j));
    if atxt{seid(id(end))}(3) == '='
        a = sscanf(atxt{seid(id(end))},'%*2s=%f');
    else
        a = sscanf(atxt{seid(id(end))},'%*2s%f');
    end
    Trials(j).se = a(1);
    Trials(j).Trial = j;
end

tt = TimeMark(tt, 'Made Trials', profiling);

if isfield(Trials,'nFr')
    for j = 1:length(Trials)
        nfa(j) = Trials(j).nFr(1);
        nfb(j) = Trials(j).nFr(2);
    end
    nrcframes = max([min(nfa) min(nfb)]);
    
    
    for j = 1:length(sqid)
        if isfield(Trials,rcvars{j})
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
end
tt = TimeMark(tt, 'Counted nFr', profiling);

k = 1;
if length(exptlist)
    exid = exid(exptlist);
end

lastid = 1;
%use only idid(idline), since these have matches in Trials.tid
a = bsxfun(@ge, idid(idline),exid');
%a(row,col) = idid(row) >= exid(col) 
%
[b,c] = find(diff(a) > 0);
idlist(c) = 1+b; %list of last id greater than exid, for each exid
if a(1,1) == 1 %first id is already past expt start
    if idlist(1) ==0
        idlist(1) = 1;
    else
        idlist = [1 idlist];
    end
end
idlist = idline(idlist(idlist>0));
idxlist = zeros(size(idlist));
%tidx = bsxfun(@eq,ids',[Trials.id]);
%tidx(row,coc) = idid(row) == Trials(col).id
%[b,c] = find(tidx);
%tidxlist(c) = b; %idid index that matcches trial
%id2trial(b) = c; %Trial index that mathces id
for j = (1+length(idlist)):length(exid)
    id = find(a(:,j));
    if ~isempty(id)
        idlist(j) = id(1);
    else
        idlist(j) = NaN;
    end
end
    
for j = 1:length(exid) %exvals lines
%    id = find(idid > exid(j));
    id = idlist(j);
    tid = [];
    if ~isempty(id) && ~isnan(id(1))
        n = find(ids(id) > 0);
        if ~isempty(n)
            tid = find([Trials.tid] == id(n(1)));
            idxlist(j) = ids(id(n(1)));
        end
        if isempty(tid)
            tid = find([Trials.tid] > id((1)));
            if length(tid)
                tid = tid(1);
                if exid(j) > sqid(1)
                    result = AddError(result,'Id missing for line %d. Using %d instead of %d\n',exid(j),Trials(tid).id,ids(id((1))));
                end
            else
                result = AddError(result, 'No Trials for Id > %d\n',ids(id((1))));
            end
        end
    else
        result = AddError(result,'Ids missing\n');
    end
    if length(tid) == 1
        Header.BlockStart(k) = Trials(tid).Trial;
        k = k+1;
    end
    result.exptlist(j).stim = atxt{exid(j)};
    id = find(etid < exid(j));
    if length(id)
        result.exptlist(j).et = atxt{etid(id(end))}(3:end);
    end
    id = find(e2id < exid(j));
    if length(id)
        result.exptlist(j).e2 = atxt{e2id(id(end))}(3:end);
    end
    id = find(e3id < exid(j));
    if length(id)
        result.exptlist(j).e3 = atxt{e3id(id(end))}(3:end);
    end
    id = find(incid < exid(j));
    if length(id)
        result.exptlist(j).ei = sscanf(atxt{incid(id(end))}(3:end),'%f');
    end
    id = find(inc2id < exid(j));
    if length(id)
        result.exptlist(j).i2 = sscanf(atxt{inc2id(id(end))}(3:end),'%f');
    end
    result.exptlist(j).name = [result.exptlist(j).et 'X' result.exptlist(j).e2];
end
tt = TimeMark(tt, 'Sorted exid', profiling);

if ~exist('seedseq') && strcmp(exptype,'ORBW')
    return;
end
result.seedseq = seedseq;
result.seqid= seqid;
result.exptype = exptype;
result.totalrw = sum([Trials([Trials.result] == 1).rw]);


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
    tt = TimeMark(tt, 'Checking Expt', profiling);
    if readexpts == 2
        Expt = SplitExpts(Expt, result);
    else
        Expt = FixSerialExpt(Expt, result);
    end
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
tt = TimeMark(tt, 'Finished', profiling);
return;
end

if checktrials
    openid = find(strncmp('Reopened',atxt,8));
    for j = 1:length(openid)
        sessionstart(j) = datenum(atxt{openid(j)}(14:end));
    end
    fxid = find(strncmp('fx',atxt,2));
    fyid = find(strncmp('fy',atxt,2));
    allid = sort([fid' gid' bid']);
    lastgood = 0;
    for j = 1:length(allid)
        a = sscanf(atxt{allid(j)},'R%*c %*2s=%f %*2s=%f');
        x = split(atxt{allid(j)});
        s = regexprep(atxt{allid(j)},'.*st=none ','');
        if length(x) == 10
        d(1) = sscanf(x{7},'%f');
        d(2) = sscanf(x{8},'%f');
        else            
        d(1) = sscanf(x{6},'%f');
        d(2) = sscanf(x{7},'%f');
        end
        tstart = sessionstart(1) + (d(1)-d(2)) ./(24 * 60 * 60);
        Trials(j).date = tstart;
        Trials(j).Fa = a(1);
        Trials(j).Fs = a(2);
        if atxt{allid(j)}(2) == 'B'
            if lastgood
                gd = Trials(lastgood).End -d(2);
            else
                gd = 1;
            end
            if gd > 0.9
                Trials(j).Result = 0;
            else
                Trials(j).Result = -2; %bad fixes shortly after completed trial
            end
        else
            lastgood = j;
            Trials(j).Result = 1;
        end
        Trials(j).dur = d(2);    
        Trials(j).End = d(1);
        Trials(j).Start = d(1)-d(2);
        xid = find(idid < allid(j));
        if ~isempty(xid)
            d = sscanf(atxt{idid(xid(end))},'id%f');
            Trials(j).id= d;
        end

        xid = find(fyid < allid(j));
        if ~isempty(xid)
            d = sscanf(atxt{fyid(xid(end))},'fy%f');
            Trials(j).fy = d(1);
        end
        xid = find(fxid < allid(j));
        if ~isempty(xid)
            d = sscanf(atxt{fxid(xid(end))},'fx%f');
            Trials(j).fx = d(1);
        end
        xid = find(rwid < allid(j));
        if ~isempty(xid)
            d = sscanf(atxt{rwid(xid(end))},'rw%f');
            Trials(j).rw = d(1);
        end
    end
    offset = 0;
    skips = find(diff([Trials.End]) < 0);
    for j = 1:length(skips)
        id = find(openid < allid(skips(j)+1));       
        offset(skips(j)+1:length(Trials)) = sessionstart(id(end))-sessionstart(1);
    end
    for j = 1:length(offset)
        Trials(j).date = Trials(j).date + offset(j);
    end
    if j
        nt = j;
    else
        nt = 0;
    end
    result.Trials = Trials;
    result.totalrw = sum([Trials([Trials.Result] == 1).rw]);
    result.sessions = sessionstart;
     return;
end

if checkfx
fxid = strmatch('#fx',txt);
for j = 1:length(fxid)
    line = atxt{fxid(j)};
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
    line = atxt{fxid(j)};
    a = sscanf(line,'%*2s%f');
    fy(j) = a;
end
plot(fy,'g.');
return;
end
    

if checkepos
fxid = strmatch('epos',txt);
for j = 1:length(fxid)
    line = atxt{fxid(j)};
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
       Fr(j) = sscanf(atxt{fxid(j)}(3:end),'%f');
       
       stimid(j) = sscanf(atxt{fid(j)}(3:end),'%f');
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
        md = sscanf(atxt{mdid(j)},'MD %d %d %f %f %d/%d %d (%d,%d,%d)');
%md(5) is n correct, md(6) is n completed        
        id = find(remid > mdid(j));
        a = find(pid < mdid(j));
        result.scores(j) = scores(a(end));
        ngood = sum(good(gs:a(end)));
        if length(id) && length(md) > 6
            sti = md(7);
            stdone(sti+1,1) = md(6);
            b = find(idid < mdid(j));
            fprintf('%d(%d):%s done(%d) = %d, N%d %s\n',sti,scores(a(end)),atxt{remid(id(1))},sti,md(5),ngood,atxt{idid(b(end))});
            [stimid, n, err, pos] = sscanf(atxt{remid(id(1))},'Remaining(%d):');
            rem = sscanf(atxt{remid(id(1))}(pos:end),'%d');
            n = length(stdone);
            totals(j,1:n) = stdone+rem(1:n);
        else
            
        end
        %WIth Late resps, RL comes after the stimuli  are swapped, chaning
        %id. So go back before last Enstim
        id = find(startid < mdid(j));
        id = find(idid < startid(id(end)));
        if length(id)
            result.ids(j) = sscanf(atxt{idid(id(end))},'%*2s%f');
        end
        id = find(seid < mdid(j));
        if length(id)
            result.ses(j) = sscanf(atxt{seid(id(end))},'%*2s%f');
        end
    end
    plot(totals);
    result.totals = totals;
    return;
end


if checkTimes
    esid = find(strncmp('O 3',atxt,3));
    if isempty(esid)
    esid = find(strncmp('O 0',atxt,3));
    end
for j = 1:length(esid)
    line = atxt{esid(j)};
    a = sscanf(line,'O %d %f %f');
    estimes(j) = a(2);
    id = find(idid < esid(j));
    ids(j) = sscanf(atxt{idid(id(end))},'%*2s%f'); 
end
    result.estimes = estimes;
    result.ids = ids;
    return;
end
for j = 1:length(pid)
    line = atxt{pid(j)};
    a = sscanf(line,'%*2s %*3s%f%*3s%f');
    x(j) = a(1);
    id = strfind(line,'st=');
    if length(id)
    in = sscanf(line(id(1):end),'%*s %f %f %f');
    times(j) = in(1);
    durs(j) = in(2);
    end
    id = find(rwid < pid(j));
    rwsz(j) = sscanf(atxt{rwid(id(end))},'%*2s%f'); 
    id = find(idid < pid(j));
    ids(j) = sscanf(atxt{idid(id(end))},'%*2s%f'); 
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
tt = TimeMark(tt, 'Finished', profiling);

function Expts = SplitExpts(Expt,result)

starts = Expt.Header.BlockStart;
ends = starts(2:end);
ends(end+1) = Expt.Trials(end).Trial +1;
trials = [Expt.Trials.Trial];
for j = 1:length(starts)
    id = find(trials >= starts(j) & trials < ends(j));
    Expts{j}.Trials = Expt.Trials(id);
    Expts{j}.Header = Expt.Header;
    Expts{j}.Header.Stimline = result.exptlist(j).stim;
    Expts{j}.Stimvals = Expt.Stimvals;    
    Expts{j}.Stimvals.ei = result.exptlist(j).ei;
    Expts{j}.Stimvals.i2 = result.exptlist(j).i2;
    Expts{j}.Stimvals.et = result.exptlist(j).et;
    Expts{j}.Stimvals.e2 = result.exptlist(j).e2;
    Expts{j}.Stimvals.e3 = result.exptlist(j).e3;
    f = fields(Expts{j}.Trials);
    id = find(ismember(f,{'px' 'py' 'vd' 'sf' 'tf' 'co'}));
    for k = 1:length(id) 
        vals = unique([Expts{j}.Trials.(f{id(k)})]);
        if length(vals) == 1
            Expts{j}.Stimvals.(f{id(k)}) = vals;
            Expts{j}.Trials = rmfield(Expts{j}.Trials,f{id(k)});
        end
    end
end


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


function FindNext(x,y)
%not quite working yet
%for each element of vector x, find first element in y
% that is >= x(j);
a = bsxfun(@ge, x(:),y(:)');
%a(row,col) = idid(row) >= exid(col) 
%
[b,c] = find(diff(a) > 0);
