function [result, Expt] = ReadSerialFile(name, varargin)

maxlen = 12000;
checkfx =0;
checkFr = 0;
checkseq = 0;
checkrpts = 0;
forcelen = 0;
checkTimes = 1;
readexpts = 0;
checkepos=0;
Expt = [];
exptlist = [];

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j}, 'Trials')
        Expt = varargin{j};
    elseif strncmpi(varargin{j},'chkFr',5)
        checkFr= 1;
    elseif strncmpi(varargin{j},'chkepos',4)
        checkepos = 1;
    elseif strncmpi(varargin{j},'chkfix',4)
        checkfx = 1;
    elseif strncmpi(varargin{j},'chkseq',6)
        checkseq = 1;
    elseif strncmpi(varargin{j},'chkrpt',6)
        checkrpts = 1;
    elseif strncmpi(varargin{j},'forcelen',4)
        j = j+1;
        forcelen = varargin{j};
    elseif strncmpi(varargin{j},'maxl',4)
        j = j+1;
        maxlen = varargin{j};
    elseif strncmpi(varargin{j},'readexpt',6)
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            exptlist = varargin{j};
            checkrpts = 1;
        end
        readexpts = 1;
    end
    j = j+1;
end

if iscell(name)
    Expt = {};
    for j = 1:length(name)
        fprintf('Reading %s\n',name{j});
        [result{j}, Expt{j}] = ReadSerialFile(name{j},varargin{:});
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
txt = char(a{1});
if forcelen
    txt = txt(:,1:forcelen);
end
elseif iscellstr(name)
    txt = name;
    
end


gid = strmatch('RG',txt);
wid = strmatch('RW',txt);
bid = strmatch('RB',txt);
lid = strmatch('RL',txt);
rwid = strmatch('rw',txt);
idid = strmatch('id',txt);
exid = strmatch('Stimulus',txt);
types(gid) = 1;
types(wid) = 2;
types(bid) = 3;
types = types(find(types > 0));
pid = union(gid,wid);
pid = union(pid,bid);
pid = union(pid,lid);
scores = zeros(size(pid));
scores(ismember(pid, gid)) =1;
scores(ismember(pid, wid)) = -1;
scores(ismember(pid, lid)) = -2;
for j = 1:length(pid)
    if strfind(txt(pid(j),:),'Corloop')
        scores(j) = scores(j)-10;
    end
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
    trid = strmatch('tr',txt);
    etid = strmatch('et',txt);
    e2id = strmatch('e2',txt);
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
if isempty(sqid)
    sqid = pid+1;
end
imid = strmatch('imve',txt);
%RB lines come before imse
%R[WG] com after imse.
% so first find the line that comes
for j = 1:length(sqid)
    id = find(nfid < sqid(j));
    nframes(j) = nfs(id(end))-1;
    line = txt(sqid(j),:);
    a = sscanf(line(6:end),'%f ');
    if length(a) > nframes(j)
        a = a(1:nframes(j));
    end
    while a(end) == 0
        a = a(1:end-1);
    end
    Trials(j).Seedseq = a;
    seedseq{j} = a;
    id = find(idid < sqid(j));
    seqid(j) = ids(id(end));
%Response string is after seedseq, by 10 lines on completed trials
% but before seqid on bad trials
%can be much later on Late resp trials
    id = find(pid < sqid(j)+100 & pid > sqid(j) - 10);
    if isempty(id)
        id = find(idid < sqid(j));
        fprintf('Missing Response String for line %d %s\n',sqid(j),txt(idid(id(end)),:));
        id = find(pid < sqid(j));
        d = 100;;
        t = id(end);
    else
        t = id(1);
        d = pid(t) - sqid(j);
    end
    if d > 20 & pid(t) ~= -2 && length(id) > 2%only late trials should be this late
        fprintf('%s Suspicious Trial (no sequence) id %d\n',name,seqid(j));
        t = id(end-1);
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
    
    
    
    a = textscan(txt(pid(t),:),'%2s %[^=]=%f %[^=]=%f %[^=]=%f');
    Stimvals.et = a{2}{1};
    Stimvals.e2 = a{4}{1};
    Trials(j).(a{2}{1}) = a{3}(1);
    Trials(j).(a{4}{1}) = a{5}(1);
    if scores(t) == 0 || scores(t) == -2 || length(seedseq{j}) > 250
        Trials(j).result = -1;
        Trials(j).RespDir = 0;
    elseif ~isempty(a{7})
        Trials(j).result = 1;
        Trials(j).RespDir = scores(t) .* sign(a{7}(1)-0.5);
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
k = 1;
if length(exptlist)
    exid = exid(exptlist);
end
for j = 1:length(exid)
    id = find(idid > exid(j));
    tid = [];
    if length(id)
        n = find(ids(id) > 0);
        tid = find([Trials.id] == ids(id(n(1))));
    end
    if isempty(tid)
        tid = find([Trials.id] > ids(id(n(1))));
        if length(tid)
        tid = tid(1);
        fprintf('Id missing for line %d. Using %d instead of %d\n',exid(j),Trials(tid).id,ids(id(n(1))));
        else
            fprintf('No Trials for Id > %d\n',ids(id(n(1))));
        end
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
result.seedseq = seedseq;
result.seqid= seqid;

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
