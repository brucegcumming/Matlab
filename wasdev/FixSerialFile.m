function [result, Expt] = FixSerialFile(name, varargin)
%
%Use a serial File to fix missing info for a file
%
%
maxlen = 12000;
checkfx =0;
checkFr = 0;
checkbd = 0;
checkseq = 0;
forcelen = 0;
checkTimes = 1;
readexpts = 0;
Expt = [];

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j}, 'Trials')
        Expt = varargin{j};
    elseif strncmpi(varargin{j},'chkbd',5)
        checkbd= 1;
    elseif strncmpi(varargin{j},'chkfix',4)
        checkfx = 1;
    elseif strncmpi(varargin{j},'chkseq',6)
        checkseq = 1;
    elseif strncmpi(varargin{j},'forcelen',4)
        j = j+1;
        forcelen = varargin{j};
    elseif strncmpi(varargin{j},'maxl',4)
        j = j+1;
        maxlen = varargin{j};
    elseif strncmpi(varargin{j},'readexpt',6)
        readexpts = 1;
    end
    j = j+1;
end

if ischar(name)
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
txt = a{1};
if forcelen
    for j = 1:length(txt)
        if length(txt{j}) > forcelen
        txt{j} = txt{j}(1:forcelen);
        end
    end
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


if checkbd
    bdid = strmatch('bd-1005',txt);
    hxid = strmatch('hx',txt);
    dxid = strmatch('dx',txt);
    siid = strmatch('#id',txt);
    idrpt = [];
    for j = 1:length(bdid)
        a = find(hxid < bdid(j));
        b = find(dxid < bdid(j));
        c = find(idid > bdid(j));
        d = find(siid < bdid(j));
        dx(j) = sscanf(txt{dxid(b(end))},'dx%f');
        hx(j) = sscanf(txt{hxid(a(end))},'hx%f');
        ids(j) = sscanf(txt{idid(c(1))},'id%f');
        sid(j) = sscanf(txt{siid(d(end))},'#id %f');
        if j > 1 & ids(j) == ids(j-1)
            idrpt = [idrpt j];
            hx(j) = hx(j-1);
            sid(j) = sid(j-1);
            dx(j) = dx(j-1);
        end
    end
%can get repeat Ids at end expt. Then last one in file is used, but 
%firts one is correct. Need to fix this.
    bd = dx;
    id = find(dx == 0);
    bd(sid == 131072) = -0.04;
    bd(sid == 131074) = -0.04;
    bd(sid == 131073) = 0.04;
    bd(sid == 5) = -0.04;
    bd(sid == 12) = 0.04;
    fixTrials.id = ids;
    fixTrials.bd = bd;
    
    Expt = fixTrials;
    return;
end

if checkfx
fxid = strmatch('#fx',txt);
for j = 1:length(fxid)
    line = txt{fxid(j)};
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
    line = txt{fxid(j)};
    a = sscanf(line,'%*2s%f');
    fy(j) = a;
end
plot(fy,'g.');
return;
end
    
if checkFr
    fxid = strmatch('Fr',txt);
    for j = 1:length(fxid)
        id = find(idid < fxid(j));
        pre(j) = fxid(j)-idid(id(end));
        post(j) = idid(id(end)+1)-fxid(j);
       fid(j) = idid(id(end)+1); %id for this trial send AFTER Fr
       Fr(j) = sscanf(txt{fxid(j)}(3:end),'%f');
       
       stimid(j) = sscanf(txt{fid(j)}(3:end),'%f');
    end
    result.Fr = Fr;
    result.id = stimid;
    result.txt = txt;
    return;
end


if checkseq
    needid = strmatch('Repneed',txt);
    remid = strmatch('Remaining',txt);
    mdid = strmatch('MD',txt);
    seid = strmatch('se',txt);
    for j =1:length(mdid)
        md = sscanf(txt{mdid(j)},'MD %d %d %f %f %d/%d %d (%d,%d,%d)');
        id = find(remid > mdid(j));
        if length(id) && length(md) > 6
        sti = md(7);
        stdone(sti+1,1) = md(6);
        [stimid, n, err, pos] = sscanf(txt{remid(id(1))},'Remaining(%d):');
        rem = sscanf(txt{remid(id(1))}(pos:end),'%d');
        n = length(stdone);
        totals(j,1:n) = stdone+rem(1:n);
        end
        id = find(idid < mdid(j));
        if length(id)
            result.ids(j) = sscanf(txt{idid(id(end))},'%*2s%f');
        end
        id = find(seid < mdid(j));
        if length(id)
            result.ses(j) = sscanf(txt{seid(id(end))},'%*2s%f');
        end
    end
    plot(totals);
    result.totals = totals;
    return;
end


if checkTimes
    esid = strmatch('O 3',txt);
for j = 1:length(esid)
    line = txt{esid(j)};
    a = sscanf(line,'O %d %f %f');
    estimes(j) = a(2);
    id = find(idid < esid(j));
    ids(j) = sscanf(txt{idid(id(end))},'%*2s%f'); 
end
    result.estimes = estimes;
    result.ids = ids;
    return;
end
for j = 1:length(pid)
    line = txt{pid(j)};
    a = sscanf(line,'%*2s %*3s%f%*3s%f');
    x(j) = a(1);
    id = strfind(line,'st=');
    if length(id)
    in = sscanf(line(id(1):end),'%*s %f %f %f');
    times(j) = in(1);
    durs(j) = in(2);
    end
    id = find(rwid < pid(j));
    rwsz(j) = sscanf(txt{rwid(id(end))},'%*2s%f'); 
    id = find(idid < pid(j));
    ids(j) = sscanf(txt{idid(id(end))},'%*2s%f'); 
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
