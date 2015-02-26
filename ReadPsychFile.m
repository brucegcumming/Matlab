function [Expts, Summary] = ReadPsychFile(name, varargin)
%[Expts, Summary] = ReadPsychFile(name, varargin)
%   .   ,'useallexpts');  makes expts for psych and non-psych blocks
%   .... ,'nmin',n) sets min # trials to make expt block
%File format
% R%d  0 = wrong
%      1 = correct
%      2 = late/foul
%      3 = BadFix
%      >4 special lines with stimulus/expt info
%      11 Badfix caused by a microsaccade

% R%d xx=%f yy=%f time trialdur rwsize
% xx = value for expt 1
% yy = value for expt 2
%
% R = 0 = WRONG, 1 = CORRECT, 2 = FOUL Choice 3 = BAD-FIX
% R = 100 or 101 are 0,1 but in a correction loop
% 50,51 means saccade was not required
% 53 = fixation only trial, but bad fixation
% R=109 indicates start of expt
% R9 Expt Start for Human Psych 
% R4 stimulus properties at start of expt, 
% R5 other stimulus properties
% R7 stimulus properties not in strict format - don't send these lines to
%                         textscan
% R27 Cancel Exp
% R10 = Expt End
% R8 = Expt Finished by verg


buildexpts = 1;
startday = 0;
DATA.verbose = 0;
DATA.plot.round = [0 0]; 
DATA.filename = name;
DATA.plot.mintrials = 10;
DATA.useallexpts = 0;
Expts = {};
Summary.ntrials = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noexpts',8)
        buildexpts = 0;
    elseif strncmpi(varargin{j},'useallexpts',8)
        DATA.useallexpts = 1;
    elseif strncmpi(varargin{j},'nmin',4)
        j = j+1;
        DATA.plot.mintrials = varargin{j};
    end
    j = j+1;
end

if iscellstr(name)
    for j = 1:length(name)
        fprintf('%d:reading %s\n',j,name{j});
        Expts{j} = ReadPsychFile(name{j});
    end
    return;
end
tic;
txt = scanlines(name);
gid = find(~strncmp('R7',txt,2) & ~strncmp('testflag',txt,7));
xid = setdiff(1:length(txt),gid);
if isempty(gid)
    cprintf('red','No Data in %s\n',name);
    return;
end

DATA.gid = gid;

if ~isempty(xid)
    vstr = 've=';
    ids = strfind(txt(xid),vstr);
    gotve = find(CellToMat(ids));
    
    if isempty(gotve)
        vstr = 'binoclean=';
        ids = strfind(txt(xid),vstr);
        gotve = find(CellToMat(ids));
    end
    if isempty(gotve)
        DATA.binocversion(1) = 0;
    else
    ve = sscanf(txt{xid(gotve(1))}(ids{gotve(1)}(1):end),[vstr '%f.%f']);
    DATA.binocversion(1) = ve(1);
    if length(ve) > 1
        DATA.binocversion(2)= ve(2);
    end
    end
else
    DATA.binocversion(1) = 0;
end


a = textscan(char(txt(gid))','R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
if length(a{1}) < length(gid)
    f = split(txt{gid(1)});
    nf = length(f);
    if nf == 7
        a = textscan(char(txt(gid))','R%f %[^=]=%f %2s=%f %2s=%f %f %f %f','bufsize',2048);
    elseif nf == 9
        a = textscan(char(txt(gid))','R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f');
    else
        
    endid = find(strncmp('R10 ',txt,4));
    for j = 1:length(endid)
        spaces = regexp(txt{endid(j)},'\s+');
        if length(spaces) < nf
            fprintf('Fixing Line %d\n',endid(j));
        end
        for k = length(spaces):11
            txt{endid(j)} = [txt{endid(j)} ' x=0'];
        end
        
    end
    bid = find(strncmp('R8',txt,2));
    for j = 1:length(bid)
        spaces = regexp(txt{bid(j)},'\s+');
        if length(spaces) < 11
            fprintf('Fixing Line %d\n',bid(j));
        end
        for k = length(spaces):11
            txt{bid(j)} = [txt{bid(j)} ' x=0'];
        end
        
    end
    bid = find(strncmp('R4 ',txt,3));
    for j = 1:length(bid)
        if strfind(txt{bid(j)},':.')
            fprintf('Fixing Line %d\n',bid(j));
            txt{bid(j)} = regexprep(txt{bid(j)}, '[0-9]+:[0-9,.]+:', '0');
        end
    end
    bid = find(strncmp('R5  ',txt,4));
    for j = 1:length(bid)
            fprintf('Fixing Line %d\n',bid(j));
            txt{bid(j)} = strrep(txt{bid(j)}, 'R5  ', 'R5 ve=1.1 ');
    end
    a = textscan(char(txt(gid))','R%f %[^=]=%[^ ] %[^=]=%[^ ] %[^=]=%[^ ] %f %f %f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
    a{3} = str2double(a{3});
    a{5} = str2double(a{5});
    a{7} = str2double(a{7});
    if length(a{1}) < length(gid)
        cprintf('red','Error scanning lines in %s line %d\n',name,length(a{1}));
    end
    end
end
saved = strncmp('S ',a{2},2);
a{2} = strrep(a{2},'S ','');
toc
id = regexp(name,'[0-9][0-9][A-Z]');
if length(id)
    startday = datenum(name(id(1):id(1)+8));
end
if DATA.verbose
    fprintf('\n%d lines,', length(a{1}));
end
if isnumeric(a) & a == -1 %textscan returned an error
    DATA.score = [];
    return;
end
if isempty(a{8})
    return;
end
d = dir(name);
DATA.filedate = d.datenum;
DATA.score = a{1};
DATA.saved = saved;
nxval = (length(a)-9)/2;

%
% find experiment actual date
id = find(CellToMat(strfind(txt(xid),'date=')));
progtimeoffset = 0;
ptime = 0;
bwtime = 0;
if ~isempty(id)
    s = txt{xid(id(1))};
    id = strfind(s,'date=');
    DATA.CreationDate = datenum(s(id(1)+[9:23]));
    id = strfind(s,'progtime=');
    if ~isempty(id)
        ptime = sscanf(s(id(1):end),'progtime=%f');
    end
    id = strfind(s,'bt=');
    if ~isempty(id)
        bwtime = sscanf(s(id(1):end),'bt=%f');
    end
    if ptime > 0 && bwtime > 0
        progtimeoffset = ptime-bwtime;
    end
    
else
    id = find(CellToMat(strfind(txt(xid),'time=')));
    if ~isempty(id)
        s = txt{xid(id(1))};
        id = strfind(s,'time=');
        h = sscanf(s(id(1):end),'time=%d:%d');
        DATA.CreationDate = startday + (h(1) + h(2)./60)./24;
    end
end
for j = 1:nxval
    x = a{9+j*2};
    tid = find(CellToMat(x)); %non-empty 
    id = strfind(x(tid),'(');  %fields with (..) encode extra values. removethe braces from field names
    if ~isempty(id)
        x = regexprep(x,'\(.*\)','');
    end
    y = a{9+j*2+1};
    if DATA.verbose && j == 1
        fprintf('%d trials,', length(x));
    end
    for k = 1:min([length(x) length(y)])
        if ~isempty(x{k})
            DATA.trialvals(k).(x{k})=y(k);
        end
    end
end
DATA.x = a{3};
if length(a) > 10
id = strmatch('Sa',a{11});  %%trials terminated by mirosaccade
if max(id) <= length(DATA.score)
tid = find(DATA.score(id) == 3);
DATA.score(id(tid)) = 7;
end
end
if DATA.verbose
    fprintf('%d xvals,', nxval);
end
if length(DATA.x) < length(DATA.score)
    DATA.score = DATA.score(1:length(DATA.x));
end
DATA.xtype = a{2};
if DATA.plot.round(1)
    DATA.x = round(DATA.x/DATA.plot.round(1)) .* DATA.plot.round(1);
end
DATA.ytype = a{4};
DATA.y = a{5};
if DATA.plot.round(2)
    DATA.y = round(DATA.y/DATA.plot.round(2)) .* DATA.plot.round(2);
end



id = strmatch('e0',a{4});
DATA.y(id) = 0;
DATA.sign = a{7};
DATA.times = a{8};
%DATA.times(find(DATA.score == 4)) = 0;
if length(a) > 11
    DATA.ztype = a{11};
    DATA.z = a{12};
end

if DATA.times(1) == 0
    id = find(DATA.times > 0 & DATA.score ~= 6);
    DATA.times(1) = DATA.times(id(1));
end
esid = find(a{1} ==4);
for j = 1:length(esid)
    if a{8}(esid(j)) > 1000
        DATA.times(esid(j)) = ctime2datenum(a{8}(esid(j)));
    end
end

id = find(DATA.times ==0);
while ~isempty(id)
    DATA.times(id) = DATA.times(id-1);
    id = find(DATA.times ==0);
end
DATA.times = DATA.times-progtimeoffset;
    
if DATA.verbose
    fprintf('%d jumps,', j);
end

DATA.rwszs = a{10};
DATA.DURS = a{9};
DATA.readtime = now;

tid = find(ismember(DATA.score, [4])); %start expts
oid = strmatch('tr',a{4}(tid));
trs = a{5}(tid(oid));
DATA.Blockid.tr = tid(oid);
DATA.Block.tr = trs;
DATA.Stimvals.tr = median(trs);  

  
tid = find(DATA.score == 5); %stim values
oid = strmatch('or',a{2}(tid));
ors = a{3}(tid(oid));
DATA.Blockid.or = tid(oid);
DATA.Block.or = ors;
DATA.Stimvals.or = median(ors);  
oid = strmatch('sz',a{2}(tid));
ors = a{3}(tid(oid));
DATA.Blockid.sz = tid(oid);
DATA.Block.sz = ors;
DATA.Stimvals.sz = mean(ors);  

oid = strmatch('co',a{4}(tid));
cos = a{5}(tid(oid));
DATA.Blockid.co = tid(oid);
DATA.Block.co = cos;
DATA.Stimvals.co = mean(cos);  
oid = strmatch('sf',a{4}(tid));
ors = a{5}(tid(oid));
DATA.Blockid.sf = tid(oid);
DATA.Stimvals.sf = median(ors);  
DATA.Block.sf = ors;

if length(a) > 15
oid = strmatch('dd',a{15}(tid));
if length(oid)
    xs = a{16}(tid(oid));
    DATA.Blockid.dd = tid(oid);
    DATA.Block.dd = xs;
    DATA.Stimvals.dd = mean(xs);  
end
end

if length(a) > 17
oid = strmatch('c2',a{17}(tid));
if length(oid)
    xs = a{18}(tid(oid));
    DATA.Blockid.c2 = tid(oid);
    DATA.Block.c2 = xs;
    DATA.Stimvals.c2 = mean(xs);  
end
end


DATA.Stimvals.wi = round(DATA.Stimvals.sz * 10)/10;
DATA.Stimvals.hi = round(DATA.Stimvals.sz * 10)/10;

oid = find(strncmp('ve',a{2}(tid),2));
stid = tid(oid);

if ~isempty(oid) && DATA.binocversion(1) == 0
    DATA.binocversion = mean(a{3}(stid));
end

if DATA.verbose
    fprintf('%d ves,', length(oid));
end

Summary.ntrials = sum(ismember(DATA.score,[51 1 0]));
rid = find(ismember(DATA.score,[51 1]));
Summary.goodtrials = length(rid);
Summary.TotalReward = sum(DATA.rwszs(rid));

[sid,eid] = FindBlocks(DATA);

Comments = [];
if length(oid)
        
    opid = find(CellToMat(strfind(txt(xid),'op=')));
    ves = a{3}(tid(oid));
    DATA.Blockid.ve = tid(oid);
    DATA.Block.ve = ves;
    DATA.Stimvals.ve = mean(ves);
    
    bos = a{9}(tid(oid));
    DATA.Blockid.bo = tid(oid);
    DATA.Block.bo = bos;
    DATA.Stimvals.bo = mean(bos);
    bcs = a{14}(tid(oid));
    DATA.Blockid.Bc = tid(oid);
    DATA.Block.Bc = bcs;
    DATA.Stimvals.bc = mean(bcs);
    
    bhs = a{10}(tid(oid));
    DATA.Blockid.bh = tid(oid);
    DATA.Block.bh = bhs;
    DATA.Stimvals.bh = mean(bhs);
    for j = 1:length(oid)
    end
    OptionCode = regexprep(txt(xid(opid)),'.*op=','op=');
    OptionCode = regexprep(OptionCode,' .*','');
    opid = find(CellToMat(strfind(txt(xid),'op=')));
    id = MatchInd(opid,oid);
    for j = 1:length(id)
        if id(j) > 0
            DATA.Block.OptionCode{id(j)} = OptionCode{j};
            DATA.Blockid.OptionCode(id(j)) = opid(j);
        end
    end
    ns = 0;
    for j = 1:length(sid)
        lines = sid(j):eid(j);
        id = find(xid > gid(sid(j)) & xid < gid(eid(j))); %R7 lines that apply to this expt
        ns(j) = length(id);
        if ns(j) > 400
            id
        end
        for k = 1:length(id)
            s = split(txt{xid(id(k))});
            if sum(strncmp('st=',s,3));
                l = find(strncmp('st=',s,3));
                DATA.Block.st(j) = StimulusName(s{l}(4:end));
                DATA.Blockid.st(j) = gid(sid(j));
            end
        end            
    end
end

id = find(CellToMat(strfind(txt(xid),'cm=')));
for j = id(:)'
    tid = find(gid < xid(j));
    if isempty(tid)
        tid = 1;
    else
        tid = tid(end);
    end
    p = strfind(txt{xid(j)},'cm=');
    Comments(end+1).text = txt{xid(j)}(p:end);
    Comments(end).time = DATA.times(tid);
end

if size(a,2) > 13
oid = strmatch('xo',a{11}(tid));
ors = a{12}(tid(oid));
DATA.Block.xo = ors;
DATA.Blockid.xo = tid(oid);

DATA.Stimvals.xo = median(ors);  
oid = strmatch('yo',a{13}(tid));
ors = a{14}(tid(oid));
DATA.Block.yo = ors;
DATA.Blockid.yo = tid(oid);
DATA.Stimvals.yo = median(ors);  
else
    ors = a{9}(tid(oid));
    DATA.Stimvals.xo = median(ors);  
    ors = a{10}(tid(oid));
    DATA.Stimvals.yo = median(ors);  
end

tid = find(DATA.score == 8); %background values
if length(tid)
oid = strmatch('xo',a{2}(tid));
    
end
oid = strmatch('xo',a{2}(tid));
if length(oid)
    xs = a{3}(tid(oid));
    DATA.Blockid.backxo = tid(oid);
    DATA.Block.backxo = xs;
    DATA.Stimvals.backxo = mean(xs);  
end
oid = strmatch('yo',a{4}(tid));
if length(oid)
    xs = a{5}(tid(oid));
    DATA.Blockid.backyo = tid(oid);
    DATA.Block.backyo = xs;
    DATA.Stimvals.backyo = mean(xs);  
end




tid = find(DATA.score == 6); %seed offset for images values
if length(tid)
    DATA.Stimvals.seedoffset = a{12}(tid(end)); %kludge assumes 4x2 stims. Need to handle more carefully one day
end
if DATA.verbose
    fprintf('%d rws,', length(DATA.rwszs));
end
tid = find(DATA.score == 4); %seed offset for images values
ts = a{8}(tid);
ts = 719528.8+ts./(60 * 60 * 24);  %convert sec to days
for j = 1:length(tid)
    id = find(stid < tid(j));
    if length(id)
    DATA.Blockid.Start(id(end)) = tid(j);
    DATA.Block.Start(id(end)) = ts(j);
    end
end
if isfield(DATA.Block,'Start') && max(DATA.Block.Start) > datenum('01-Jan-2000')
    for j = 2:length(DATA.Block.Start)
        if DATA.Block.Start(j) == 0
            DATA.Block.Start(j) = DATA.Block.Start(j-1);
        end
    end
elseif startday > 0
    DATA.Block.Start = startday; 
    DATA.Blockid.Start = 1; 
end

if length(DATA.rwszs) ~=  length(DATA.score)
  n = min([length(DATA.rwszs)  length(DATA.score)]);
  DATA.rwszs = DATA.rwszs(1:n);
  DATA.score = DATA.score(1:n);
end

DATA.rwsum = cumsum(DATA.rwszs .* (ismember(DATA.score,[1 51])));
if DATA.verbose
    fprintf('%d rws,', length(DATA.rwszs));
end

a = diff(find(DATA.score > 100));  %Corr loops
nconsec = 1;
maxconsec = 1;
for j = 1:length(a)
    if a(j) == 1
        nconsec = nconsec+1;
        if nconsec > maxconsec
            maxconsec = nconsec;
        end
    else
        nconsec = 1;
    end
end
if maxconsec > 20
    fprintf('Max corr loop length %d\n',maxconsec);
end

ids = find(ismember(DATA.score,[0 1 2 3 7 100 101 102 103]));
if length(ids) < length(DATA.score)/5  %not psychophysics
    DATA.ispsych = 0;
else
    DATA.ispsych = 1;
end
if DATA.binocversion == 0
    DATA.Stimvals.ve = 0;
end
if buildexpts
    DATA.Comments = Comments;
    Expts = MakeAllExpts(DATA);
else
    Expts = DATA;
end

function [sid, eid] = FindBlocks(DATA)
    
    sid = [];
    eid = [];
    if ~isfield(DATA,'score')
    return;
    end

if DATA.binocversion(1) < 0.2 || DATA.binocversion(1) > 4
    id = find(ismember(DATA.score,[4 109]));  %block starts.
    eid = find(DATA.score(1:end-1) == 8 & DATA.score(2:end) == 5); %Expt End or Cancel
else
    id = find(ismember(DATA.score,[9 109]));  %block starts.
    %Get R8 when Stim2Psych is called with flag 0. When Verg forces end
    
    eid = find(ismember(DATA.score,[8 10 27])); %Expt End or Cancel
end
if isempty(eid)
    eid = length(DATA.score);
end



if isempty(id)
    fprintf('No Expts found\n');
    return;
end

while eid(1) < id(1)
    eid = eid(2:end);
end
blkstart = id;
blkend = eid;

if isempty(id)
    id(1) = 1;
    eid(1) = length(DATA.score)-1;
elseif eid(end) < id(end)
    eid(end+1:length(id)) = length(DATA.score)-1;
end
nx = 1;
if length(eid) > length(id) && eid(1) < id(1)
    eid = eid(2:end);
end

j = 2;
while j < length(id)
    if length(eid) >= j-1
        if id(j) < eid(j-1)
            fprintf('Blocks misaligned at line %d\n',DATA.gid(id(j)));
            ng = sum(ismember(DATA.score(id(j-1):id(j)),[0 1 51]));
            if ng < 3
                fprintf('Ignoring Empty Expt Lines  %d - %d\n',DATA.gid(id(j-1)),DATA.gid(id(j)));
                id(j-1) = [];
                j = j-1;
            else
                fprintf('Cant ingore lines %d-%d - has %d trials\n',DATA.gid(id(j-1)),DATA.gid(id(j)),ng);
            end
        end
    else
        fprintf('Blocks misaligned (missing end) at %d\n',DATA.gid(id(j)));        
    end
    j = j+1;
end

sid = id;
if length(eid) < length(sid) %crashes left expts empty
    fprintf('Last end Block at %d, last start at %d\n',DATA.gid(eid(end)),DATA.gid(sid(end)));
    eid(1:length(sid)-1) = sid(2:end);
    eid(length(sid)) = length(DATA.score);
end

function Expts = MakeAllExpts(DATA, varargin)

Expts = {};
psychonly = 0;
if ~isfield(DATA,'score')
    return;
end
[id, eid] = FindBlocks(DATA);
if psychonly
    sid = find(DATA.score <= 5);
else
    sid = find(DATA.score <= 5 | ismember(DATA.score,[51 50]));
end

%Check for unterminated blocks e.g. becuase of crash
for j = 1:length(id)    
    if j < length(id) && eid(j) > id(j+1)
        eid(j+1:end+1) = eid(j:end);
        eid(j) = id(j+1)-1;
    end
end
nx = 1;
for j = 1:length(id);
    Expt = [];
    tmp = DATA;
%this union seems to conflate all the blocks.  ? not sid intented   ?
%intersect? Changed to intersect Jan 2015. But tests in MakeExpt anyway
    ids = intersect([id(j):eid(j)],sid);
    ids = id(j):eid(j);
    tmp.score = DATA.score(ids);
    if isfield(DATA,'trialvals')
        tmp.trialvals = DATA.trialvals(ids);
    end
    tmp.x = DATA.x(ids);
    tmp.y = DATA.y(ids);
    if isfield(DATA,'ztype')
        tmp.z = DATA.z(ids);
        tmp.ztype = DATA.ztype(ids);
    end
    tmp.xtype = DATA.xtype(ids);
    tmp.ytype = DATA.ytype(ids);
    tmp.sign = DATA.sign(ids);
    tmp.times = DATA.times(ids);
    tmp.rwszs = DATA.rwszs(ids);
    tmp.rwsum = DATA.rwsum(ids);
    tmp.DURS = DATA.DURS(ids);
    tmp.saved = DATA.saved(ids);
    Expt = MakeExpt(tmp, varargin{:});
    f = fields(DATA.Block);
    n = (j * 2)-1;
    if isfield(Expt,'Trials') && length(Expt.Trials) > DATA.plot.mintrials
        goodexpt(j) = 1;
    else 
        goodexpt(j) = 0;
    end
    if DATA.score(eid(j)) == 27 && ~DATA.useallexpts
        goodexpt(j) = 0;
    end
    if goodexpt(j)
        for k = 1:length(f)
            bid = find(DATA.Blockid.(f{k}) <= id(j));
            if length(bid) & bid(end) <= length(DATA.Block.(f{k}))
                Expt.Stimvals.(f{k}) = DATA.Block.(f{k})(bid(end));
            else
                Expt.Stimvals.(f{k}) = NaN;
            end
        end
        Expts{nx} = Expt;
        if DATA.score(eid(j)) == 27
           Expts{nx}.Header.result = 19;
        else
           Expts{nx}.Header.result = 2;
        end
        nx = nx+1;
    end
end
if ~isempty(DATA.Comments)
    ts = 0;
    for j = 1:length(Expts)
        te = Expts{j}.Trials(end).End(end)+10000;
        cid = find([DATA.Comments.time] > ts & [DATA.Comments.time] <= te);
        Expts{j}.Comments.text = {DATA.Comments(cid).text};
        Expts{j}.Comments.times = [DATA.Comments(cid).time];
        ts = te;
    end
end

if DATA.verbose
    fprintf('%d Expts Made\n',nx);
end

function Expt = MakeExpt(DATA, varargin)

useall = 0;
skipblock = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'skipblock',5)
        j = j+1;
        skipblock = varargin{j};
    elseif strncmpi(varargin{j},'useall',5)
        useall = 1;
    end
    j = j+1;
end
if DATA.useallexpts
    sid = ismember(DATA.score, [1 2 51]);
    useall = 1;
else
    sid = ismember(DATA.score, [1 2]);
end
Expt = [];

if sum(sid) == 0 %no psych trials
    if useall
        sid = find(DATA.score ==51);
        if isempty(sid)
            return;
        end
    else
        return;
    end
end
[et, p] = GetType(DATA.xtype(sid));
[e2, p] = GetType(DATA.ytype(sid));
if isfield(DATA,'ztype')
[e3, p] = GetType(DATA.ztype(sid));
else
    e3 = 'e0';
end

Expt.Stimvals.et = et;
if ~isfield(Expt.Stimvals,'e2')
Expt.Stimvals.e2 = 'e0';
end

if ~isfield(Expt.Stimvals,'e3')
    if sum(strncmp(e3,{'se'},2))
        Expt.Stimvals.e3 = 'e0';
    else
        Expt.Stimvals.e3 = e3;
    end
end
chkf = {'ei' 'i2' 'i3'};
defaults = {0 0 0};
for j = 1:length(chkf)
if ~isfield(Expt.Stimvals,chkf{j})
    Expt.Stimvals.(chkf{j}) = defaults{j};
end
end

fn = fields(DATA.Stimvals);
for j = 1:length(fn)
    Expt.Stimvals.(fn{j}) = DATA.Stimvals.(fn{j});
end

Expt.Header.rc = 0;
Expt.Header.expname = DATA.filename;
if useall
    usetrials = [0 1 3 51];
else
    usetrials = [0 1 3];
end
id = find(ismember(DATA.score,usetrials));
if length(id) > 1 && (length(unique(DATA.y(id))) > 1 || DATA.Stimvals.ve > 5 || DATA.Stimvals.ve < 2)
    Expt.Stimvals.e2 = e2;
    blocks(1) = id(1);
else
    Expt = [];
    return;
end
nb = 1;
if isfield(DATA,'trialvals')
f = fields(DATA.trialvals);
else 
    f = [];
end

if DATA.verbose
    fprintf('Filling %d trials..',length(id));
end
if id(end) > length(DATA.y)
    fprintf('Y too short!!\n');
end
if id(end) > length(DATA.x)
    fprintf('X too short!!\n');
end
if id(end) > length(DATA.sign)
    fprintf('Sign too short!!\n');
end
Expt.Header.ntrials = length(id);
  for j = 1:length(id)
      if id(j) > 1 && DATA.score(id(j)-1) == 5
          blocks(nb) = id(j);
          nb = nb+1;
      end
      Expt.Trials(j).(et) = DATA.x(id(j));
      if ~strcmp(Expt.Stimvals.e2,'e0')
          Expt.Trials(j).(e2) = DATA.y(id(j));
      end
          
      Expt.Trials(j).rwdir = DATA.sign(id(j));
      if DATA.score(id(j)) == 1
          Expt.Trials(j).RespDir = DATA.sign(id(j));
      elseif DATA.score(id(j)) == 0
          Expt.Trials(j).RespDir = -DATA.sign(id(j));
      else
          Expt.Trials(j).RespDir = 0;
          Expt.Trials(j).rwdir = 0;
      end
      Expt.Trials(j).Start = DATA.times(id(j)).*10000;
      Expt.Trials(j).End = (DATA.times(id(j)).*10000)+20000;
      Expt.Trials(j).Trial = id(j);
      Expt.Trials(j).rw = DATA.rwszs(id(j));
      Expt.Trials(j).rwsum = DATA.rwsum(id(j));
      for k = 1:length(f)
          if length(DATA.trialvals(id(j)).(f{k}))
              Expt.Trials(j).(f{k}) = DATA.trialvals(id(j)).(f{k});
          end
      end
      if isfield(Expt.Trials,'TwoCylDisp') && isfield(Expt.Trials,'hx')
          Expt.Trials(j).dx = Expt.Trials(j).TwoCylDisp;
          Expt.Trials(j).bd = Expt.Trials(j).TwoCylDisp;
          if Expt.Trials(j).TwoCylDisp == -1011  %flip
              Expt.Trials(j).dx = abs(Expt.Trials(j).hx) .* DATA.sign(id(j));
              Expt.Trials(j).bd = -Expt.Trials(j).dx;
          elseif Expt.Trials(j).TwoCylDisp == 0  %flip
              Expt.Trials(j).bd = Expt.Trials(j).hx;
          end
          Expt.Trials(j).rd = Expt.Trials(j).dx - Expt.Trials(j).bd;
      end
      Expt.Trials(j).OptionCode = '+2a';
      Expt.Trials(j).Spikes = 0;
  end
  blocks(nb) = max([Expt.Trials.Trial]);
if DATA.verbose
    fprintf('..Done\n');
end
  if skipblock
      t = blocks(skipblock+1)-1;
      Expt.Header.BlockStart = blocks(skipblock+1:end);
      Expt.Trials = Expt.Trials(t:end);
  else
      Expt.Header.BlockStart = blocks;
  end
 id = find([Expt.Trials.rwdir] == 0);
 Expt.Header.Start = Expt.Trials(1).Start(1);
 Expt.Header.End = Expt.Trials(end).End(end);
 Expt.Header = CopyFields(Expt.Header,DATA,{'CreationDate'});
 %When reading an onlien file, don't know if smr files were split or not
 % treat them as if not, sp keep one creationdate, then references all
 % times to that.
 if 0
 id = find(DATA.score == 4 & DATA.times > 1000);
 if ~isempty(id)
     Expt.Header.CreationDate = DATA.times(id(1));
 end
 end
 
 function [type,frac] = GetType(x)

v = unique(x);
for j = 1:length(v)
    n(j) = length(strmatch(v{j},x));
end
[a,b] = max(n);
type = v{b};
frac = a./sum(n);

function dn = ctime2datenum(x)
         
    
    dn = x./(24 * 60 * 60);
%dont understand where 20/24 comes from. Someone doesn't use midnight    
    dn = dn + 719528 + (20/24); %add 1 Jan 1970
    