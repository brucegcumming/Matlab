function out = PsychMon(varargin)
% X = PsychMon(filename)
% Read and plots psych files made by binoc for remote monitoring
%
% X = PsychMon(filename,'getexpts')
% returns a cell array of Expt structures, one for each block.
% X = PsychMon(filename,'getexpts','useall')
%              includes trials that are fixation only
%
%PsychMon(filename,'name', name) 
%  Adds name to figure tags, so that you can monitor more than one
%  animal/file
%
% binoc writes a file summarizing psych in a textscan friendly format
%
% R%d xx=%f yy=%f time trialdur rwsize
% xx = value for expt 1
% yy = value for expt 2
%
% R = 0 = WRONG, 1 = CORRECT, 2 = FOUL Choice 3 = BAD-FIX
% R = 100 or 101 are 0,1 but in a correction loop
% 50,51 means saccade was not required
% R=4 means start of expt, xx = Covary Xpos, yy = target ratio 
name = 'PsychMon';
fname = 'PsychMon';
%
%
%If strings is empty, no list is shown. Otherwise a listbox is included.
%
strings = { 'ruf915', 'ruf922'};
strings = {};

tag = [name 'Panel']; %need to change these.
init = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'name',4)
        j = j+1;
        name = varargin{j};
        tag = [name 'Panel']; %need to change these.
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        tag = varargin{j}; %need to change these.
    end
    j = j+1;
end

if ishandle(varargin{1})
    toplevel = varargin{1};
    tmp = get(toplevel,'UserData');
    while ~isfigure(toplevel) & ~isfigure(tmp)
        toplevel = get(toplevel,'Parent');
        tmp = get(toplevel,'UserData');    
    end
    if isfigure(tmp)
        toplevel = tmp;
        DATA = get(toplevel,'UserData');
    end
    
    varargin = varargin(3:end);
    init = 0;
else
    toplevel = findobj('Tag',tag); %%GUI already up
end


if length(varargin) > 1 && strncmpi(varargin{2},'getexpt',7)
        DATA.filename = varargin{1};
elseif ~isempty(toplevel)
  if strncmpi(varargin{1},'store',5)
    set(toplevel,'UserData',varargin{2});
    DATA = varargin{2};
  else
    DATA = get(toplevel,'UserData');
  end
else
    if iscell(varargin{1})
        DATA.fstrings = varargin{1};
        DATA.filename = DATA.fstrings{1};
    else
    DATA.filename = varargin{1};
    end
    DATA.figtag = [name 'Online'];
    DATA.plot.noplot = 0;
    DATA.online = 1;
end


if nargin
    if ishandle(varargin{1})
        varargin{1} = varargin{3}; %temp kludge for new callback syntax
    end
    if strncmpi(varargin{1},'update',5)
        update(DATA);
    elseif strncmpi(varargin{1},'close',5)
        CloseTag(DATA.figtag);
        CloseTag(tag);
        stop(DATA.timerobj);
    elseif strncmpi(varargin{1},'getexp',5) %% Get expt from DATA currently associated with figure
        out = MakeExpt(DATA,varargin{:});
        return;
    elseif strncmpi(varargin{1},'newblock',5) %% tell binoc to run another block
        [a,b] = fileparts(DATA.filename);
        cmdfile = [a '/exptcmd'];
        fid = fopen(cmdfile,'w');
        fprintf(fid,'!expt\n');
        fclose(fid);
        
    elseif strncmpi(varargin{1},'reloadifnew',8)
        if DATA.plot.autoplot
        d = dir(DATA.filename);
        if d.datenum > DATA.readtime 
            DATA = ReadFile(DATA, DATA.filename);
            set(DATA.toplevel,'UserData',DATA);
        end
        end
    elseif strncmpi(varargin{1},'reload',5)
        ReadFile(DATA, DATA.filename);
    elseif strncmpi(varargin{1},'prev',4)
        DATA.current = DATA.current-1;
        set(DATA.gui.filename,'string',DATA.fstrings{DATA.current});
        DATA.filename = DATA.fstrings{DATA.current};
        ReadFile(DATA, DATA.fstrings{DATA.current});
    elseif strncmpi(varargin{1},'next',4)
        DATA.current = DATA.current+1;
        set(DATA.gui.filename,'string',DATA.fstrings{DATA.current});
        DATA.filename = DATA.fstrings{DATA.current};
        ReadFile(DATA, DATA.fstrings{DATA.current});
    elseif strncmpi(varargin{1},'setplot',5)
        ReadFile(DATA, DATA.filename);
    else
        j = 1;
        init = 1;
        while(j <= length(varargin))
            if(strncmpi(varargin{j},'name',3))
                j = j+1;
                name = varargin{j};
    elseif strncmpi(varargin{j},'getexpts',8)
          DATA.name = name;
          DATA.online = 0;

        DATA = SetDefaults(DATA);
        DATA.plot.noplot = 1;
        DATA = ReadFile(DATA,DATA.filename);
        out = MakeAllExpts(DATA,varargin{j+1:length(varargin)});
        return;
    elseif strncmpi(varargin{j},'getexpt',4)
        DATA.online = 0;
          DATA.name = name;

        DATA = SetDefaults(DATA);
        DATA.plot.noplot = 1;
          DATA = ReadFile(DATA,DATA.filename);
          if ~isempty(DATA.score)
          out = MakeExpt(DATA,varargin{j:end});
          else
              out = [];
          end
        return;
    elseif strncmpi(varargin{j},'mintrials',5)
        j = j + 1;
        DATA.plot.mintrials = varargin{j};
    elseif strncmpi(varargin{j},'nmin',4)
        j = j+1;
        DATA.plot.nmin = varargin{j};
        elseif strncmpi(varargin{j},'noplot',5)
            DATA.plot.noplot = 1;
            end
            j = j+1;
        end
    end
else
    init = 1;
end

if init & isempty(findobj('Tag',tag))
    DATA.tag = tag;
    DATA.name = name;
    DATA = SetDefaults(DATA);
  if DATA.plot.noplot == 0
      DATA.fcn =  'PsychMon';
  DATA = InitInterface(DATA);
     DATA.timerobj = timer('timerfcn',{@updateplot, DATA.tag},'period',2,'executionmode','fixedspacing');
     start(DATA.timerobj);
  end
  DATA = ReadFile(DATA,DATA.filename);
  out = DATA;
end


function DATA = SetDefaults(DATA)
  DATA.plot.smooth = 10;
  DATA.plot.round = [0 0];
  DATA.plot.autoplot = 0;
  DATA.plot.flip = 0;
  DATA.plot.byrw = 0;
  DATA.plot.timerange = [0 0];
  DATA.plot.shown = 1;
  DATA.plot.rwtype = 1;  
  DATA.plot.alttype = 1;
  DATA.current = 1;
  DATA.verbose = 0;
  if ~isfield(DATA.plot,'nmin')
  DATA.plot.nmin = 5;
  end
  if ~isfield(DATA.plot,'mintrials')
  DATA.plot.mintrials = 10;
  end
  DATA.readtime = 0;

function updateplot(tag,varargin)

 
PsychMon('reloadifnew','tag', varargin{2});

function [type,frac] = GetType(x)

v = unique(x);
for j = 1:length(v)
    n(j) = length(strmatch(v{j},x));
end
[a,b] = max(n);
type = v{b};
frac = a./sum(n);

function Expts = MakeAllExpts(DATA, varargin)

Expts = {};
if ~isfield(DATA,'score')
    return;
end
id = find(ismember(DATA.score(2:end),[0 1 2 3 7 100 101 102 103]) & ismember(DATA.score(1:end-1),[4 5])); %block starts.
eid = find(ismember(DATA.score(1:end-1),[0 1 2 3 7 100 101 102 103]) & ismember(DATA.score(2:end),[4 5])); %block end.
sid = find(DATA.score == 5);
if isempty(id)
    id(1) = 1;
    eid(1) = length(DATA.score)-1;
else
eid(length(id)) = length(DATA.score)-1;
end
nx = 1;
if length(eid) > length(id) && eid(1) < id(1)
    eid = eid(2:end);
end
for j = 1:length(id);
    Expt = [];
tmp = DATA;
ids = union([id(j):eid(j)+1],sid);
tmp.score = DATA.score(ids);
if isfield(DATA,'trialvals')
tmp.trialvals = DATA.trialvals(ids);
end
tmp.x = DATA.x(ids);
tmp.y = DATA.y(ids);
tmp.xtype = DATA.xtype(ids);
tmp.ytype = DATA.ytype(ids);
tmp.sign = DATA.sign(ids);
tmp.times = DATA.times(ids);
tmp.rwszs = DATA.rwszs(ids);
tmp.rwsum = DATA.rwsum(ids);
tmp.DURS = DATA.DURS(ids);
Expt = MakeExpt(tmp, varargin{:});
f = fields(DATA.Block);
n = (j * 2)-1;
if isfield(Expt,'Trials') && length(Expt.Trials) > DATA.plot.mintrials
    for k = 1:length(f)
        bid = find(DATA.Blockid.(f{k}) <= id(j));
        if length(bid) & bid(end) <= length(DATA.Block.(f{k}))
        Expt.Stimvals.(f{k}) = DATA.Block.(f{k})(bid(end));
        else
            Expt.Stimvals.(f{k}) = NaN;
        end
    end
    Expts{nx} = Expt;
    nx = nx+1;
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
sid = ismember(DATA.score, [1 2]); 
if sum(sid) == 0 %no psych trials
    if useall
        sid = find(DATA.score ==51);
    else
    Expt = [];
    return
    end
end
[et, p] = GetType(DATA.xtype(sid));
[e2, p] = GetType(DATA.ytype(sid));

Expt.Stimvals.et = et;
Expt.Stimvals.e2 = 'e0';
Expt.Stimvals.e3 = 'e0';
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
  
function DATA = PlotData(DATA)



if DATA.verbose
    fprintf('Plotting DATA...');
end
if ~DATA.plot.noplot
h = GetFigure(DATA.figtag);
figure(h);
hold off;
it = findobj('Tag','plottype','Parent',DATA.toplevel);
plottype = get(it,'value');
else
    plottype = 0;
end

if DATA.plot.timerange(1) > 0
    tid = find(DATA.times >= DATA.plot.timerange(1));
else
    tid = 1:length(DATA.times);
end
if DATA.plot.timerange(2) > 0
    rid = find(DATA.times <= DATA.plot.timerange(2));
else
    rid = 1:length(DATA.times);
end
tid = intersect(tid,rid);
if isempty(tid)
    return;
end
if DATA.verbose
    fprintf('Maxtid %d, score %d,',max(tid),length(DATA.score));
end
tid = tid(find(tid < length(DATA.score)));
PLOT.score = DATA.score(tid);
PLOT.y = DATA.y(tid);
PLOT.x = DATA.x(tid);
PLOT.times = DATA.times(tid);
PLOT.sign = DATA.sign(tid);
PLOT.rwszs = DATA.rwszs(tid);
PLOT.dur = DATA.DURS(tid);
sids = find(DATA.score == 4); %% Start of Expt
eids = find(DATA.score == 5); %% End of Expt
stime = DATA.DURS(sids(end));
stime = floor(stime) + mod(stime,1)/60; % Expt start time in hours
stime = stime - DATA.times(sids(end))/(60 * 60); % time = 0 in hours
stime = floor(now) + stime/24;
% New way to calculate start time. Should give better estimate of time of
% the last trial.....
stime = DATA.filedate - DATA.times(end) .* 1/(24 * 60 * 60);
PLOT.times = stime + (PLOT.times  .* 1/(24 * 60 * 60));
nexpstim = DATA.rwszs(sids(end));
ndone = sum(DATA.score(sids(end):end) < 2);
ids = find(ismember(PLOT.score,[0 1 2 3 7 100 101 102 103]));
rwsum = sum(PLOT.rwszs(find(ismember(PLOT.score,[1 51]))));
if length(ids) < length(PLOT.score)/5  %not psychophysics
ids = find(ismember(PLOT.score,[51 53]));
ndone = sum(DATA.score(sids(end):end) == 51);
b = histc(PLOT.score(ids),[51 53]);
perfstr = sprintf('%s %d/%d correct  Expt %d/%d %.1fml',datestr(DATA.readtime,'HH:MM'),b(1),b(2)+b(1),ndone,nexpstim,rwsum)
PLOT.ispsych = 0;
else
b = hist(PLOT.score(ids),[0 1 2 3 7 101 102 103]);
perfstr = sprintf('%s %d/%d correct + %d Bad/Foul/Late %d microsacc Expt %d/%d %.1fml',datestr(DATA.readtime,'HH:MM'),b(2),b(2)+b(1),b(3)+b(4)+b(5)+b(6)+b(7)+b(8),b(5),ndone,nexpstim,rwsum);
PLOT.ispsych = 1;
end
shortperfstr = sprintf(' %d/%d %.1fml',ndone,nexpstim,rwsum);
    PLOT.exptover = 0;
if eids & eids(end) > sids(end) & ~DATA.plot.noplot
    perfstr = [perfstr ' Over']
    shortperfstr = [shortperfstr ' Over']
    PLOT.exptover = 1;
elseif PLOT.score(end) > 99  %in correction loop
    perfstr = [perfstr '*']
    shortperfstr = [shortperfstr '*']
end

if ismember(plottype,[0 1])
    ids = find(PLOT.score < 3  & PLOT.sign ~= 0);
    if isempty(ids) %no psych. Use fixation trials
        ids = find(PLOT.score == 51);
    end
    rws = unique(PLOT.rwszs(ids));
    if DATA.plot.rwtype > 1 & DATA.plot.rwtype <= length(rws)+1
        ids = find(PLOT.score < 3 & PLOT.rwszs == rws(DATA.plot.rwtype-1));
        perfstr = [perfstr sprintf('  rw %.2f',rws(DATA.plot.rwtype-1))];
    end
    b = hist(PLOT.score(ids),[0 1 2]);
    if DATA.plot.flip
        e2v = unique(PLOT.x(ids));
    else
        e2v = unique(PLOT.y(ids));
    end
    np = 0;
    if plottype > 0
    hold off;
    end
    colors = mycolors;
            nl = 1;
            pp = [];
            labels = {};
    for j = 1:length(e2v)
        yid = find(PLOT.y(ids) == e2v(j) & PLOT.score(ids) < 2);
        if(yid)
            xv = unique(PLOT.x(ids(yid)));
            nnp = np+1;
            for k = 1:length(xv)
                xid = find(PLOT.x(ids(yid)) == xv(k));
                xid = ids(yid(xid));
                if length(xid) >= DATA.plot.nmin
                    np = np+1;
                    pp(np).x = xv(k);
                    pp(np).y = e2v(j);
                    pp(np).n = length(xid);
                 pp(np).ncorr  = sum(PLOT.score(xid));
% if score -sign > 0, means choice was up
                 pp(np).resp  = sum(PLOT.score(xid) == (PLOT.sign(xid) > 0));
                 pp(np).p = pp(np).resp/pp(np).n;
                end
            end
            if length(xv) > 1 & plottype == 1 & np >= nnp;
            plot([pp(nnp:np).x],[pp(nnp:np).p],'o-','color',colors{j});
            if DATA.plot.shown
                for k = [nnp:np]
                   h =  text(pp(k).x,pp(k).p+0.05,sprintf('%d/%d',pp(k).resp,pp(k).n));
                   set(h,'color',colors{j});
                end
            end
            labels{nl} = sprintf('%.3f',e2v(j));
            nl = nl+1;
            hold on;
            end
        end
    end
    DATA.pp = pp;
    if plottype > 0
    set(gca,'ylim',[0 1]);
    title(perfstr);
    legend(labels);
    end
elseif plottype == 3  %reward sizes (actual)
    ids = find(PLOT.score < 2);
    hold off;
    plot(PLOT.times(ids),PLOT.rwszs(ids),'ro');
    hold on;
    plot(PLOT.times(ids),PLOT.rwszs(ids) .* PLOT.score(ids),'o');
    title(perfstr);

elseif plottype == 5  %result types
    ids = find(PLOT.score < 8);
    fprintf('%d Trails\n',length(ids));
    plot(PLOT.times(ids),DATA.score(ids),'o');
    ms = PLOT.score(ids) == 7;
    ns(1) = sum(ms);
    ms = smooth(ms.*10,DATA.plot.smooth);
    hold on;
    h(1) = plot(PLOT.times(ids),ms,'r');

    ms = PLOT.score(ids) == 0;
    ns(2) = sum(ms);
    ms = smooth(ms.*10,DATA.plot.smooth);
    hold on;
    h(2) = plot(PLOT.times(ids),ms,'g');
    
    ms = PLOT.score(ids) == 3;
    ns(3) = sum(ms);

    ms = smooth(ms.*10,DATA.plot.smooth);
    hold on;
    h(3) = plot(PLOT.times(ids),ms,'m');
    
    ms = PLOT.score(ids) == 1;
    ns(4) = sum(ms);
    ms = smooth(ms.*10,DATA.plot.smooth);
    hold on;
    h(4) = plot(PLOT.times(ids),ms,'b');

    legend(h,{sprintf('Saccade %d',ns(1)) sprintf('Wrong %d',ns(2)) sprintf('BadFix %d',ns(3)) sprintf('Correct %d',ns(4))});
    b = hist(DATA.score(ids),[0 1 2 3 7]);
    title(perfstr);

elseif plottype == 6  %cumulative reward
    ids = find(PLOT.score < 2 | PLOT.score == 51 | PLOT.score ==53); %rewared
    hold off;
    score = PLOT.score(ids);
    score(score ==51) = 1;
    score(score ==53) = 0;
    fprintf('%d data',length(ids));
    rwt = cumsum(PLOT.rwszs(ids) .* score);
    
    plot([PLOT.times(ids)' PLOT.times(end)],[rwt' rwt(end)]); %make sure last time is plotted, even if not rewarded
    hold on;
    cid = find(PLOT.score > 99); %correction loops
    rws = interp1(PLOT.times(ids), rwt, PLOT.times(cid));
    plot(PLOT.times(cid),rws,'r.');
    fprintf('Max %d vs %d',ids(end), length(PLOT.times));

    if DATA.score(end) == 5
        hold on;
        plot([PLOT.times(end) PLOT.times(end)],[0 rwt(end)],'r');
    end
    title(perfstr);
    datetick('x','HH:MM');
    fprintf('Done at %s\n',datestr(PLOT.times(end)));
    
elseif plottype == 7  % Trial signal/score
    if PLOT.ispsych == 0
        gid = find(PLOT.score ==51);
        bid = find(PLOT.score ==53);
        plot(PLOT.times(gid),PLOT.dur(gid),'o');
        hold on;
        plot(PLOT.times(bid),PLOT.dur(bid),'ro');        
        datetick('x','HH:MM');
        exid = find(PLOT.score == 4);
        for j = 1:length(exid)
            line([PLOT.times(exid(j)) PLOT.times(exid(j))],get(gca,'ylim'),'linestyle',':');
        end
    else
    respdir = zeros(size(PLOT.score));
    scores = 2.*(PLOT.score -0.5);
    subplot(1,1,1);
    if DATA.plot.alttype == 1
    ti = PLOT.times;
    else
    ti = 1:length(DATA.times);
    end
    if DATA.plot.timerange(1) > 0 && diff(DATA.plot.timerange) > 0
        tid = (ti >= DATA.plot.timerange(1) & ti <= DATA.plot.timerange(2))';
    else 
        tid = ones(size(DATA.times));
    end
        tid = ones(size(PLOT.times));
    %sign = 0 for fixation only trials
    ids = find(PLOT.sign == -1 & PLOT.score == 1);
    respdir(ids) = -1;
    ids = find(PLOT.sign == 1 & PLOT.score == 0);
    respdir(ids) = -1;
    ids = find(PLOT.sign == 1 & PLOT.score == 1);
    respdir(ids) = 1;
    ids = find(PLOT.sign == -1 & PLOT.score == 0);
    respdir(ids) = 1;
    
    [a,b] = Counts(DATA.ytype)
    [c,d] = max(a);
    if strcmp(b{d},'ob')
        signal = sd2cv(PLOT.y);
        signal(signal < 0.001) = 0;
        ycrit = 200;
    else
        ycrit = [];
        signal = PLOT.x;
    end
    ids = find(PLOT.sign == -1);
    signal(ids) = signal(ids) .* -1;

    ids = find(PLOT.score < 4 & PLOT.score > 1 & PLOT.sign ~= 0 & tid); 
    hold off;
    if DATA.plot.alttype == 1
    plot(DATA.times(ids),signal(ids),'o');
    else
    plot(ids,signal(ids),'o');
    end
    
    hold on;
    ids = find(respdir == 1  & PLOT.sign ~= 0 & tid & signal ~=0); 
    if DATA.plot.alttype == 1
    plot(DATA.times(ids),signal(ids),'o','markerfacecolor','b');
    else
    plot(ids,signal(ids),'o','markerfacecolor','b');
    end

%red symbol = respdir -1, = choose negative
    ids = find(respdir == -1  & PLOT.sign ~= 0 & tid & signal ~= 0);
    if DATA.plot.alttype == 1
        plot(DATA.times(ids),signal(ids),'ro','markerfacecolor','r');
    else
        plot(ids,signal(ids),'ro','markerfacecolor','r');
    end

%Check Proportion of choices for 0 signal cases, by nominal signal
% plot trials with red and blue circles
    ids = find(PLOT.sign ~= 0 & tid); 
    if isempty(ycrit)
    ycrit = mean(PLOT.y(ids));
    end
    
    did = find(signal == 0 & PLOT.y < ycrit & PLOT.sign ~= 0 & tid);
    pid = find(respdir(did) > 0);
    nid = find(respdir(did) < 0);
    
    if DATA.plot.alttype == 1
        plot(DATA.times(did), signal(did),'o');
        [x,y] = xysmooth(DATA.times(did), respdir(did),DATA.plot.smooth);
        plot(DATA.times(did(pid)), signal(did(pid))-0.02,'ro','markerfacecolor','r');
        plot(DATA.times(did(nid)), signal(did(nid))-0.01,'o','markerfacecolor','b');
    else
        [x,y] = xysmooth(did, respdir(did),DATA.plot.smooth);
        plot(did(pid), signal(did(pid))-0.02,'ro','markerfacecolor','r');
        for j = 1:length(pid)
            k = did(pid(j));
            plot([k-1 k],[signal(k-1) -0.02],'r')
        end
        plot(did(nid), signal(did(nid))-0.01,'o','markerfacecolor','b');
    end
    plot(x,y,'linewidth',2);
    %uid is noominally + 0 signal, includes bad trials 
    uid = find(signal == 0 & PLOT.y >= ycrit  & PLOT.sign ~= 0 & tid);
    pid = find(respdir(uid) > 0);
    nid = find(respdir(uid) < 0);
    if DATA.plot.alttype == 1
    [x,y] = xysmooth(DATA.times(uid), respdir(uid),DATA.plot.smooth);
        plot(DATA.times(uid(pid)), signal(uid(pid))+0.01,'r^','markerfacecolor','r');
        plot(DATA.times(uid(nid)), signal(uid(nid))+0.02,'^','markerfacecolor','b');
    else
    [x,y] = xysmooth(uid, respdir(uid),DATA.plot.smooth);
        plot(uid(pid), signal(uid(pid))+0.01,'r^','markerfacecolor','r');
        plot(uid(nid), signal(uid(nid))+0.02,'^','markerfacecolor','b');
        for j = 1:length(nid)
            k = uid(nid(j));
            plot([k-1 k],[signal(k-1) 0.02],'b')
        end
    end
    plot(x,y,'r','linewidth',2);
    
%ytype is whether nomimal stimulus is + or -
    ytype = PLOT.y < ycrit;
    ytype = (ytype -0.5) * 2;

%plot proportion of nominal stimlus signs for 0 signal
    ids = find(signal == 0 & PLOT.sign ~= 0 & tid);
    if DATA.plot.alttype == 1
    [x,y] = xysmooth(DATA.times(ids), ytype(ids),DATA.plot.smooth);
    else
    [x,y] = xysmooth(ids, ytype(ids),DATA.plot.smooth);
    end
    plot(x,y,'g');
% and proportion of responses for these
    if DATA.plot.alttype == 1
    [x,y] = xysmooth(DATA.times(ids), respdir(ids),DATA.plot.smooth);
    else
    [x,y] = xysmooth(ids, respdir(ids),DATA.plot.smooth.*2);
    end
    plot(x,y,'k');
%also get proportion nominal stimuli in one direction all types
    ids = find(PLOT.sign ~= 0 & tid);
    [x,y] = xysmooth(ids, ytype(ids),DATA.plot.smooth.*2);
    plot(x,y,'g:');
    
    b = hist(DATA.score(ids),[0 1 2 3]);
    title(perfstr);
    GetFigure('PsychHist');
    zid = [uid; did];
    id = find(scores(zid) <=1); %remove BAD trials
    zid = zid(id);
    pid = find(respdir(uid) > 0);
    nid = find(respdir(uid) < 0);
    preres = respdir(zid-1) .* scores(zid-1);
    rshift = respdir(zid) .* respdir(zid-1);
    id = find(scores(zid-1) == 1)  %win
    pshift(1) = mean(rshift(id));
    id = find(scores(zid-1) == -1)  %lose
    pshift(2) = mean(rshift(id));
    fprintf('Winshift %.3f, LoseShift %.3f\n',pshift(1),pshift(2));
    nomsignal = signal;
    nomsignal(uid) = 0.01;
    nomsignal(did) = -0.01;
    if DATA.plot.alttype == 3
    presig = signal(zid-1);
    subplot(2,1,1);
    hold off;
    [a,b] = hist(presig,50);
    maxs(1) = max(a);
    bar(b,a);
    hold on;
    presig = signal(uid(nid)-1);
    [a,b] = hist(presig,50);
    else
    presig = nomsignal(uid-1);
    subplot(2,1,1);
    hold off;
    [a,b] = hist(presig,100);
    maxs(1) = max(a);
    bar(b,a);
    hold on;
    
    presig = nomsignal(uid(nid)-1);
    [a,b] = hist(presig,100);
    bar(b,a,'r');
    title('Signal before 0, nominal +. Red = + choice');

    pid = find(respdir(did) > 0);
    nid = find(respdir(did) < 0);
    presig = nomsignal(did-1);
    subplot(2,1,2);
    hold off;
    [a,b] = hist(presig,100);
    bar(b,a);
    maxs(2) = max(a);
    hold on;
    presig = nomsignal(did(nid)-1);
    [a,b] = hist(presig,100);
    bar(b,a,'r');
    title('Signal before 0, nominal -');
    end
    if length(uid) 
    set(gca,'ylim',[0 max(maxs)]);
    subplot(2,1,1);
    set(gca,'ylim',[0 max(maxs)]);
    end
    end
elseif plottype == 8  %Trial signal/score
    Expts = MakeAllExpts(DATA);
    ExptPsych(Expts, 'sum','sequence','shown');
elseif plottype == 108   %Check 0 Bias
    respdir = zeros(size(PLOT.score));
    %sign = 0 for fixation only trials
    ids = find(PLOT.sign == -1 & PLOT.score == 1);
    respdir(ids) = -1;
    ids = find(PLOT.sign == 1 & PLOT.score == 0);
    respdir(ids) = -1;
    ids = find(PLOT.sign == 1 & PLOT.score == 1);
    respdir(ids) = 1;
    ids = find(PLOT.sign == -1 & PLOT.score == 0);
    respdir(ids) = 1;
    signal = PLOT.x;
    ids = find(PLOT.sign == -1);
    signal(ids) = signal(ids) .* -1;

%Check Proportion of choices for 0 signal cases, by nominal signal
    hold off;
    ids = find(PLOT.sign ~= 0); 
    ycrit = mean(PLOT.y(ids));
    ids = find(signal == 0 & PLOT.y < ycrit & PLOT.sign ~= 0);
    [x,y] = xysmooth(ids, respdir(ids),DATA.plot.smooth);
    plot(x,y,'linewidth',2);
    ids = find(signal == 0 & PLOT.y > ycrit  & PLOT.sign ~= 0);
    [x,y] = xysmooth(ids, respdir(ids),DATA.plot.smooth);
    hold on;
    plot(x,y,'r','linewidth',2);
elseif plottype == 9   %Check 0 Bias
    fprintf('Making Expt...');
    Expts = MakeAllExpts(DATA);
    fprintf('plotting Expt   ');
    epargs = {};
    if DATA.plot.byrw
        epargs = {epargs{:} 'byrw'};
    end
    ExptPsych(Expts{end},'verbose','shown',epargs{:});
    t = get(gca,'title');
    title(sprintf('%s %s',t,perfstr));
    fprintf('\n%s Done at %s\n',perfstr,datestr(now));
elseif plottype == 10   %plot eye pos
    for j = 1:length(DATA.trialvals)
        if isempty(DATA.trialvals(j).lv)
            good(j) = 0;
        else
            good(j) = 1;
        end
    end
    scores = DATA.score(find(good));
    ttimes = DATA.times(find(good));
    trialvals = DATA.trialvals(find(good));
    badt = ismember(scores,[3 7 103 107]);
    bid = find(badt);
    gid = find(~badt);
    plot(ttimes(gid),[trialvals(gid).lh],'o');
    hold on;
    plot(ttimes(gid),[trialvals(gid).lv],'ro');
    plot(ttimes(gid),[trialvals(gid).rh],'go');
    plot(ttimes(gid),[trialvals(gid).rv],'co');
    plot(ttimes(bid),[trialvals(bid).rv],'cx');
    plot(ttimes(bid),[trialvals(bid).lh],'x');
    plot(ttimes(bid),[trialvals(bid).lv],'rx');
    plot(ttimes(bid),[trialvals(bid).rh],'gx');
    legend('LH','LV', 'Rh','RV');
    a(1,1) = prctile([trialvals.lh],1);
    a(1,2) = prctile([trialvals.lh],99);
    a(2,1) = prctile([trialvals.lv],1);
    a(1,2) = prctile([trialvals.lv],99);
    yl(1) = min(a(:,1));
    yl(2) = max(a(:,2));
    yr = diff(yl);
    yl(2) = yl(2)+yr/2;
    yl(1) = yl(1)-yr/2;
    set(gca,'ylim',yl);
elseif plottype == 2
    ids = find(PLOT.score < 4 | ismember(PLOT.score,[51 53])); %ITI all
    [x,y] = xysmooth(PLOT.times(ids(2:end)),diff(PLOT.times(ids)).*(24 *60 * 60),DATA.plot.smooth);
    plot(x,y);
    ids = find(PLOT.score < 3 | PLOT.score == 51); %ITI for completed trials
    [x,y] = xysmooth(PLOT.times(ids(2:end)),diff(PLOT.times(ids)).*(24 *60 * 60),DATA.plot.smooth);
    hold on;
    plot(x,y,'r');
    set(gca,'ylim',[0 10]); %10 sec max
    datetick('x','HH:MM');
    b = hist(PLOT.score(ids),[0 1 2 4]);
    ylabel('Inter-Trial Interval (sec)');
    legend('All Trials','Good Trials');
    title(perfstr);
elseif plottype == 4 %%correct
    id = find(ismember(PLOT.score,[0 1]));
    if DATA.plot.alttype == 2
        plot(smooth(PLOT.score(id),DATA.plot.smooth))
    else
        [x,y] = xysmooth(PLOT.times(id),PLOT.score(id),DATA.plot.smooth);
        plot(x,y);
        [a,b] = xcounts(PLOT.y);
        id = find(a > max(a)/2);
        [a, c] = min(a(id));
        ybest = b(id(c));
        id = find(ismember(PLOT.score,[0 1]) & PLOT.y == ybest);
        if length(id) > DATA.plot.smooth
        hold on;
        [x,y] = xysmooth(PLOT.times(id),PLOT.score(id),DATA.plot.smooth);
        plot(x,y,'r');
            datetick('x','HH:MM');
        legend('All',sprintf('y = %.2f',ybest));
        end
    end
    title(perfstr);
end
if PLOT.exptover
    set(gca,'color','k');
else
    set(gca,'color','w');
end    

function DATA = ReadFile(DATA,fname)
%File format
% R%d  0 = wrong
%      1 = correct
%      2 = late/foul
%      3 = BadFix
%      >4 special lines with stimulus/expt info
%      7 Badfix caused by a microsaccade

startday = 0;

if ~exist(fname,'file')
    fprintf('Can''t read %s\n',fname);
    DATA = [];
    return;
end
oname = fname;
if DATA.online
fname = [oname '.txt'];
delete(fname);
try
copyfile(oname, fname);
catch
    fname = oname;
end
end

fid = fopen(fname,'r');
if fid > 0
    a = fgetl(fid);
    id = strfind(a,' ');
    frewind(fid);
    if length(id) <= 7
        a = textscan(fid,'R%f %[^=]=%f %2s=%f %2s=%f %f %f %f','bufsize',2048);
        nxval = 0;
    elseif ismember(length(id),[8 9]) %later version, gives seed or id
        a = textscan(fid,'R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f');
        nxval = 2;
    elseif ismember(length(id),[12 13]) %later version, gives seed or id
        a = textscan(fid,'R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
        nxval = 6;
    end
    fclose(fid);
else
    fprintf('Cant Read %s\n',fname);
    DATA.score = [];
    return;
end
id = regexp(DATA.filename,'[0-9][0-9][A-Z]');
if length(id)
    startday = datenum(DATA.filename(id(1):id(1)+8));
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
d = dir(fname);
DATA.filedate = d.datenum;
DATA.score = a{1};
for j = 1:nxval
    x = a{9+j*2};
    id = strfind(x,'(');  %fields with (..) encode extra values. removethe braces from field names
    if ~isempty(id)
        x = regexprep(x,'\(.*\)','');
    end
    y = a{9+j*2+1};
    if DATA.verbose && j == 1
        fprintf('%d trials,', length(x));
    end
    for k = 1:min([length(x) length(y)])
    DATA.trialvals(k).(x{k})=y(k);
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
DATA.times(find(DATA.score == 6)) = 0;
DATA.times(find(DATA.score == 4)) = 0;
if DATA.times(1) == 0
    id = find(DATA.times > 0);
    DATA.times(1) = DATA.times(id(1));
end
id = find(DATA.times == 0);
while length(id) %can get two in a row.
    DATA.times(id) = DATA.times(id-1);
    id = find(DATA.times == 0);
end
DATA.times(id) = DATA.times(id-1);
steps = diff(DATA.times);
jumps = find(steps < 0 & DATA.times(2:end) > 0)+1;
steps = cumsum(steps(jumps-1));
jumps(end+1) = length(DATA.times)+1;
for j =1:length(jumps)-1
    DATA.times(jumps(j):jumps(j+1)-1) = DATA.times(jumps(j):jumps(j+1)-1) - steps(j);
end
if DATA.verbose
    fprintf('%d jumps,', j);
end

DATA.rwszs = a{10};
DATA.DURS = a{9};
DATA.readtime = now;

tid = find(DATA.score == 4); %start expts
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

oid = strmatch('dd',a{15}(tid));
if length(oid)
    xs = a{16}(tid(oid));
    DATA.Blockid.dd = tid(oid);
    DATA.Block.dd = xs;
    DATA.Stimvals.dd = mean(xs);  
end

oid = strmatch('c2',a{17}(tid));
if length(oid)
    xs = a{18}(tid(oid));
    DATA.Blockid.c2 = tid(oid);
    DATA.Block.c2 = xs;
    DATA.Stimvals.c2 = mean(xs);  
end


DATA.Stimvals.wi = round(DATA.Stimvals.sz * 10)/10;
DATA.Stimvals.hi = round(DATA.Stimvals.sz * 10)/10;

oid = strmatch('ve',a{2}(tid));
stid = tid(oid);

if DATA.verbose
    fprintf('%d ves,', length(oid));
end
if length(oid)
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
DATA.Block.bh = bhs;DATA.Stimvals.bh = mean(bhs);  
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
if ~DATA.plot.noplot
    set(DATA.toplevel,'UserData',DATA);
    PlotData(DATA);
else
    DATA = PlotData(DATA);
end

function DATA = InitInterface(DATA)


    name = DATA.name;
    fname = DATA.fcn;
    scrsz = get(0,'Screensize');
    wsc = scrsz(3) /1280;  %scale everything as is 1280 wide
    if scrsz(3) > 2000
        wsc = 1.5;
    end
    size(1) = 380 * wsc;
    size(2) = 200 * wsc;
    listtag = [DATA.tag 'List'];
   
    HSPACE = 3;
    VSPACE = 2;

    cw = 9;
    ch = 18*wsc + VSPACE;
    nlines = 4;
    cntrl_box = figure('Position', [10 scrsz(4)-220*wsc 300*wsc 200*wsc],...
        'NumberTitle', 'off', 'Tag',DATA.tag,'Name',DATA.name);
    DATA.toplevel = cntrl_box;
    
    if( ~isempty(DATA.filename))
        lst = uicontrol(gcf, 'Style','text','String', DATA.filename,...
            'Position',[10 10 size(1) ch]);
        DATA.gui.filename = lst;
        
    end

    
    bp(1) = HSPACE; bp(3) = 10*cw; bp(2) = size(2)-ch; bp(4) = 22;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [fname '(''reload'',''Tag'',''' DATA.tag ''')'],...
        'String', 'reload', 'Position', bp);
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [fname '(''next'',''Tag'',''' DATA.tag ''')'],...
        'String', '>>', 'Position', bp);
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [fname '(''prev'',''Tag'',''' DATA.tag ''')'],...
        'String', '<<', 'Position', bp);
    
 %New row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw * wsc;
    uicontrol(gcf,'Style', 'checkbox','String', 'AutoPlot', 'Tag', 'AutoPlot', 'Callback', [fname '(''update'',''Tag'',''' DATA.tag ''')'],...
        'Position', bp);
   
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'FlipPlot', 'Tag', 'Flip', 'Callback', [fname '(''update'',''Tag'',''' DATA.tag ''')'],...
        'Position', bp);

    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'checkbox','String', 'byrw', 'Tag', 'byrw', 'Callback', [fname '(''update'',''Tag'',''' DATA.tag ''')'],...
        'Position', bp);

    bp(2) = bp(2) - ch;
    bp(3) = 5 * cw * wsc;
    uicontrol(gcf,'Style', 'text','String','Plot','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','PSF|ITI|Rewards|%Correct|Types|Total Reward|Trials|Bias|Expt|EMpos', ...
        'Callback', [fname '(''setplot'',''Tag'',''' DATA.tag ''')'], 'Tag','plottype',...
        'position',bp,'value',1);
    
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','Time|Trial#', ...
        'Callback', [fname '(''update'',''Tag'',''' DATA.tag ''')'], 'Tag','AltPlot',...
        'position',bp,'value',DATA.plot.alttype);

    bp(1) = bp(1)+bp(3)+HSPACE;
    bp(3) = 5 * cw * wsc;
    uicontrol(gcf,'Style', 'text','String','Rws','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','All|1|2|3|4|5|6', ...
        'Callback', [fname '(''update'',''Tag'',''' DATA.tag ''')'], 'Tag','rwtype',...
        'position',bp,'value',1);

    
    bp(2) = bp(2) - ch;
    bp(1) = HSPACE;
    uicontrol(gcf,'Style', 'text','String','Smooth','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',DATA.plot.smooth),'Position', bp,'Tag','Smooth','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','Round X','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%f',DATA.plot.round(1)),'Position', bp,'Tag','xround','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','Y','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%f',DATA.plot.round(2)),'Position', bp,'Tag','yround','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);
        
    bp(2) = bp(2) - ch;
    bp(1) = HSPACE;
    uicontrol(gcf,'Style', 'text','String','Time from','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',DATA.plot.timerange(1)),'Position', bp,'Tag','TStart','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','End','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%f',DATA.plot.timerange(2)),'Position', bp,'Tag','TEnd','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);

    
  bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'Style', 'text','String','Nmin','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.0f',DATA.plot.nmin),'Position', bp,'Tag','Nmin','Callback', ...
	    [fname '(''Update'',''Tag'',''' DATA.tag ''')'],'Backgroundcolor',[1 1 1]);
  bp(1) = bp(1)+bp(3)+HSPACE;

  hm = uimenu(gcf,'Label','File');
  uimenu(hm,'Label','Close','Callback',{@PsychMon, 'close', 'name', name});
  uimenu(hm,'Label','Run Another Block','Callback',{@PsychMon, 'newblock'});
    set(gcf,'Menubar','none');
    
  set(cntrl_box,'UserData',DATA);
  
  h = GetFigure(DATA.figtag);
  hm = uimenu(gcf,'Label','Run');
  uimenu(hm,'Label','Close','Callback',{@PsychMon, 'close', 'name', name});
  uimenu(hm,'Label','Run Another Block','Callback',{@PsychMon, 'newblock'});
  uimenu(hm,'Label','Reload','Callback',{@PsychMon, 'reload'});
  set(hm,'UserData',DATA.toplevel);  

  
function update(DATA,varargin)
str = get(findobj(DATA.toplevel,'Tag','Smooth','Parent',DATA.toplevel),'string');
if(str) 
     DATA.plot.smooth = sscanf(str,'%f');
end
str = get(findobj(DATA.toplevel,'Tag','xround'),'string');
if(str) 
     DATA.plot.round(1) = sscanf(str,'%f');
end
str = get(findobj(DATA.toplevel,'Tag','yround'),'string');
if(str) 
     DATA.plot.round(2) = sscanf(str,'%f');
end
str = get(findobj(DATA.toplevel,'Tag','TStart'),'string');
if(str) 
     DATA.plot.timerange(1) = sscanf(str,'%f');
end
str = get(findobj(DATA.toplevel,'Tag','TEnd'),'string');
if(str) 
     DATA.plot.timerange(2) = sscanf(str,'%f');
end

str = get(findobj(DATA.toplevel,'Tag','Nmin'),'string');
if(str) 
     DATA.plot.nmin = sscanf(str,'%f');
end

DATA.plot.rwtype = get(findobj(DATA.toplevel,'Tag','rwtype'),'value');
DATA.plot.alttype = get(findobj(DATA.toplevel,'Tag','AltPlot'),'value');

DATA.plot.flip = get(findobj(DATA.toplevel,'Tag','FlipPlot','Parent',DATA.toplevel),'value');
DATA.plot.byrw = get(findobj(DATA.toplevel,'Tag','byrw','Parent',DATA.toplevel),'value');
auto = DATA.plot.autoplot;
DATA.plot.autoplot = get(findobj(DATA.toplevel,'Tag','AutoPlot','Parent',DATA.toplevel),'value');
if DATA.plot.autoplot & auto
    PlotData(DATA);
end
if isvalid(DATA.timerobj)
if auto == 1 && DATA.plot.autoplot == 0 %turned off
    stop(DATA.timerobj);
    DATA.verbose = 0;
end
if auto == 0 && DATA.plot.autoplot == 1 & strcmp('off',get(DATA.timerobj,'Running'))  %turned on
    start(DATA.timerobj);
    DATA.verbose = 1;
end
end
set(DATA.toplevel,'UserData',DATA);


function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
    close(it);
end

  
function [counts, u] = xcounts(x)

u = unique(x);
for j = 1:length(u)
    counts(j) = sum(x == u(j));
end

