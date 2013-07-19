function [Expt, Expts, AllData] = PlaySpkFile(name, varargin)

SpkDefs;
playspikes = 0;
defaults.fz = 96;
starttrial = 0;
idxfile = [];
thecluster = 1;
mkidx = 0;
if ischar(name)
    load(name);
    Text = Ch30;
    Spks = Ch5;
    Events = Ch31;
    idxfile = strrep(name,'.mat','idx.mat');
    if idxfile & ~exist(idxfile,'file')
        mkidx = 1;
    else
        load(idxfile);
    end
end

if isstruct(name)
    if isfield(name,'title') & isfield(name,'values') &... 
            size(name.values,2) == 46
        Spks = name;
        name = 'Unnamed'
    end
end

j = 1;
while j <= nargin-1
    vg = varargin{j};
    if isstruct(vg) 
        if isfield(vg,'text')
            Text = vg;
        elseif isfield(vg,'codes');
            Events = vg;
        end
    elseif strncmpi(vg,'Defaults',4)
        j = j+1;
        defaults = varargin{j};
        if isfield(defaults,'starttrial')
            starttrial = defaults.starttrial;
        end
    elseif strncmpi(vg,'name',4)
        j = j+1;
        name = varargin{j};
    elseif strncmpi(vg,'cluster',2)
        j = j+1;
        thecluster = varargin{j};
    elseif strncmpi(vg,'play',4)
        playspikes = 1;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    end
    j = j+1;
end

Events.times = Events.times * 10000;
Spks.times = Spks.times * 10000;
Text.times = Text.times * 10000;
Header.Name = name;
if idxfile & exist(idxfile,'file')
    load(idxfile);
    AllData.Spikes = Spks;
    AllData.Events = Events;
    return;
end
nt = 0;
nx = 0;
Expts = [];

%id = find(ismember(Spks.codes(:,1),thecluster));

frametimes = [];
bstimes = [];

if exist('Ch7','var')
    frametimes = Ch7.times * 10000;
end
if exist('Ch8','var')
%    id = find(Ch8.level == 1);
    bstimes = Ch8.times(Ch8.level == 1) * 10000;
end
tic;
opid = strmatch('op',Text.text);
for j = 1:length(opid)
    str= sscanf(Text.text(opid(j),:),'op%d');
    if ~isempty(str)
        storing(j) = str;
    else
        storing(j) = 0;
    end
end
storing = bitand(storing,8192);
storeonoff = 1+find(abs(diff(storing)) > 0); 
if isempty(storeonoff)
    storeonoff = [1 1];
end

ids = strmatch('fz',Text.text);
if ids
fzs = textscan(Text.text(ids,:)','fz%n');
framerate = mean(fzs{1});
%    framerate = sscanf(Text.text(ids(1),:),'fz%n');
else
    if isfield(defaults,'fz')
        framerate = defaults.fz;
    else
        framerate = 96;
    end
end

ids = strmatch('nf',Text.text);
if ids
    fstr = Text.text(ids,1:5);
    fstr(fstr==0) = ' ';
    fzs = textscan(fstr','nf%n');
    nomdur = prctile(fzs{1},90) .* 10000/framerate;
end
instim = 0;
inexpt = 0;
nextonoff = 2;
storestate = storeonoff(1);
tonoff = Text.times(opid(storeonoff(nextonoff)));
for j = 1:size(Events.codes,1)
    if Events.codes(j,1) == FRAMESIGNAL
        nowt = Events.times(j);
        while nowt > tonoff & nextonoff < length(storeonoff)
            storestate = storing(storeonoff(nextonoff));
            nextonoff = nextonoff+1;
            tonoff = Text.times(opid(storeonoff(nextonoff)));
        end
        if storestate
            nt = nt+1;
        Trials(nt).Start = Events.times(j);
        if nt > 1 & isempty(Trials(nt-1).End)
            fprintf('Missing end Trial %d\n',nt-1);
        end
        Trials(nt).Trial = nt+starttrial;
        instim = 1;
        id = find(bstimes < Events.times(j));
        if ~isempty(id)
            if(bstimes(id(end)+1) - Trials(nt).Start < 10) %< 1ms to next = probably early
                Trials(nt).Start = bstimes(id(end));
            else
                Trials(nt).Start = bstimes(id(end));
            end
            Trials(nt).stored = storestate;
            if ~isempty(frametimes)
            id = find(frametimes > Trials(nt).Start);
            if ~isempty(id)
                Trials(nt).delay = frametimes(id(1)) - Trials(nt).Start;
                Trials(nt).Start = frametimes(id(1));
                Trials(nt).TrueStart = 1;
            else
                Trials(nt).Start = 0;
                Trials(nt).TrueStart = 0;
            end
            end
        else
            Trials(nt).Start = Events.times(j);
            Trials(nt).TrueStart = 0;
        end
    end
    elseif Events.codes(j,1) == ENDSTIM & storestate & nt
        Trials(nt).End = Events.times(j);
        Trials(nt).Result = 1;
        if (Trials(nt).End - Trials(nt).Start) < 1000
            instim = 0;
        end
        instim = 0;
    elseif Events.codes(j,1) == ENDTRIAL & storestate
        if instim
            fprintf('End Trial without End stim: %d (%.2f)\n',nt-1,Events.times(j)/10000);
            Trials(nt).End = Events.times(j);
            if Events.times(j) - Trials(nt).Start > 0.9 * nomdur
                Trials(nt).Result = 1;
            else
                Trials(nt).Result = 0;
            end
        end
    elseif Events.codes(j,1) == BADFIX & storestate
        Trials(nt).End = Events.times(j);
        Trials(nt).Result = 0;
    elseif Events.codes(j,1) == STARTEXPT & storestate
        if inexpt
            Expts(nx).end = Events.times(j);
            Expts(nx).lasttrial = nt;
        end
        nx = nx+1;
        Expts(nx).start = Events.times(j);
        Expts(nx).firsttrial = nt+1;
        inexpt = 1;
    elseif Events.codes(j,1) == ENDEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        inexpt = 0;
        Expts(nx).result = ENDEXPT;
    elseif Events.codes(j,1) == CANCELEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        Expts(nx).result = CANCELEXPT;
        inexpt = 0;
    end
end

if inexpt
    Expts(nx).lasttrial =nt;
end
toc;
tic;
trial = 1;
k = 1;
ix = 1;
% this should work, but need to check with some xxx= strings.
%aText.text = cellstr(Text.text);
%aText.times = Text.times;
%aText.codes = Text.codes;
%ids = strmatch('xxx=',aText.text);
trynew = 1;
if trynew
    tic;
aText.text = cellstr(Text.text);
aText.times = Text.times;
aText.codes = Text.codes;
ids = strmatch('xxx=',aText.text);
    for j = length(ids):-1:1
        aText.text{ids(j)-1} = [aText.text{ids(j)-1} aText.text{ids(j)}];
    end
    fprintf('New method takes %.2f\n',toc);
end
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
    else
    aText.text{k} = Text.text(j,1:end-1);
    aText.times(k) = Text.times(j);
    aText.codes(k,:) = Text.codes(j,:);
    end
    if length(aText.text{k})
        k = k+1;
    end
end
for j = 1:length(aText.text)
    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);
%    id = find(txt == 0);
%    txt = txt(1:id(1)-1);
    if ~isempty(txt);
        ss = txt(1:2);
        val = txt(3:end);
        if strncmp(txt,'st',2)
            Stimulus.st = strmatch(val, stimnames,'exact');
            Stimulus.st = Stimulus.st -1;
        elseif strncmp(txt,'mtop=op',7)
            Stimulus.OptionCode = txt(8:end);
        elseif strncmp(txt,'mtrS=',5)
            Stimulus.Stimseq = sscanf(txt(6:end),'%d');
            if ~instim & trial > 1
                Trials(trial-1).Stimseq = Stimulus.Stimseq;
            end
        elseif strncmp(txt,'mtei=',5)
            Expts(ix).e1vals = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'mte2=',5)
            Expts(ix).e2vals = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'op',2)
            Stimulus.op = sscanf(txt(3:end),'%f');
            if isempty(Stimulus.op)
                fprintf('Missing op stim %d\n',trial);
                Stimulus.op = 0;
            end
        else
%            if strmatch(ss,{'0' '1' '2' '3' '4'})
            if regexp(ss,'^[0-9]')
                ss = ['x' ss];
            end
            [Stimulus.(ss), ok] = sscanf(val,'%f');
            if ~ok
                Stimulus.(ss) = val;
            end
        end
    end
    if t > Trials(trial).End & instim
        fn = fieldnames(Stimulus);
        for k = 1:length(fn)
            Trials(trial).(fn{k}) = Stimulus.(fn{k});
        end
        trial = trial+1;
        instim = 0;
        if trial > length(Trials)
            break;
        end
    elseif t > Trials(trial).Start
        instim = 1;
    end
end
toc;

trial = 1;
instim = 0;
colors = mycolors;
lastspk = 1;
maxspk = 1;
tpause = 0;
if playspikes
    GetFigure('SpikeV');
end
postdur = 500;

tic;
for trial = 1:length(Trials)
    if isempty(Trials(trial).End)
        Trials(trial).End = Trials(trial).Start + nomdur;
    end
    spkids = find(Spks.times > Trials(trial).Start & Spks.times < Trials(trial).End);
    Trials(trial).Spikes = (Spks.times(spkids) - Trials(trial).Start);
    if spkids
        spkid(trial,:) = [spkids(1) spkids(end)];
        Trials(trial).Cluster = Spks.codes(spkids,1);
        frames = [];
        if frametimes
            for k = spkids
                id = find(frametimes < Spks.times(k));
                if id
                    frames = [frames frametimes(id(end))-Spks.times(k)];
                end
            end
        end
            Trials(trial).Frames = frames;
    else
        spkid(trial,:) = [0 0];
        Trials(trial).Frames = [];
        Trials(trial).Cluster = [];
    end
    if playspikes
        for spk = spkids; 
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            drawnow;
            hold on;
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            if tpause
                pause(tpause);
            end
        end
        %        GetFigure('SpikeXY');
        subplot(2,1,1);
        plot(energy(lastspk:maxspk),vw(lastspk:maxspk),'.');
        hold on;
        %       GetFigure('SpikeV');
        subplot(2,1,2);
        hold off;
    end
  
end
fprintf('Spikes  take %.3f\n',toc);

tic;
for spk = 1:size(Spks.times)
    t = Spks.times(spk);
    if t > Trials(trial).End+postdur & instim
        Trials(trial).Spikes = (Spks.times(lastspk:maxspk) - Trials(trial).Start);
        spkid(trial,:) = [lastspk, maxspk];%
%        Trials(trial).Frames = frames;
        Trials(trial).Cluster = Spks.codes(lastspk:maxspk,1);
        trial = trial+1;
        instim = 0;
        lastspk = spk;
        elseif t > Trials(trial).Start
            instim = 1;
        end
        if instim & playspikes
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            hold on;
        end
        maxspk = spk;
        if trial > length(Trials)
            break;
        end
end
fprintf('Spikes Took %.3f\n',toc);
%Trials = Trials(2:end);
AllData.Events = Events;
AllData.Text = aText;
AllData.Spikes = Spks;
stored = find(bitand([Trials.op],8192));
% bitand 8192 identifies if storage was on;
%Trials = Trials(stored);

%Expts #20 seems to be empty, but storage was on...
Expt.Trials = Trials;
Expt.Spkid = spkid;

Header.Name = name;

if isempty(Expts) %nothing in this file
    return;
end
if mkidx
    fprintf('Saving index %s',idxfile);
    save(idxfile,'Expt');
end
Expts = SortExpts(Expts, Trials, Header, thecluster);
if mkidx
    save(idxfile,'Expt','Expts');
end

function Expts = SortExpts(AllExpts, AllTrials, Header, thecluster)
SpkDefs;

frameperiod = 167;
fn = fieldnames(AllTrials);
% make a list of fileds that are NOT automatically set 
ids = strmatch('Spikes',fn);
ids = [ids strmatch('OptionCode',fn)];
ids = [ids strmatch('Start',fn)];
ids = [ids strmatch('End',fn)];
ids = [ids strmatch('Result',fn)];
ids = [ids strmatch('op',fn)];
ids = [ids strmatch('Stimseq',fn)];
ids = [ids strmatch('Cluster',fn)];
ids = [ids strmatch('Frames',fn)];
fn = fn(setdiff([1:length(fn)],ids));
%fn is now a list of fields that are set for each trial automatically.
nexpts = 1;
for nx = 1:length(AllExpts)
    
    clear Trials;
    a = AllExpts(nx).firsttrial:AllExpts(nx).lasttrial;
    igood = find([AllTrials(a).Result] > 0);
    nt = length(igood);
    if nt> 3 & igood(1) < length(AllTrials)
    [Trials(1:nt).Start] = deal(AllTrials(igood).Start);
    [Trials(1:nt).TrialStart] = deal(AllTrials(igood).Start);
    [Trials(1:nt).End] = deal(AllTrials(igood).End);
    [Trials(1:nt).op] = deal(AllTrials(igood).op);
    [Trials(1:nt).Trial] = deal(AllTrials(igood).Trial);
    for nf = 1:length(fn)
        nv = unique([AllTrials(igood).(fn{nf})]);
        if isnumeric(nv) & length(nv) > 1
            [Trials(1:nt).(fn{nf})] = deal(AllTrials(igood).(fn{nf}));
            Expt.Stimvals.(fn{nf}) = prctile([AllTrials(igood).(fn{nf})],50);
        else
            Expt.Stimvals.(fn{nf}) = AllTrials(a(1)).(fn{nf});
        end
    end
    duration = mean([Trials.End] - [Trials.Start]);
    et = Expt.Stimvals.et;
    if isfield(AllExpts(nx),'e1vals')
        Expt.e1vals = AllExpts(nx).e1vals;
    end
    if isfield(AllExpts(nx),'e2vals')
        Expt.e2vals = AllExpts(nx).e2vals;
    end
    allev = [];
    for k = 1:nt
        sid = find(ismember(AllTrials(igood(k)).Cluster, thecluster));
        Trials(k).Spikes = round(AllTrials(igood(k)).Spikes(sid));
        Trials(k).Count = sum(find(Trials(k).Spikes > 500 & Trials(k).Spikes < duration+500));
        Trials(k).OptionCode = AllTrials(igood(k)).OptionCode;
        if ~isempty(AllTrials(igood(k)).Frames)
            Trials(k).Frames = AllTrials(igood(k)).Frames(sid);
        end
        if bitand(LMONOC,Trials(k).op)
            Trials(k).me  = -1;
        elseif bitand(RMONOC,Trials(k).op)
            Trials(k).me  = 1;
        else
            Trials(k).me  = 0;
        end
        if isfield(AllTrials(igood(k)),'Stimseq')
            ev = Expt.e1vals(AllTrials(igood(k)).Stimseq+1);
            Trials(k).st = ones(size(ev)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(ev)) * Expt.Stimvals.ce;
            Trials(k).me = zeros(size(ev)) * Trials(k).me;
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod;
            Trials(k).End = Trials(k).Start + frameperiod;
            Trials(k).(et) = ev;
            eb = Expt.e2vals(AllTrials(igood(k)).Stimseq+1);
            Trials(k).(Expt.Stimvals.e2)= eb;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            Trials(k).ce(find(ev == IUNCORR)) = 0; %% uncorr
            Trials(k).me(find(ev == ILEFTMONOC)) = -1; 
            Trials(k).me(find(ev == IRIGHTMONOC)) = 1;
        end
    end
    Expt.Trials = Trials;
    Expt.Header = Header;
    if strmatch(Expt.Stimvals.e2, 'e0')
        expname = [Expt.Stimvals.et];
    else
        expname = [Expt.Stimvals.et 'X' Expt.Stimvals.e2];
    end
    if isempty(expname)
        Expt.Header.expname = 'None';
    else
        Expt.Header.expname = expname;
    end
    Expts{nexpts} = Expt;
    nexpts = nexpts+1;
    end
end
    
  