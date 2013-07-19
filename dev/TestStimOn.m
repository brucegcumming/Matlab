function Result = TestStimOn(name, varargin)
%TestStimOn
% Try to find a better way to detect trials in a matlab file.

% [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%
% Builds/retnes lists of trials/expts from a matlab file
% APlaySpkFile(name, 'relist')  Re-builds the index filt
% APlaySpkFile(name, 'setprobe', n)  sets which probes data is loaded
% initially
% APlaySpkFile(name, 'setprobe', -1) rebuilds probe list


SpkDefs;
playspikes = 0;
defaults.fz = 96;
starttrial = 0;
onlinedata = 0;
idxfile = [];
Expts = [];
Expt = [];
AllData = [];
s2version = 0;
testframe = 0;
findtrial = 0;
spkch = 'Ch5';
stimch = 'Ch8';
framechname = 'Ch7';
mainsname = 'Ch24';
dstimch = 'Ch18';
ustimmarkname = '';
ustimmarkch = [];
mainsch = [];
stimlvl = [];
logfid = 0;
UstimV = 0;
setprobe = 1;
nospikes = 0;
nerr = 0;
preperiod=2000;
postperiod = 2000;
savedvdt = 0;
dvfile = [];
Raw = [];
timeoffset = 0;
plottype = 0;
ignoreSpikeO = 0;  %% for files with full voltage records elsewhere, ignore these
j = 1;
while j <= length(varargin)
    vg = varargin{j};
    if ischar(vg) & strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'setprobe',4)
        j = j+1;
        setprobe = varargin{j};
    elseif strncmpi(vg,'findprobes',6)
       mkidx = 1;
       setprobe = -1; %force relisting of probes
    elseif strncmpi(vg,'rfs',3)
       nospikes = 1;
    elseif strncmpi(vg,'bysuffix',7)
       ignoreSpikeO = 2;
    elseif strncmpi(vg,'scatter',7)
        plottype = 1;
    elseif strncmpi(vg,'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    
    j = j+1;
end
    
thecluster = 1;
mkidx =0;

if iscell(name)
    for j = 1:length(name)
        if isfield(name{j},'Trials')
        T = name{j}.Trials;
        if ~isfield(T,'Expts') && isfield(T,'ExptStart')
            for k = 1:length(T.ExptStart)
                T.Expts(k).Start = T.ExptStart(k);
                T.Expts(k).End = T.ExptEnd(k);
            end
        end
        if length(name{j}.badisi) > 0
            fprintf('%s %d Badisis\n',name{j}.name,length(name{j}.badisi));
        end
        nt = [];
        for e = 1:length(T.Expts)
           id = find([T.Expts.End] > T.Expts(e).Start);
           exptend(e) = T.Expts(id(1)).End;
            ts{e} = find(T.Start > T.Expts(e).Start & T.Start < exptend(e));
            nt(e) = length(ts{e});
        end
        [a,e] = max(nt);
        id = find(T.Result == 1 & T.Start > T.Expts(e).Start & T.Start < exptend(e));
        
        if isempty(id)
            dur(j) = NaN;
            durdiff(j) = NaN;
            sdiff(j) = NaN;
            dsd(j) = NaN;
        else
            durs = T.End(id)-T.Start(id);
            dsd(j) = std(durs);
            durdiff(j) = max(durs)-min(durs);
            dur(j) = mean(durs);
            sdiff(j) = diff(name{j}.startcounts);
        end
        end
    end
    hold off;
    for j  = 1:length(dur)
    plot(dur(j),dsd(j),'o','buttondownfcn',{@HitPopPoint, j});
    hold on;
    end
    id = find(sdiff ~= 0);
    for j = 1:length(id)
        plot(dur(id(j)),dsd(id(j)),'ro','buttondownfcn',{@HitPopPoint, id(j)});
    end
    id = find(durdiff > 150); %more than one frame
    for j = 1:length(id)
        plot(dur(id(j)),dsd(id(j)),'go','markerfacecolor','g','buttondownfcn',{@HitPopPoint, id(j)});
    end
    DATA.dsd = dsd;
    DATA.dur = dur;
    DATA.sdiff = sdiff;
    DATA.Result = name;
    DATA.plottype = plottype;
    DATA.toplevel = gcf;
    Result = DATA;
    set(gcf,'UserData',DATA);
    return;
elseif isfield(name,'Trials')
    GetFigure('TestStim');
    T = name.Trials;
    durs = T.End-T.Start;
    gid = find(T.Result > 0);
    if plottype == 1
        plot(T.nf(gid),durs(gid),'o');
    else
    hist(durs(gid));
    end
    return;
end
if ischar(name) & exist(name,'dir')
    Result = {};
    list = TreeFind(name,'name','[0-9]\.[0-9]*\.mat');
    for j = 1:length(list)
        fprintf('%s:\n',list{j});
        Result{j} = TestStimOn(list{j});
        if ~isempty(Result{j})
        Result{1}.Summary.means(j) = mean(Result{j}.bsdelay);
        Result{1}.Summary.latebs(j) = sum(Result{j}.bsdelay > 300);
        Result{1}.Summary.earlybs(j) = sum(Result{j}.bsdelay < -300);
        Result{1}.Summary.badisi(j) = length(Result{j}.badisi);
        end
    end
    return;
end
Result = [];
if ischar(name)
    if ~exist(name,'file')
        fprintf('No file %s\n',name);
        return;
    end
    
    np = 0;
    nlfp = 0;
    nspkt = 0;
    probes = [];
    logname = strrep(name,'.mat', '.log');
%    fprintf('Log %s\n',logname);
   logfid = 0;
   Oprobe = 0;
   if onlinedata
       %        oname = strrep(name,'online','online2');
       oname = strrep(name,'/Expt','A/Expt');
       oname = regexprep(oname,'(\.[0-9]*.mat)','A$1');
       if exist(oname,'file')
           af = load(oname);
           f = fields(af);
       else
           f = {};
       end
       for j = 1:length(f)
           if ~isempty(regexp(f{j},'Ch[0-9]*'))
               ch = af.(f{j});
               if strncmpi(ch.title,'Spike',5)
                   np = np+1;
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               elseif strncmpi(ch.title,'4Trode',5)
                   np = np+1;
                   probe = sscanf(ch.title,'4Trode%d');
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               end
           end
       end
   end
%    fprintf('Reading %s\n',name);
   mkmatver = 0;

 
   if ignoreSpikeO == NaN % don't need this any more, for online at least
       oname = regexprep(name,'.([0-9]*.mat)','A.$1');
       if exist(oname,'file')
           load(oname);
       end
       avars = who('Ch[0-9]*');
       for j = 1:length(avars)
           eval([avars{j} 'A = ' avars{j} ';']);
           clear(avars{j});
       end

   end

   load(name);
    
    vars = who('Ch*');
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
            if chn > 400  %a memory/extra offline channel
            elseif strncmpi(ch.title,'SpikeO',6) && ignoreSpikeO;
            elseif strncmpi(ch.title,'Spike',5)
                np = np+1;  
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                probes(np).traces = ch.traces;
                probes(np).source = 1;


            elseif strncmpi(ch.title,'4Trode',5)
                np = np+1;  
                probe = sscanf(ch.title,'4Trode%d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                if probe == setprobe(1)
                    Chspk = ch;
                    spkch = 'Chspk';
                end
                probes(np).traces = ch.traces;
                probes(np).source =1;

            elseif strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            elseif strncmpi(ch.title,'uStim',5)
                np = np+1;  
                probe = 100;
                UstimV = ch;
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'StimChan',8) % stim change detector
                dstimch = vars{j};
                stimchange = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framechname = vars{j};
                framech = ch;
   %             fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            elseif strncmpi(ch.title,'DigMark',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            end
        end
    end
  vnames = {'Ch30' 'Ch31'  framechname stimch};
    vlabels = {'Text' 'Events'  'Frames' 'Stim ON/OFF'};
    for j = 1:length(vnames)
        missing(j) = ~exist(vnames{j},'var');
    end
    if sum(missing)
        msgbox(sprintf('%s Missing %s',name,vlabels{find(missing)}),'Error!!');
        fprintf('%s Missing %s\n',name,vlabels{find(missing)});
        if logfid
            fprintf(logfid, '%s Missing %s\n',name,vlabels{find(missing)});
        
        fclose(logfid);
        end
        return;
    end
    Text = Ch30;
    Events = Ch31;
  
end



forcefix = 0;
argon = {};
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
    elseif strncmpi(vg,'fixlfp',5)
        if strncmpi(vg,'fixlfpforce',8)
            forcefix = 1;
        end
        lfpfile = strrep(name,'.mat','A.lfp.mat');
        fixfile = strrep(name,'.mat','.lfp.mat');
% with 8 channels, everything is in one file - there is no 'A.lfp.mat' 
        if ~exist(lfpfile,'file') && exist(fixfile,'file')
            lfpfile = fixfile;
            forcefix = 1;
        end
        if exist(lfpfile,'file') && (~exist(fixfile,'file') || forcefix) ...
                && exist(mainsname,'var')
            load(lfpfile);
            if isfield(LFP.Header,'MainsNoise')
                fprintf('%s Already fixed\n',lfpfile);
                if logfid
                fprintf(logfid,'%s Already fixed\n',lfpfile);
                end
            else
            [LFP, avgs, NoiseAmp] = FixLFPMains(LFP,mainsch.times .* 10000);
            LFP.Header.amps = LFPGains(LFP);
            a = LFP.Header.amps ./ max(LFP.Header.amps);
            nch = sum(LFP.Header.chanlist > 0);
            if std(a(find(a > 0.1))) > 0.2 && nch <= 8  %% 8 channel probe tends to have mixed LFP gain
                LFP.Header.needscale = 1;
            else
                LFP.Header.needscale = 0;
            end
            save(fixfile,'LFP','NoiseAmp');
           if logfid > 0
               fprintf(logfid, '%s Fixed LFP in %s\n',datestr(now),fixfile);
           end
            end
        else
            fprintf('No LFP file %s\n',lfpfile);
           if logfid > 0
               fprintf(logfid, 'No LFP file %s\n',lfpfile);
           end
            
        end
        if strncmpi(vg,'fixlfponly',8)
            return;
        end
            
    elseif strncmpi(vg,'findtrial',5)
        j = j+1;
        argon = {argon{:} varargin{j-1} varargin{j}};
        findtrial = varargin{j};

    elseif strncmpi(vg,'name',4)
        j = j+1;
        name = varargin{j};
    elseif strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif strncmpi(vg,'cluster',2)
        j = j+1;
        thecluster = varargin{j};
    elseif strncmpi(vg,'play',4)
        playspikes = 1;
    elseif strncmpi(vg,'relistonly',8)
        mkidx = 2;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    end
    j = j+1;
end

Header.Name = BuildName(name);

Events.times = Events.times * 10000;
%spk times need to be ints for trigsdf.
Text.times = Text.times * 10000;

bsidx = zeros(size(Events.times));
if strncmp(Text.comment,'GridData',8)
    Expt.DataType = Text.comment;
else
    Expt.DataType = 'Spike2';
end





nt = 0;
nx = 0;
Expts = [];

%id = find(ismember(Spks.codes(:,1),thecluster));
Expt.setprobe = setprobe; 
frametimes = [];
bstimes = [];
estimes = [];

if exist(framechname,'var')
    if isfield(framech,'level')
        id = find(framech.level == 1);
        frametimes = framech.times(id) .* 10000
    else
    frametimes = eval([framechname '.times * 10000']);
    end
    Header.frameperiod = median(diff(frametimes));
else
    Header.frameperiod = 167;
end
if ~isempty(ustimmarkch)
    ustimmarkch.times = ustimmarkch.times .* 10000;
else
    ustimes = [];
end
if exist(stimch,'var') & isfield(stimlvl,'level') & length(stimlvl.times) > 1
    maxt = max([max(Text.times) max(Events.times) max(stimlvl.times)]);
%    id = find(Ch8.level == 1);
% need to add an extra event at the end of each list to avoid issues with
% empty finds on the last trial. But make it long after to avoid any
% confusion with the real one;
    if isfield(stimlvl,'inverted') && stimlvl.inverted
        bstimes = stimlvl.times(stimlvl.level == 0) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 1) * 10000;
        estimes(end+1) = maxt+50010;
    else
        bstimes = stimlvl.times(stimlvl.level == 1) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 0) * 10000;
        estimes(end+1) = maxt+50010;
    end
end

if length(bstimes) > length(estimes)
    bstimes = bstimes(1:length(estimes));
end
bsid = strmatch('bss',Text.text);
bsstimes = Text.times(bsid);
esid = strmatch('ess',Text.text);
esstimes = Text.times(esid);

Trials.bsstimes = bsstimes;
Trials.esstimes = zeros(size(bsstimes));
for j = 1:length(esstimes)
    id = find(bsstimes < esstimes(j));
    Trials.esstimes(id(end)) = esstimes(j);
end
iid = strmatch('id',Text.text);
Trials.idtimes = Text.times(iid);
showidstarts = 1;
 if showidstarts
     id = union(cat(1,bsid,esid),iid);
     for  j = 1:length(id)
         fprintf('%.3f: %s\n',Text.times(id(j))./10000,Text.text(id(j),:));
     end
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
storing = bitand(storing,STOREBIT);
%storeonoff is a list of times where storing is toggled
storeonoff = 1+find(abs(diff(storing)) > 0); 
if isempty(storeonoff)
    storeonoff = [1 1];
end

ve = strmatch('BGCS Version',Text.text);
if ve
    version = sscanf(Text.text(ve(1),:),'BGCS Version %f');
else
    version = 1.1;
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
nextonoff = 1;
if onlinedata
    storestate = 1;
else
    storestate = storing(1);
    tonoff = Text.times(opid(storeonoff(nextonoff)));
end

if onlinedata
Events.store = ones(size(Events.times));
else
tic;
Events.store = zeros(size(Events.times));
if ~storing(storeonoff(1)) & storing(1) %% first event is an off
    onid = find(Events.times < Text.times(opid(storeonoff(1))));
    Events.store(onid) = 1;
end

for j = 1:length(storeonoff)
    if storing(storeonoff(j))
        if length(storeonoff) > j
        onid = find(Events.times >= Text.times(opid(storeonoff(j))) ...
            & Events.times < Text.times(opid(storeonoff(j+1))));
        else
        onid = find(Events.times >= Text.times(opid(storeonoff(j))));
        end
        Events.store(onid) = 1;
    end
end
fprintf('Store Index takes %.2f\n',toc);
end

if onlinedata && sum(Events.codes(:,1)==STARTEXPT) == 0
    fprintf('Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    if logfid > 0
        fprintf(logfid,'Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    end
    id = find(Text.codes(:,1) == 1 & Text.times > 0.1);
    if ~isempty(id)
        Events.times = [Text.times(id(1)); Events.times];
        Events.store = [1; Events.store];
        Events.codes = [[STARTEXPT 0 0 0 ]; Events.codes];
        fprintf('True Start at %.1f\n',Text.times(id(1))./10000);
    end
end
nonstore = 0;
lastend = 0;
ntrials = sum(Events.codes(:,1) == FRAMESIGNAL);
trynew = 0;


%New bit
if length(estimes) == length(bstimes)+1
    estimes = estimes(2:end);
end

isi = bstimes(2:end)-estimes(1:end-1);

id = find(isi < 400); %shouldnt happen, but see lemM209.5
if length(id)
    Result.badisi = id;
    fprintf('%d impossible isis\n',length(id));
    useid = setdiff(1:length(bstimes), id+1);
    bstimes = bstimes(useid);
    estimes = estimes(useid);
    isi = bstimes(2:end)-estimes(1:end-1);
else
    Result.badisi = [];
end
exendid = find(Events.codes(:,1) == ENDEXPT | Events.codes(:,1) == CANCELEXPT);
exstartid = find(Events.codes(:,1) == STARTEXPT);
evid = [];
bsid = [];
ExptStart = Events.times(exstartid);
ExptEnd = Events.times(exendid);
ExptRes = Events.codes(exendid,1);
for j = 1:length(exstartid)
    id = find(ExptEnd > ExptStart(j));
    if isempty(id)
        Trials.Expts(1).End = Events.times(end);
        Trials.Expts(j).EndCode = 0;
    else
        Trials.Expts(j).End = ExptEnd(id(1));
        Trials.Expts(j).EndCode = ExptRes(id(1));
    end
    Trials.Expts(j).Start = ExptStart(j);
    ts = Trials.Expts(j).Start;
    te = Trials.Expts(j).End;
    
    id = find(Events.times > ts & Events.times < te);
    evid = cat(1,evid,id);
    fsid = find(Events.codes(id,1) ==  5);
    bid = find(bstimes > ts & bstimes < te);
    bsid = cat(1, bsid, bid);
    Result.startcounts = [length(fsid) length(bid)];
    Result.name = name;
    if isempty(fsid)
        fprintf('No Trials in Expt %d\n',j);
    elseif length(fsid) == length(bid)
        fprintf('Ex %d Stimlevel and FrameSignal Lengths match %d)\n',j,length(fsid));
        bsidx(id(fsid)) = bid;
        Result.bsdelay = Events.times(id(fsid)) - bstimes(bid);
    else
        fprintf('Length Mismatch for Stimlevel (%d) and FrameSignal (%d)\n',length(bid),length(fsid));
        Result.bsdelay = [];
    end
    Expts(j).start = Events.times(j);
    if length(bid)
    Expts(j).firsttrial = bid(1);
    Expts(j).lasttrial = bid(end);
    else
    Expts(j).firsttrial = NaN;
    Expts(j).lasttrial = NaN;
    end
    Expts(j).result = Trials.Expts(j).EndCode;
end

evid = unique(evid);
bsid = unique(bsid);
bfid = find(Text.codes(:,1) ==  11);
fsid = find(Events.codes(evid,1) ==  5);
fstimes = Events.times(evid(fsid));
bftimes = Text.times(bfid);
for j = 1:length(estimes)
    bt = bstimes(j);
    et = estimes(j);
    id = find(bftimes > bt & bftimes < et);
    if isempty(id)
        TRes(j) = 1;
    else
        TRes(j) = 0;
    end
    durs(j) = et-bt;
    if j > 1
        id = find(fstimes < et & fstimes  > estimes(j-1));
        if length(id) > 1
%            fprintf('%d FrameSignals at %.2f\n',length(id),et);
        elseif isempty(id)
  %          fprintf('No FrameSignals at %.2f\n',et);
        end
            
    end
end



if trynew %didn't help
Trials.Start(1:ntrials) = 0;
Trials.End(1:ntrials) = 0;
Trials.Trial(1:ntrials) = 0;
Trials.TrueEnd(1:ntrials) = 0;
Trials.Startev(1:ntrials) = 0;
Trials.stored(1:ntrials) = 0;
Trials.Result(1:ntrials) = 0;
Trials.serdelay(1:ntrials) = 0;
Trials.bstimes(1:ntrials) = 0;
Trials.delay(1:ntrials) = 0;
Trials.endelay(1:ntrials) = 0;
Trials.estimes(1:ntrials) = 0;
Trials.id(1:ntrials) = 0;
Trials.FalseStart(1:ntrials) = 0;
end

Trials.badbad = [];
Trials.latebad = [];
settrials = 0;
if length(fsid) == length(bsid)
   Result.bsdelay = fstimes - bstimes(bsid);

  for j = 1:length(fsid)
      Trials.Start(j) = bstimes(bsid(j));
      Trials.End(j) = estimes(bsid(j));
      Trials.bstimes(j) = Trials.Start(j);
      Trials.estimes(j) = Trials.End(j);
      Trials.Trial(j) = j+starttrial;
      Trials.Result(j) = 1;
%remove trials where bss is followed by another bss, not ess
%badfix can happen just before stimon if its after trial start
      id = find(esstimes > fstimes(j));
      if length(id) == 0 %happens if last trial is an early badfix
          Trials.Result(j) = -1;
      elseif j <length(fsid)
          if esstimes(id(1)) > fstimes(j+1)
              Trials.Result(j) = -1;
          end
      end
      if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(j) & frametimes < Trials.Start(j)+500);
                if ~isempty(id)
                    Trials.delay(j) = frametimes(id(1)) - Trials.Start(j);
                    Trials.Start(j) = frametimes(id(1));
                    Trials.FalseStart(j) = 0;
                else
                    Trials.FalseStart(j) = 1;
                    Trials.delay(j) = NaN;
                end

          
          id = find(frametimes > Trials.End(j));
          if ~isempty(id) & frametimes(id(1))-Trials.End(j) < 500
              Trials.endelay(j) = Trials.End(j) - frametimes(id(1));
              Trials.End(j) = frametimes(id(1));
              Trials.TrueEnd(j) = frametimes(id(1));
          else
              Trials.TrueEnd(j) = NaN;
          end
      end
  end
  %badfix is always (?? what about Sa?) detected by spike2, so the stimoff
  %marker will always be after this.
  bid = find(Events.codes(:,1) == 11); %bad fix
  sid = find(Events.codes(:,1) == 12) %start trial
  eid = find(Events.codes(:,1) == 16) %End trial
  for j = 1:length(bid)
      id = find(Trials.bstimes < Events.times(bid(j)));
      eid = find(Trials.bstimes < Events.times(bid(j)));
      if length(id)
          Trials.Result(id(end)) = 0;
          Trials.bfdelay(id(end)) = Trials.estimes(id(end))-Events.times(bid(j));
      else
          Trials.badbad = [Trials.badbad Events.times(bid(j))];
      end
      id = find(Events.times(sid) < Events.times(bid(j)));
      ied = find(Events.times(eid) < Events.times(bid(j)));
      if length(id) & length(ied) && Events.times(eid(ied(end))) > Events.times(sid(id(end))) %badfix after end
          Trials.latebad = [Trials.latebad Events.times(bid(j))];
      end
  end
  
  bid = find(Text.codes(:,1) == 11 & Text.codes(:,4) == 2); %bad fix from spike2
  gbid = strmatch('BAD Saccade',Text.text(bid,:));  %Did complete trial, but invalid psych sacc
  bid = setdiff(bid, bid(gbid));
  for j = 1:length(bid)
      t = Text.times(bid(j));
      id = find(Trials.bstimes < Text.times(bid(j)));
      if length(id)
      Trials.Result(id(end)) = 0;
      Trials.bfdelay(id(end)) = Trials.estimes(id(end))-Text.times(bid(j));
      else
          Trials.badbad = [Trials.badbad Text.times(bid(j))];
      end
      ied = find(Trials.estimes < Text.times(bid(j)));
      if length(id) & length(ied) && ied(end) >= id(end) ...
              && t > Trials.estimes(ied(end)) + 1000 %badfix after end
          Trials.latebad = [Trials.latebad Text.times(bid(j))];
      end
  end
      


  settrials = 1;
  nt = length(fsid);
end

fsid = find(Events.codes == FRAMESIGNAL);
lastbsid = 0;

if settrials == 0
for j = 1:size(Events.codes,1)

% if storage turned off mid-stim, don't want to miss ENSTIM marker
    if Events.codes(j,1) ~= ENDSTIM
        storestate = Events.store(j);
    end        
    if Events.codes(j,1) == FRAMESIGNAL
        nowt = Events.times(j);
        id = find(fsid > j);        
        if length(id)
            nextfsignal = Events.times(fsid(id(1)));
            ntev = id(1)-1;
        else
            nextfsignal = nowt + 300;
        end
        if storestate
            if nt & Trials.Result(nt) < 0
                nt = nt;
            end
            nt = nt+1;
        Trials.Start(nt) = Events.times(j);
        Trials.End(nt) = Trials.Start(nt)+nomdur; %% just in case an online file is missing end
        if findtrial  & Trials.Start(nt) > findtrial
            findtrial = 0;
        end
        Trials.Startev(nt) = Events.times(j);
        Trials.stored(nt) = storestate;
        Trials.Result(nt) = 1;
        if nt > 1 & length(Trials.End) < nt-1
            nerr = nerr+1;
            errs{nerr} = sprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
            fprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
%            Trials.End(nt-1) = NaN;
        end
        Trials.Trial(nt) = nt+starttrial;
        instim = 1;
%the event time can be just before the stim signal since it is not delayed
%to the vertical retrace (? correct explanation - it should be delayed...)
%
% the digital step can be >20ms after the serial signal is received, e.g.
% in ruf1989 at 209.43 sec
% in lem017 at 120.6 it is nearly 80ms late. at 289.6 its 300ms late. Looks
% also lemM209.3 at 179.5
%how is this possible?  Spike2 must be using an old timestamp to but with
%the event marker. 
% like we could check for immediately preceding STARTSTIM (6) being before 
% StimON to check for this

        id = find(bstimes < Events.times(j)+300 & bstimes < nextfsignal);
 %must mean that this Event belongs to a stimon > 30ms in the future
 %deal with problem above by checking to see if bstimes(id(end) was used
 %already. If this is so, must need one > 30ms in the future. 
        if id(end) == lastbsid 
            id = id(end)+1;
        end
        if Events.times(j) > 238170000  && 0 %was this a test? or is this the max time?
            bstimes(id(end));
        end
        if id
            Trials.serdelay(nt) = Events.times(j) - bstimes(id(end));
            Trials.bstimes(nt) = bstimes(id(end));
        end
        if (~isempty(id) & bstimes(id(end)) > lastend) || bsidx(j) > 0
            if bsidx(j) > 0
                Trials.Start(nt) = bstimes(bsidx(j));
            elseif id(end) < length(bstimes) & bstimes(id(end)+1) - Trials.Start(nt) < 10 %< 1ms to next = probably early
                Trials.Start(nt) = bstimes(id(end)+1);
            elseif bstimes(id(end)) > lastend
                Trials.Start(nt) = bstimes(id(end));
            end
            if Trials.serdelay(nt) > 10000
                nt = nt;
            end
            Trials.stored(nt) = storestate;
            if id(end) -lastbsid > 1 
                if Events.times(j) < bstimes(id(end))
                    id = id(1:end-1);
                else
                    fprintf('Skipped Digital StimOn %d - %d at %.1f\n',lastbsid,id(end),bstimes(id(end)));
                end
            end
            a = find(fsid == j);
            if bsidx(j) > 0
                bsid(a) = bsidx(j);
                lastbsid = bsidx(j);
            else
            bsid(a) = id(end);
            lastbsid = id(end);
            end
%here, Trials.Start is the time of the Digital event marker
%If vertical retrace was recorded, set start time to next one of these
%To a first approximation, it is the find(framemtimes ......  calls (here
%and for ENDSTIM) that take all the time in this loop. 
            if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(nt) & frametimes < Trials.Start(nt)+500);
                if ~isempty(id)
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                    Trials.FalseStart(nt) = 0;
                else
                    Trials.FalseStart(nt) = 1;
                    Trials.delay(nt) = NaN;
                end
            else
%                Trials.Start(nt) = Events.times(j);
                Trials.FalseStart(nt) = 2;
                Trials.delay(nt) = NaN;
            end
%        elseif  % bstimes can be > lastend
            
        else
            Trials.Start(nt) = Events.times(j);
            Trials.delay(nt) = NaN;
            if isempty(id)
                Trials.FalseStart(nt) = 1;
            elseif id(end) < length(bstimes) & bstimes(id(end)+1) - Events.times(j) < 800 ...
%                    & Events.codes(j-1) == STARTSTIM ... %< 1ms to next = probably early
            nerr = nerr+1;
                errs{nerr} = sprintf('StimON at %.2f is %.1f ms late but STARTSTIM at %.2f',...
                bstimes(id(end)+1),(bstimes(id(end)+1)-Events.times(j))./10,Events.times(j-1));
            fprintf('%s\n',errs{nerr}); 
                Trials.FalseStart(nt) = 0;
                
            else
%Serial input can be very late if Spike2 got busy. Use the DIO stimon -
%this is the true start
                Trials.FalseStart(nt) = Events.times(j) - bstimes(id(end));
                Trials.Start(nt) = bstimes(id(end));
                 nerr = nerr+1; errs{nerr} = sprintf('Missing StimON at %.2f %.2f), but STARTSTIM at %.2f',Events.times(j),bstimes(id(end)),Events.times(j-1));
                fprintf('%s\n',errs{nerr});
                id = find(frametimes > Trials.Start(nt));
                if ~isempty(id) 
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                end
            end
            Trials.stored(nt) = storestate;
        end 
        else
            nonstore = nonstore+1;
        end %if storestate
  
    elseif Events.codes(j,1) == ENDSTIM & storestate & nt
    if abs(Events.times(j) - 62740381) < 200
            Trials.End(nt) = Events.times(j);
            Trials.endelay(nt) = NaN;
    end
        if estimes
            id = find(estimes < Events.times(j)+500);
            if(id)
                Trials.TrueEnd(nt) = estimes(id(end));
                Trials.End(nt) = estimes(id(end));
                Trials.endelay(nt) = NaN;
                Trials.estimes(nt) = estimes(id(end));
  %if this is out by 400ms, probably failed to find correct end mark
                if Trials.TrueEnd(nt) < Events.times(j) - 4000  
                    fprintf('End event %.3f but marker %.3f\n',...
                        Events.times(j)./10000,Trials.TrueEnd(nt)./10000);
                    if Trials.End(nt) < Trials.Start(nt) & length(estimes) > id(end) & ...
                        estimes(id(end)+1) - Events.times(j) < 10000
                        Trials.TrueEnd(nt) = estimes(id(end));
                        Trials.End(nt) = estimes(id(end));
                        Trials.endelay(nt) = NaN;
                        Trials.estimes(nt) = estimes(id(end));
                    end
                end

                if ~isempty(frametimes) % have VTR channel
                    id = find(frametimes > Trials.End(nt));
                    if ~isempty(id) & frametimes(id(1))-Trials.TrueEnd(nt) < 500
                        Trials.endelay(nt) = Trials.End(nt) - frametimes(id(1));
                        Trials.End(nt) = frametimes(id(1));
                        Trials.TrueEnd(nt) = frametimes(id(1));
                    else
                    end
                end
            else
                Trials.TrueEnd(nt) = 0;
            end
            
        end
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 1;
        if (Trials.End(nt) - Trials.Start(nt)) < 1000
            instim = 0;
        end
        instim = 0;
        if Trials.TrueEnd(nt)
            lastend = Trials.TrueEnd(nt);
        else
            lastend = Trials.End(nt);
        end
    elseif Events.codes(j,1) == ENDTRIAL & storestate
        if instim
%  can't figure this out here because the BADFIX is only recorded in text,
%  not SampleKey (becuase this is send from Spike2, not received by, and
%  setting codes for sample keys is such a pain. But maybe should make all
%  of these events with code2 set to indicate it is from Spike2?
% Seems like this happens when fixation is broken just BEFORE stimulus on,
% but Spike2 has not registered this yet e.g. ruf2000 at 8867.9
%            fprintf('End Trial without End stim: %d (%.2f)\n',nt-1,Events.times(j)/10000);
            Trials.End(nt) = Events.times(j);
            Trials.Result(nt) = -1;  % this will be set to 0 if a BadFix is found.
            Trials.TrueEnd(nt) = NaN;
        end
    elseif Events.codes(j,1) == BADFIX & storestate %% Doesn't happen. Badfix is in Text, because it is sent, not received
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 0;
        instim = 0;
    elseif Events.codes(j,1) == STARTEXPT & storestate
        if inexpt %close an existing expt (e.g. if crashed out)
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
    elseif Events.codes(j,1) == ENDEXPT
        nx = nx;        
    end
end
end


Trials.id = zeros(size(Trials.Trial))'; %needs to be a row 
if nt == 0
    Expts = [];
    return;
end
if inexpt
    Expts(nx).lasttrial =nt;
end
ntrials = nt;
tic;

trial = 1;
k = 1;
ix = 1;
nx = 1;
% trynew should work, but need to check with some xxx= strings.
%problem with new method was cell2str deblanked. mat2cell works better
%
trynew = 1;
if trynew
    tic;
    
aText.text = mat2cell(Text.text(:,1:end-1),ones(1,size(Text.text,1)),size(Text.text,2)-1);
id = strmatch('xxx=',Text.text);
ids = setdiff(1:size(Text.text,1),id);
aText.text = aText.text(ids);
xid = ones(size(id));
nl = size(Text.text,1)-length(id);
k=1;
m = 1;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
        xid(m) = j;
        m = m+1;
% if we ever go back to this, probably need deblank(Text.....
    end
    if length(aText.text{k}) == 0 
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
    end
    k = k+1;
end

aText.times = Text.times(ids);
aText.codes = Text.codes(ids,:);
aText.text = deblank(aText.text);
end



ix = 1;
if trynew == 0
tic;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
    else
% if we ever go back to this, probably need deblank(Text.....
    aText.text{k} = Text.text(j,1:end-1);
    aText.times(k) = Text.times(j);
    aText.codes(k,:) = Text.codes(j,:);
    end
    if j == 206991
        k
    end
    if length(aText.text{k}) %%? safe may remove blank strings that have codes.
        k = k+1;
    else
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
        k = k+1;
    end
end
aText.text = deblank(aText.text);
fprintf('Reading Text %.2f\n',toc);
end
%
%AddTxtFile allows problems with data files to be fixed by adding
%additional lines written by hand. Fornmat is
%time  text
%where time is an int in timestamp units (0.1ms)
AddTxtFile = strrep(name,'.mat','Add.txt');
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%d %s','delimiter','\n');
    fclose(fid);
    t = a{1};
    s = a{2};
    for j = 1:length(t)
        %?why do we look for whitespace?? removed Aug 2010.
        id = findstr(s{j},' ');
        id = [];
        if length(id)
            txt = s{j}(id(1)+1:end);
        else
            txt = s{j};
        end
        if t(j) < 0 %special case for fixing lines
            if strncmp(s{j},'cm=rf',5)
                id = strmatch('cm=rf',aText.text);
                for k = 1:length(id)
                    aText.text(id,:) = s(j);
                end
            end
        elseif t(j) == 0 || t(j) <= aText.times(1)
            aText.text = {txt aText.text{:}};
            aText.times = [0; aText.times];
            aText.codes = [0 0 0 0; aText.codes];
        else
            id = find(aText.times < t(j));
            id = id(end);
            aText.text = {aText.text{1:id} txt aText.text{id+1:end}};
            aText.times = [aText.times(1:id); t(j); aText.times(id+1:end)];
            aText.codes = [aText.codes(1:id,:); 0 0 0 0; aText.codes(id+1:end,:)];
        end
    end
end
Trials.Stimseq = {};
intrial = 0;
Peninfo.trode = '';
nrw = 0;
tic;
tstart = now;
lasttook = 0;
Stimulus.CorLoop = 0;
Stimulus.SpikeGain = 50; %default
Stimulus.id = 0;
gotend = 0;
txtid = [];
Stimulus.Flag = '';

for j = 1:length(aText.text)
    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);

%    id = find(txt == 0);
%    txt = txt(1:id(1)-1);
    if ~isempty(txt);
        ss = txt(1:2);
        if length(txt) > 2 & txt(3) == '='
            val= txt(4:end);
        else
            val= txt(3:end);
        end
        if aText.codes(j,4) == 2  % this was FROM spike 3
            if aText.codes(j,1) == 3 && instim == 1 %end stim
                instim = 2;
            end

        elseif aText.codes(j,1) == 5 %stim start
            instim = 1;
            Stimulus.Seedseq = {};  %% these must be set for each stim
            Stimulus.Stimseq = {};
            Stimulus.Phaseseq = [];
            Stimulus.cLseq = [];
            Stimulus.cRseq = [];
            gotend = 0;
        elseif aText.codes(j,1) == 3 %end stim
            instim = 2;
        elseif strncmp(txt,'{}',2) %bug!!
        txt = aText.text{j};
            
        elseif strncmp(txt,'EndStim',7) %finished reading all text related to last stim
            gotend = 1;
        elseif strncmp(txt,'exvals',6)
        elseif strncmp(txt,'mixac',5)
            Stimulus.mixac = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'Off at',6) %Storage turned off - should be outside trial
            if intrial
                fprintf('Storage Off in Trial at %.1f',atext.times(j));
            end
            
        elseif strncmp(txt,'RightHemi',9) || strncmp(txt,'Electrode',8)
            txtid = [txtid j];
            Peninfo.trode = txt;
            a = InterpretLine(txt);
            Peninfo = CopyFields(Peninfo, a);
            id = strfind(txt,'Contact');
            if length(id)
                x = id(1);
                id = strfind(txt(id:end),' ');
                sscanf(txt(id+x:end),'%d',x);
                Peninfo.probesep = x;
            end

        elseif strncmp(txt,'cm=rf',5)
            a = sscanf(txt,'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
            Stimulus.rf = a;
        elseif strncmp(txt,'StartDepth',10)
            Stimulus.StartDepth = str2num(txt(11:end));
        elseif strncmp(txt,'CLOOP',5)
            Stimulus.CorLoop = 1;
%        elseif strncmp(txt,'id',2)
%           Stimulus.id = sscanf(txt(3:end),'%d')
        elseif strncmp(txt,'bt',2)
             g = sscanf(txt,'bt%d spkgain %f');
             if length(g) > 1 & g(2) > 1
                 Stimulus.SpikeGain = g(2);
             end
             if length(g) > 0
                 ExptStartTime = g(1);
             end
        elseif strncmp(txt,'rw',2)
            nrw = nrw+1;
            [Trials.rws(nrw), ok] = sscanf(val,'%f');
            Trials.rwset(nrw) = t;
        elseif strncmp(txt,'st',2)
            Stimulus.st = strmatch(val, stimnames,'exact');
            Stimulus.st = Stimulus.st -1;
        elseif strncmp(txt,'mtop=op',7)
            Stimulus.OptionCode = txt(8:end);
        elseif strncmp(txt,'fl+',3)
            Stimulus.Flag = txt(3:end);
        elseif strncmp(txt,'mtrP=',5)
            Stimulus.Phaseseq = sscanf(txt(6:end),'%d');
            if trial > length(Trials.Start)
                Trials.Phaseseq{trial} = Stimulus.Phaseseq;
            elseif instim ~= 1  && trial > 1 && t < Trials.Start(trial)
                Trials.Phaseseq{trial-1} = Stimulus.Phaseseq;
            end
            if Stimulus.id == 6136 || Stimulus.id > 570
                trial;
            end
        elseif strncmp(txt,'mtco=',5)
            Stimulus.Stimseq = sscanf(txt(6:end),'%d');
        elseif strncmp(txt,'mtcL=',5)
            Stimulus.cLseq = sscanf(txt(6:end),'%x');
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.cLseq{trial-1} = Stimulus.cLseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.cLseq{trial} = Stimulus.cLseq;
                end
            end
        elseif strncmp(txt,'mtcR=',5)
            Stimulus.cRseq = sscanf(txt(6:end),'%x');
            if sum(Stimulus.cRseq < 0) > 1
                Stimulus.cRseq = sscanf(txt(6:end),'%x');
            end
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.cRseq{trial-1} = Stimulus.cRseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.cRseq{trial} = Stimulus.cRseq;
                end
            end
        elseif strncmp(txt,'mtrS=',5)
            if aText.codes(j,1) == ENDSTIM
                istim = 2;
            end
            Stimulus.Stimseq = sscanf(txt(6:end),'%d');
            if Stimulus.id >= 5430
                Stimulus.Stimseq;
            end
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.Stimseq{trial-1} = Stimulus.Stimseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.Stimseq{trial} = Stimulus.Stimseq;
                end
            else
                instim;
            end
        elseif strncmp(txt,'mtse=',5)
            Stimulus.Seedseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.Seedseq{trial-1} = Stimulus.Seedseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.Seedseq{trial} = Stimulus.Seedseq;
                end
            else
                instim;
            end
            
        elseif strncmp(txt,'Nf',2) %comes after end TRIAL, not every stim
            if ~instim & trial > 1
                Trials.Nf(trial-1) = str2num(val);
            end
        elseif strncmp(txt,'mtet=',5)
            id = strfind(txt,'Fr');
            Stimulus.Fr = sscanf(txt(id+3:end),'%d');
        elseif strncmp(txt,'mtei=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
            if ix
                Expts(ix).e1vals = sscanf(txt(6:end),'%f');            
            else
                nerr = nerr+1;
                errs{nerr} = sprintf('No Expt for mtei at trial %d',trial);
                fprintf('%s\n',errs{nerr});
                ix = lastix;
            end
            end
        elseif strncmp(txt,'mte3',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
            lastix = ix;
            ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
            if ix
                Expts(ix).e3vals = sscanf(txt(6:end),'%f');
            else
                ix = lastix;
            end
            end
        elseif strncmp(txt,'mte2=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
            if ix
                Expts(ix).e2vals = sscanf(txt(6:end),'%f');
            else
                ix = lastix;
            end
            end
        elseif strncmp(txt,'Off at',5)
        elseif strncmp(txt,'EndExpt',5)
            ix = ix+1;
        elseif strncmp(txt,'sonull',5)
        elseif strncmp(txt,'NewConnect',7)
        elseif strncmp(txt,'BGCS Version',7)
        elseif strncmp(txt,'testflag',7)
        elseif strncmp(txt,'imve ',3)
            [a,b] = sscanf(txt,'imve %f,%f %f');
            Stimulus.imver = a(1);
            Stimulus.imseed = a(2);
            if length(a) > 2 & a(3) < 1
                Stimulus.impx = a(3);
            end
        elseif strncmp(txt,'Sa:',3)
            if aText.codes(j,4) == 1
            a =  aText.codes(j,1);
            end
        elseif strncmp(txt,'op',2)
            a = sscanf(txt(3:end),'%f,%f');
            Stimulus.op = a(1);
            if length(a) > 1
                Stimulus.optionb = a(2);
            end
            if isempty(Stimulus.op)
                fprintf('Missing op stim %d\n',trial);
                Stimulus.op = 0;
            end
        elseif strncmp(txt,'backMov',7)
            Stimulus.backMov = sscanf(txt(8:end),'%f');
        elseif strncmp(txt,'annTyp',6)
            Stimulus.annTyp = sscanf(txt(7:end),'%f');
        elseif ~isstrprop(ss(1),'alphanum')
            fprintf('Non-Printing Name %s\n',txt);
        elseif strncmp(txt,'fp',2)
            a = sscanf(val,'%f');
            if instim == 1
                Stimulus.dfx = a(1);
                Stimulus.dfy = a(2);
            else
                Stimulus.fx = a(1);
                Stimulus.dfx = a(1);
                Stimulus.fy = a(2);
                Stimulus.dfy = a(2);
            end
        elseif strncmp(txt,'fx',2) || strncmp(txt,'fy',2)
%if really in a stimulus, make note of new fx but keep original also
            a = sscanf(val,'%f');
            if instim == 1
                Stimulus.(['d' ss]) = a;
            else
                Stimulus.(ss) = a;
                Stimulus.(['d' ss]) = a;
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
            if strncmp(ss,'et',2)
                Stimulus.(ss) = val;
                if(ix) Expts(ix).et = val; end
            end
            if strncmp(ss,'e2',2) & length(ix)==1
                Expts(ix).e2 = val;
            end
            if strncmp(ss,'e3',2) & length(ix) ==1
                Expts(ix).e3 = val;
            end
        end
    end
    if (aText.codes(j,4) == 2 || isempty(txt))  && exist('Stimulus','var') %empty text = code from spike2 -> binoc
        correctdir = 0;
        if isfield(Stimulus,'OptionCode') & strfind(Stimulus.OptionCode,'+2a')...
                & isfield(Stimulus,'vs')
            [a,b] = max(abs([Stimulus.vs(1) Stimulus.sq(1)]));
%
% historically negative respdir means +ve sacccade value
% for exactly oblique saccades, sign of vertical component does it. 
            if b == 1
                correctdir = -sign(Stimulus.vs(1));
            else
                correctdir = -sign(Stimulus.sq(1));
            end
            Stimulus.rwdir = correctdir;
        end
        if trial == 3057
            instim = instim;
        end
%Real ON/Off times are set from the events above. But need to know that
%text following WURTZOK applies to the next stimulus. So instim = 2 means
%text has been received ending trial, but not officially over yet. 
%trial gets incremented at end stim. So response applies to trial -1
%but check that the time is sensible. Can get one of these events when
%storage is off, resetting the last stored trials
       if trial > 1
           tdelay = aText.times(j) - Trials.End(trial-1);
       else
           tdelay  = 0;
       end
        if aText.codes(j,1) == WURTZOKW & trial > 1 & tdelay < 10000
            Trials.RespDir(trial-1) = -1 * correctdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
  %          Trials.score(trial-1) = 0;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == WURTZOK & trial > 1  & tdelay < 10000
            Trials.RespDir(trial-1) = 1 * correctdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
 %           Trials.score(trial-1) = 1;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == BADFIX
            if trial <= length(Trials.End)
%With new method of finding trial end/start this can reset a good trial
%when fixation breaks just after trial over, so that a BADFIX is followed
%by WURTZOK in the sequence. lemM195.29 at 281.6
                if aText.times(j) > Trials.End(trial) && settrials == 0
                    Trials.End(trial) = aText.times(j);
                end
                if t < Trials.Start(trial)
%                if trial > 1 && t < Trials.End(trial-1)
%                    Trials.Result(trial-1) = 0;
%                end
            else
                Trials.Result(trial) = 0;
            end
            end
            instim = 2;
       end
    end
    if length(Trials.esstimes) >= trial && t >= Trials.esstimes(trial) & instim & trial <= length(Trials.Start)
%need to read past the end of the last trial a litle way to get things like
%Stimseq which come afterwards
        if trial > length(Trials.Start)
            break;
        end
        took = (now-tstart) * 24 * 60 *60;
        if took - lasttook > 30
            fprintf('%.0fsec..',took);
            lasttook = took;
        end
        if isfield(Stimulus,'st')
        Trials = SetTrial(Stimulus, Trials, trial, ntrials);
        Stimulus.CorLoop = 0;
        Stimulus.uf = '';
        if ~isempty(ustimmarkch)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1 ...
                & ustimmarkch.codes(:,1) == 1);
            if length(marks)
            Trials.uStimt{trial} = ustimmarkch.times(marks);
            elseif isfield(Trials,'optionb') && bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
            if bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
        end
%        AllStimuli(trial) = Stimulus; %fails when extra element is added to Stimulus
        trial = trial+1;
        end
        instim = 0;
    elseif trial > length(Trials.Start)
        if gotend
          break;
        end
    elseif t > Trials.Start(trial) && instim == 0
        instim = 1;
    end
end


cmid = strmatch('cm=',aText.text);
rfid = strmatch('cm=rf',aText.text);
bkid = strmatch('cm=noback',aText.text);
cmid = setdiff(cmid,[rfid; bkid]);
cmid = union(cmid,txtid);
Expt.Comments.text = {aText.text{cmid}};
Expt.Comments.times = aText.times(cmid);
Expt.Comments.Peninfo = Peninfo;

if length(Trials.op) < length(Trials.Start) 
        Trials = SetTrial(Stimulus, Trials, length(Trials.Start),ntrials);
end 
if isfield(Trials,'RespDir') &  length(Trials.RespDir) < length(Trials.Start)%fill in final trial
        Trials = SetTrial(Stimulus, Trials, length(Trials.Start),ntrials);
end
idx = find(Trials.Result < 0);
if ~isempty(idx)
    fprintf('%d Trials missing End\n',length(idx));
end
tic;
fn = fieldnames(Trials);
ntrials = length(Trials.Start);
cellids = {};
for j = 1:length(fn)
    if iscell(Trials.(fn{j}))
        cellids = {cellids{:} fn{j}};
        if length(Trials.(fn{j})) < ntrials
            Trials.(fn{j}){ntrials} = '';
        end
    end
end
for j = 1:length(Trials.Start)
    for k = 1:length(cellids)
        if isempty(Trials.(cellids{k}){j})
            Trials.(cellids{k}){j} = '';
        end
    end
end
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

if exist(dstimch,'var') & exist('stimchange','var')
    for j = 1:length(stimchange.times)
        et = stimchange.times(j) * 10000;
        id = find(Trials.Start < et);
        if length(id)
            t = id(end);
            if Trials.End(t) > et
                id = find(frametimes < et);
                if ~isempty(id)
                    diffs(j) = et-frametimes(id(end));
                    et = frametimes(id(end));
                end
                Trials.Events{t}{1,2} = et - Trials.Start(t);
                Trials.Events{t}{1,1} = 'ns';
            end
        end
    end
else
    Trials.Events{trial} = [];
end
tic;

Result.Header = Header;
GetFigure('TrialStimOn');
hold off;
for j = 1:length(bstimes)
    plot([bstimes(j) bstimes(j) estimes(j) estimes(j)],[0.9 1 1 0.9],'k');
    hold on;
end
nt = min([length(bsid) length(fstimes)]);
for j = 1:nt
    if bsid(j) > 0
        plot([fstimes(j) bstimes(bsid(j))],[1.1 1]);
    end
end
plot(fstimes,1.1,'rx');
for j = 1:length(exendid)
    t = Events.times(exendid(j));
    if t > Events.times(exstartid(1))
        plot([t t],[0.9 1],'r-');
    end
end

durs = Trials.End - Trials.Start;
id = find(Trials.Result == 0);
if length(id)
plot(Trials.End(id),1.02,'bx');
end
id = find(Trials.Result < 0);
if length(id)
plot(Trials.End(id),1.02,'kx');
end
id = find(Events.codes(:,1) == 16);
if length(id)
plot(Events.times(id),1.17,'kx');
end


tdur = median(durs(Trials.Result == 1));
id = find(durs < 0.9 * tdur & Trials.Result == 1);
for j = 1:length(id)
    plot([Trials.End(id(j)) Trials.End(id(j))],[1 1.2],'m-');
end
id = find(abs(durs-tdur) > 150  & Trials.Result == 1);
for j = 1:length(id)
    plot([Trials.End(id(j)) Trials.End(id(j))],[1 1.2],'g-');
end

id = strmatch('bss',Text.text);
if length(id)
plot(Text.times(id),1.198,'ms','buttondownfcn',@HitTextPoint);
end
id = strmatch('ess',Text.text);
if length(id)
plot(Text.times(id),1.202,'gs','buttondownfcn',@HitTextPoint);
end
id = find(Text.times > Events.times(exstartid(1)));
if length(id)
plot(Text.times(id),1.2,'bx','buttondownfcn',@HitTextPoint);
end
id = find(Text.codes(:,1) == 11); %bad fix
if length(id)
plot(Text.times(id),1.17,'bs','buttondownfcn',@HitTextPoint);
end

id = find(Events.codes(:,1) == 11); %bad fix
if length(id)
plot(Events.times(id),1.18,'ks','buttondownfcn',@HitEventsPoint);
end
id = find(Events.codes(:,1) == 3); %end stim
if length(id)
plot(Events.times(id),1.182,'gs','buttondownfcn',@HitEventsPoint);
end
id = find(Events.codes(:,1) == 5); %frame signal
if length(id)
plot(Events.times(id),1.178,'ms','buttondownfcn',@HitEventsPoint);
end
id = find(Events.codes(:,1) == 1); %start expt
if length(id)
plot(Events.times(id),1.18,'kx','buttondownfcn',@HitEventsPoint);
end
id = find(Result.bsdelay > 300);
for j = 1:length(id)
    plot(fstimes(id(j)),1.05,'rs');
end
id = strmatch('uf',Text.text(1:10,:));
filename = strrep(Text.text(id,3:end),'\','/');
title(filename);
drawnow;
toc

fprintf('Events bs %d, fs %d\n',length(bstimes),length(fstimes));
Result.Trials = Trials;
DATA.Text = Text;
DATA.Events = Events;
DATA.Trials = Trials;
DATA.toplevel = gcf;
set(gcf,'UserData',DATA);
set(gcf,'keypressfcn',{@KeyPress});

function HitPopPoint(a,b, id)

DATA = GetDataFromFig(a);
if ~isempty(id)
    if isfield(DATA.Result{id},'name')
    fprintf('%d: %s\n',id,DATA.Result{id}.name);
    else
        fprintf('Res %d\n',id);
    end
    GetFigure('TestStim');
    hold off;
    sid = find(DATA.Result{id}.Trials.Result ==1);
    durs = DATA.Result{id}.Trials.End(sid) - DATA.Result{id}.Trials.Start(sid); 
    sdurs = DATA.Result{id}.Trials.estimes(sid) - DATA.Result{id}.Trials.bstimes(sid); 
    if DATA.plottype == 2
        plot(sdurs,durs,'o');
    elseif DATA.plottype == 1
        plot(DATA.Result{id}.Trials.nf(sid),durs,'o');
    else
    hist(durs);
    end
end

function HitTextPoint(a,b)

DATA = GetDataFromFig(a);
pos = get(gca,'currentpoint');
x = get(a,'Xdata');
id = find(DATA.Text.times == x);
if ~isempty(id)
    fprintf('%s\n',DATA.Text.text(id,:));
end
DATA.lasttxtid = id;
set(DATA.toplevel,'UserData',DATA);

function HitEventsPoint(a,b)

DATA = GetDataFromFig(a);
pos = get(gca,'currentpoint');
x = get(a,'Xdata');
id = find(DATA.Events.times == x);
if ~isempty(id)
    fprintf('%d at %.3f\n',DATA.Events.codes(id,1),x);
    if DATA.Events.codes(id,1) == 11
        id = find(DATA.Trials.estimes > x);
        fprintf('Next Off at %.0f, Last %.0f\n',DATA.Trials.estimes(id(1)),DATA.Trials.estimes(id(1)-1));
    end
end

function KeyPress(src, ks)

isshift = length(strmatch('shift',ks.Modifier));
    xl = get(gca,'xlim');
    xw = diff(xl);
if strmatch(ks.Key,'rightarrow')
    if isshift
    else
    xl = get(gca,'xlim');
    xw = diff(xl);
    xl(1) = xl(2);
    xl(2)= xl(1) + xw;
    set(gca,'xlim',xl);
    end
elseif strmatch(ks.Key,'leftarrow')
    xl(2) = xl(1);
    xl(1)= xl(2) - xw;
    set(gca,'xlim',xl);
elseif strmatch(ks.Key,{'downarrow' 'uparrow'})
   DATA = GetDataFromFig(src);
   if DATA.lasttxtid <= 0
       id = find(DATA.Text.times > xl(1))
       DATA.lasttxtid = id(1);
   end
   if strmatch(ks.Key,'downarrow')
       id = DATA.lasttxtid+1:DATA.lasttxtid+10;
   else
       id = DATA.lasttxtid-1:-1:DATA.lasttxtid-10;
   end
   
   for j = id
       fprintf('%.1f: %s\n',DATA.Text.times(j)./10000,DATA.Text.text(j,:));
   end
   DATA.lasttxtid = j;
   set(DATA.toplevel,'UserData',DATA);
elseif ks.Character == '-'    %zoom out
    xm = mean(xl);
    xl(1) = xm-xw;
    xl(2) = xm+xw;
    set(gca,'xlim',xl);
elseif ks.Character == '+' %zoom in
    xm = mean(xl);
    xl(1) = xm-xw/4;
    xl(2) = xm+xw/4;
    set(gca,'xlim',xl);
end





function ExptList = MkExList(Expts)
ExptList = [];
if ~iscell(Expts)
    return;
end
for j =1:length(Expts)
        ExptList(j).expname = Expts{j}.Header.expname;
        ExptList(j).start = Expts{j}.Header.Start;
        ExptList(j).end = Expts{j}.Header.End;
        ExptList(j).et = Expts{j}.Stimvals.et;
        ExptList(j).e2 = Expts{j}.Stimvals.e2;
        ExptList(j).e3 = Expts{j}.Stimvals.e3;
end

function Expts = AddComments(Expts, Expt)

if isfield(Expt,'Comments')

for j = 1:length(Expt.Comments.times)
    id = find(Expt.Trials.Start < Expt.Comments.times(j));
    if isempty(id)
        Expt.Comments.id(j) = Expt.Trials.id(1);
    else
        Expt.Comments.id(j) = Expt.Trials.id(id(1));
    end
end


    for j = 1:length(Expts)
        cid = find(Expt.Comments.times > Expts{j}.Header.trange(1)-100000 & ...
       Expt.Comments.times < Expts{j}.Header.trange(2)+100000);
        if ~isempty(cid)
            bid = strmatch('cm=back=',Expt.Comments.text(cid));
            for k = 1:length(bid)
                %these lines are before start expt. If after end its for
                %next expt
                if Expt.Comments.times(cid(bid(k))) < Expts{j}.Header.trange(2)
                    Expts{j} = ParseComment(Expts{j}, Expt.Comments.text{cid(bid(k))});
                end
            end
        end
        Expts{j}.Comments.text = {Expt.Comments.text{cid}};
        Expts{j}.Comments.times = {Expt.Comments.times(cid)};
        cid = strmatch('cm=VisualArea',Expt.Comments.text);
        for k = 1:length(cid)
            Expts{j}.Header.Area{k} = Expt.Comments.text{cid(k)}(15:end);
        end
    end
end



function Expt = ParseComment(Expt,s)

if strmatch('cm=back',s)
    sscanf(s,'cm=back=%s,wi=%f,hi=%f');
    c = findstr(s,',');
    Expt.Stimvals.backstim = s(9:c(1)-1);
    b = sscanf(s(c(1):end),',%2s=%*f');
    a = sscanf(s(c(1):end),',%*2s=%f');
    for j = 1:length(a)
        f = ['back' b(j*2-1) b(j*2)];
        Expt.Stimvals.(f) = a(j);
    end
    if 0 %old method. Fails if the format changes...
    a = sscanf(s(c(1):end),',wi=%f,hi=%f,ce=%f,xo=%f,yo=%f,dx=%f,co=%f');
    Expt.Stimvals.backwi = a(1);
    Expt.Stimvals.backhi = a(2);
    Expt.Stimvals.backxo = a(4);
    Expt.Stimvals.backyo = a(5);
    Expt.Stimvals.backdx = a(6);
    Expt.Stimvals.backco = a(7);
    end
end

function Trials = SetTrial(Stimulus, Trials, trial, ntrials)

fn = fieldnames(Stimulus);
if isfield(Stimulus,'Ro') && isfield(Stimulus,'dx')
    ca = cos(Stimulus.Ro * pi/180);
    sa = sin(Stimulus.Ro * pi/180);
% need to recheck sign conventions here in replay...    
    Stimulus.Op = Stimulus.yo .* ca - Stimulus.xo .* sa;
    Stimulus.Pp = Stimulus.yo .* sa + Stimulus.xo .* ca;
    Stimulus.dO = Stimulus.dx .* sa - Stimulus.dy .* ca;
    Stimulus.dP = Stimulus.dy .* sa + Stimulus.dx .* ca;
end
for k = 1:length(fn)
    F = fn{k};
    if ~isfield(Trials,F)
%        Trials.F(1:ntrials) = NaN; % pre-allocate memory
    end
    if strcmp(F,'St')
        Trials.St(trial) = Stimulus.(F);
    elseif strncmp(F,'trode',5)
    elseif ischar(Stimulus.(F)) 
        if ~isfield(Trials,F) || ischar(Trials.(F)) || iscell(Trials.(F))
        if length(Stimulus.(F)) > 0
            Trials.(F){trial} = Stimulus.(F);
        end
        else
            fn{k}
        end
    elseif strcmp(F,'Seedseq') || strcmp(F,'Stimseq') || strcmp(F,'Phaseseq')   || strcmp(F,'cLseq')  || strcmp(F,'cRseq')    
           Trials.(fn{k}){trial} = Stimulus.(fn{k});
    else
        if isempty(Stimulus.(F))
            Trials.(F)(trial) = NaN;
        else
            Trials.(F)(trial,1:length(Stimulus.(F))) = Stimulus.(F);
        end
    end
end
if isempty(Stimulus.st) || isnan(Stimulus.st) || isnan(Trials.st(trial))
    Stimulus.st
end
if isfield(Trials,'rwset') & isfield(Trials,'rws')
id = find(Trials.rwset < Trials.Start(trial));
eid = find(Trials.rwset < Trials.End(trial));
if length(id)
Trials.rw(trial) = Trials.rws(id(end));
end
end
if isfield(Trials,'RespDir') & length(Trials.RespDir) < length(Trials.Start)
    Trials.RespDir(length(Trials.Start)) = 0;
end

function [Expts, Idx] = SortExpts(AllExpts, AllTrials, Header, thecluster, Idx, varargin)
SpkDefs;
timeoffset = 0;
%stimnames = {'None', 'Gabor', 'RDS', 'Grating', 'bar', 'circle', 'rectangle', 'test', 'square', 'probe', '2grating', 'Cylinder', 'twobar', 'rls', 'annulus', 'rdssine', 'nsines'};
Expts = [];
spikid = Idx.Spkid;
if ~isfield(Idx,'newerrs')
    Idx.newerrs = 0;
end
if isfield(Header,'frameperiod')
frameperiod = Header.frameperiod;
else
frameperiod = 167;
end
findtrial = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'findtrial',5)
        j = j+1;
        findtrial = varargin{j};
    elseif strncmpi(varargin{j},'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    j = j+1;
end

STIM_BAR = 4;
STIM_GRATING = 3;
for j = 1:length(stimnames)
    if strmatch(stimnames{j},'bar')
        STIM_BAR = j-1; %starts with 0
    end
end
fn = fieldnames(AllTrials);
% make a list of fileds that are NOT automatically set 
ids = strmatch('Spikes',fn);
ids = [ids strmatch('OptionCode',fn)];
ids = [ids strmatch('Start',fn,'exact')];
ids = [ids strmatch('End',fn)];
ids = [ids strmatch('Result',fn)];
ids = [ids strmatch('op',fn)']; % has op and optionb
ids = [ids strmatch('Stimseq',fn)];
ids = [ids strmatch('Seedseq',fn)];
ids = [ids strmatch('cLseq',fn)];
ids = [ids strmatch('cRseq',fn)];
ids = [ids strmatch('Phaseseq',fn)];
ids = [ids strmatch('Cluster',fn)];
ids = [ids strmatch('Frames',fn)];
ids = [ids strmatch('StartEv',fn)];
ids = [ids strmatch('Events',fn)];
ids = [ids strmatch('rws',fn,'exact')];
ids = [ids strmatch('endelay',fn)];
ids = [ids strmatch('rwset',fn)];
ids = [ids strmatch('estimes',fn)];
ids = [ids strmatch('uStimt',fn)];
ids = [ids strmatch('Flag',fn,'exact')];
%ids = [ids strmatch('imver',fn)];
%ids = [ids strmatch('imse',fn)];

fn = fn(setdiff([1:length(fn)],ids));
%non-zero values of falsestart should indicate the time gap
if isfield(AllTrials,'TrueStart') & ~isfield(AllTrials,'FalseStart')
    AllTrials.FalseStart = AllTrials.TrueStart;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 0)) = 1;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 1)) = 0;
end
if ~isfield(AllTrials,'Nf') || length(AllTrials.Nf) < length(AllTrials.Start)
    AllTrials.Nf(length(AllTrials.Start)) = NaN;
end

if ~isfield(AllTrials,'Events') || length(AllTrials.Events) < length(AllTrials.Start)
    AllTrials.Events{length(AllTrials.Start)} = [];
end
%If there are a minority of trials with False Starts (usually missing the
%StimON channel, or a long delay to the matching FRAMESIGNAL on the serial line
%set the result to 2, so that these can be exluded if necessary.
% only do this if the result is already not 0 (otherwise Add in Badfix
% trials
if sum(AllTrials.FalseStart > 0) < length(AllTrials.FalseStart)/5
    AllTrials.Result(find(AllTrials.FalseStart > 0 & AllTrials.Result > 0)) = 2;
end
    %fn is now a list of fields that are set for each trial automatically.
nexpts = 1;
    phasevals  = [0:360];
    phasevals(2:4) = [pi pi/4 3*pi/4];
    lphasevals = [0 pi pi 0]; %%see binoc.c SetRandomPhase();
    rphasevals = [0 pi 0 pi];
    for nx = 1:length(AllExpts)
    
    clear Trials;
    frpt = 1;
    a = AllExpts(nx).firsttrial:AllExpts(nx).lasttrial;
    igood = find(AllTrials.Result(a) > 0);
    nt = length(igood);
    igood = a(igood);
    nu = 0; %number of ustim pulses
    if nt> 3 & igood(1) < length(AllTrials.Start)
       spkids = [];
       needfields = {};
       if strmatch(AllExpts(nx).e3,{'ar'})
           needfields = {needfields{:} 'wi' 'hi'};
       end
       for j = 1:nt
            if findtrial & AllTrials.Start(igood(j)) > findtrial
                findtrial = 0;
            end
            [Trials(j).Start] = AllTrials.Start(igood(j));
            [Trials(j).TrialStart] = AllTrials.Start(igood(j));
            if isfield(AllTrials,'TrueEnd') & AllTrials.TrueEnd(igood(j)) > 0
                [Trials(j).End] = AllTrials.TrueEnd(igood(j));
            else
                [Trials(j).End] = AllTrials.End(igood(j));
            end
            Trials(j).dur = Trials(j).End(end)-Trials(j).Start(1);
            if isfield(AllTrials,'optionb')
            [Trials(j).uStim] = bitand(AllTrials.optionb(igood(j)),64);
            end
            [Trials(j).op] = AllTrials.op(igood(j));
            [Trials(j).Trial] = AllTrials.Trial(igood(j));
            [Trials(j).id] = AllTrials.id(igood(j));
            if ~isempty(AllTrials.Events{igood(j)})
                Trials(j).Events = AllTrials.Events{igood(j)};
            end
            if isfield(Trials,'uStimt')
                nu = nu + legnth(Trials(j).uStimt);
            end
%            spkids = [spkids spkid(igood(j),1):spkid(igood(j),2)];
        end
 % cant do this. for some reason max spkids) is > size(Spikes.times)
%        Expt.gui.spks = unique(spkids);
%  problem is that we don't rebuilds Spkid if set s probe != 1, so save
%  Spkid no good
 %      
        for nf = 1:length(fn)
            if iscell(AllTrials.(fn{nf}))
                if strcmp(fn{nf},'uStimt')
                else
                nv = unique({AllTrials.(fn{nf}){igood}});
                Expt.Stimvals.(fn{nf}) = AllTrials.(fn{nf}){a(1)};
                end
            else
                nv = unique([AllTrials.(fn{nf})(igood)]);
            end
            
            if isnumeric(nv) & sum(~isnan(nv))> 1
                for j = 1:nt
                    Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j));
                end
                if strcmp(fn{nf},'Nf')
                    if sum(nv) == max(nv) % only 1 value + 0;
                        Expt.Stimvals.(fn{nf}) = max(nv);
                    end
                    Trials = rmfield(Trials,'Nf');
                end
                Expt.Stimvals.(fn{nf}) = prctile([AllTrials.(fn{nf})(igood)],50);
            elseif isnumeric(nv)
                Expt.Stimvals.(fn{nf}) = nv;
            end
        end
        if isfield(AllTrials,'rf');
            Expt.Stimvals.rf = median(AllTrials.rf(igood,:));
            Expt.Stimvals.st = mode(AllTrials.st(igood));
        end
    duration = mean([Trials.End] - [Trials.Start]);
    et = Expt.Stimvals.et;
    e2 = Expt.Stimvals.e2;
    if isfield(Expt.Stimvals,'Fr') && Expt.Stimvals.Fr > 0
        frpt = Expt.Stimvals.Fr;
    end
    if isfield(AllExpts(nx),'e1vals') & ~isempty(AllExpts(nx).e1vals)
        Expt.e1vals = AllExpts(nx).e1vals;
    end
    if strmatch('ce',{et e2})
        Expt.Stimvals.ce = median(abs(AllTrials.ce(igood)));
    end
     if strmatch(et,{'Op' 'Pp'},'exact') & isfield(Expt.Stimvals,'rf')
            ca = cos(Expt.Stimvals.rf(5) * pi/180);
            sa = sin(Expt.Stimvals.rf(5) * pi/180);
            rOp = Expt.Stimvals.rf(2) .* ca - Expt.Stimvals.rf(1) .* sa;
            rPp = Expt.Stimvals.rf(2) .* sa + Expt.Stimvals.rf(1) .* ca;
            Expt.Stimvals.rOp = rOp;
            Expt.Stimvals.rPp = rPp;
    if isfield(AllExpts(nx),'e1vals')
% don't mess with the values for interleaved extras
            sid = find(AllExpts(nx).e1vals > -1000);
            if strcmp(et,'Op')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rOp;
            elseif strcmp(et,'Pp')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rPp;
            end
        end
    end
    if isfield(AllExpts(nx),'e2vals')
        Expt.e2vals = AllExpts(nx).e2vals;
    end
    allev = [];
    serrid = [];
    nframes = [];
    psychtrial = 0;
    seqtrial = 0;
    crtrial = 0; %count # with contrast reversal
    nu = 0; %count trials with uStimt;
    for k = 1:nt
        sid = find(ismember(AllTrials.Cluster{igood(k)}, thecluster));
        Trials(k).Spikes = round(AllTrials.Spikes{(igood(k))}(sid));
        Trials(k).count = sum(find(Trials(k).Spikes > 500 & Trials(k).Spikes < duration+500));
        Trials(k).sz = AllTrials.wi(igood(k));
        Trials(k).OptionCode = AllTrials.OptionCode{igood(k)};
        if isfield(AllTrials,'Flag')
        Trials(k).Flag = AllTrials.Flag{igood(k)};
        if strfind(Trials(k).Flag,'+mm')
            Trials(k).flatsurf = 1;
        else
            Trials(k).flatsurf = 0;
        end
        end
        if isfield(AllTrials,'uStimt')
        Trials(k).uStimt = round(AllTrials.uStimt{(igood(k))});
        nu = nu + length(Trials(k).uStimt);
        end
        if strfind(Trials(k).OptionCode,'+2a')
            psychtrial = psychtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+fS')
            seqtrial = seqtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+cr')
            crtrial = crtrial+1;
        end
        if ~isnan(AllTrials.Frames(igood(k)))
            Trials(k).Frames = AllTrials.Frames(igood(k));
        end
        if bitand(LMONOC,Trials(k).op)
            Trials(k).me  = -1;
        elseif bitand(RMONOC,Trials(k).op)
            Trials(k).me  = 1;
        else
            Trials(k).me  = 0;
        end
        if isfield(Trials,'Fr')
            frpt = Trials(k).Fr;
        end
% for image seqs, Stimseq records the order of seeds, and there will not
% be a conversion to stimulus type
        if isfield(AllTrials,'Seedseq') && length(AllTrials.Seedseq{igood(k)}) > 1
            Trials(k).Seedseq = AllTrials.Seedseq{igood(k)};
        end
        for f = 1:length(needfields)
            Trials(k).(needfields{f}) = AllTrials.(needfields{f})(igood(k));
        end
        if isfield(AllTrials,'Stimseq') && isfield(Expt,'e1vals') && ~isempty(Expt.e1vals)
            if length(Expt.e1vals) && ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ((max(AllTrials.Stimseq{igood(k)}) < length(Expt.e1vals) && seqtrial > k/2) || strcmp(et,'backMov'))
            evid = AllTrials.Stimseq{igood(k)}+1;
            if frpt > 1
                evid = evid(1:frpt:end);
            end
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
            if strmatch(et,'Dc')
                ev = Expt.e1vals(evid);
                if frpt > 1
                    ev = ev(1:frpt:end);
                end
                Trials(k).(et) = AllTrials.Dc(igood(k));
                Trials(k).(Expt.Stimvals.e2)= ev;
                Trials(k).(e2)(find(ev == ISIGNALFRAME)) = AllTrials.(e2)(igood(k)); %% blanks
                if strcmp(Expt.Stimvals.e2,'or')
                    Trials(k).ori = AllTrials.(e2)(igood(k));
                end
            elseif strmatch(et,'backMov')
                ev = AllTrials.Stimseq{igood(k)}+1;
                Trials(k).(et) = ev;
             %   eb = Expt.e2vals(AllTrials.Stimseq{igood(k)}+1);
              %  Trials(k).(Expt.Stimvals.e2)= eb;
            else
                ev = Expt.e1vals(evid);
                Trials(k).(et) = ev;
                 if length(Expt.e2vals) > max(AllTrials.Stimseq{igood(k)})
                    eb = Expt.e2vals(evid);
                    Trials(k).(Expt.Stimvals.e2)= eb;
                end
            end
            if isfield(AllTrials,'cLseq') & length(AllTrials.cLseq{igood(k)}) > 0
                Trials(k).cL = AllTrials.cLseq{igood(k)};
                Trials(k).cL(Trials(k).cL < 0) = 0;
                Trials(k).cL(Trials(k).cL > 512) = 0;
            end
            if isfield(AllTrials,'cRseq') & length(AllTrials.cRseq{igood(k)}) > 0
                Trials(k).cR = AllTrials.cRseq{igood(k)};
                Trials(k).cR(Trials(k).cR < 0) = 0;
                Trials(k).cR(Trials(k).cR > 512) = 0;
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 0
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            Trials(k).ce(find(ev == IUNCORR)) = 0; %% uncorr
            Trials(k).me(find(ev == ILEFTMONOC)) = -1; 
            Trials(k).me(find(ev == IRIGHTMONOC)) = 1;
                if size(Trials(k).(et),2) > 1
                    size(evid)
                end
                if isfield(Trials,e2) & size(Trials(k).(e2),2) > 1
                    size(evid)
                end
%in the Lopos X Ropos expt, could exceed max stim combinatinos, then mtei
%and mte2 strings only gave the inidividual Lpos/Rpos values, not all
%combinations
        elseif isfield(AllTrials,'Stimseq') && AllTrials.ve(end) < 4.85 && ...
                ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ~isempty(strfind(Trials(k).OptionCode,'+fS')) && ...
             max(AllTrials.Stimseq{igood(k)}) > length(Expt.e1vals)
         extras = sum(Expt.e1vals < -999);
         evid = AllTrials.Stimseq{igood(k)}-extras;
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
         n2 = length(Expt.e2vals)-extras;
         n1 = length(Expt.e1vals)-extras;
         e1 = mod(evid,n2)+extras+1;
         id = find(evid <0);
         e1(id) = evid(id)+extras+1;
         eb = floor(evid./length(Expt.e1vals))+1+extras;
         eb(id) = evid(id)+extras+1;
         Trials(k).(et) = Expt.e1vals(e1);
         Trials(k).(et)(id) = 0;
         ev = Expt.e1vals(e1);
         if ~isempty(Expt.e2vals)
         Trials(k).(e2) = Expt.e2vals(eb);
         end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            if Trials(k).id == 6138
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 0
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
        elseif isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 2 && seqtrial > k/2
            serrid = [serrid k];
            end
        end
        
        durs(k) = Trials(k).End(end) - Trials(k).Start(1);
        nframes(k) = length(Trials(k).Start);
    end
    if isfield(Trials,'CorLoop')
        id = find([Trials.CorLoop] == 0);
        if length(id)
            Trials = Trials(id);
            nframes = nframes(id);
            durs = durs(id);
        end
    end
    if isfield(Trials,'flatsurf')
        fs = unique([Trials.flatsurf]);
        if length(fs) == 1
            Trials = rmfield(Trials,'flatsurf');
            Expt.Stimvals.flatsurf = fs;
        end
    end
    if mean(nframes > 2)
       id = find(nframes > 1);
       Trials = Trials(id);
    end
    if nu == 0  && isfield(Trials,'uStimt')
        Trials = rmfield(Trials,'uStimt');
    end
    if isfield(Trials,'uStim') && sum([Trials.uStim]) == 0
            Trials = rmfield(Trials,'uStim');
    end


    id = [];
    if isfield(Trials,'delay')
        id = find(~isnan([Trials.delay]));
    end
    if length(id) > length(Trials)/2
        nid = find(isnan([Trials.delay]));
        Expt.Trials = Trials(id);
        durs = durs(id);
        if length(nid)
            if isfield(Trials,'FalseStart')
            fprintf('Delay Nan at %.2f(%.2f), id%.0f\n',Trials(nid(1)).Start(1),Trials(nid(1)).FalseStart, Trials(nid(1)).id);
            else
            fprintf('Delay Nan at %.2f, id%.0f\n',Trials(nid(1)).Start, Trials(nid(1)).id);
            end
        end
    else
        Expt.Trials = Trials;
    end
    Expt.Header = Header;
    Expt.Header.trange(1) = AllExpts(nx).start;
    if isfield(AllExpts,'end')
        AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
    elseif nx == 1  %%online, unfinished, edpt
        AllExpts(nx).end = Expt.Trials(end).End(end)+10000;        
    end
    Expt.Header.trange(2) = AllExpts(nx).end;
    if psychtrial > nt/2
        Expt.Header.psych = 1;
    else
        Expt.Header.psych = 0;
    end
    if seqtrial > nt/2
        Expt.Header.rc = 1;
    else
        Expt.Header.rc = 0;
    end
    Expt.Header.Options = '';
    if crtrial > nt/2
        Expt.Header.Options = [Expt.Header.Options '+cr'];
    end

    
    [a, expname, b, stimname] = Expt2Name(Expt);
    if isempty(expname)
        Expt.Header.expname = 'None';
    else
        Expt.Header.expname = [stimname '.' expname];
    end
    Expt.Header.Start = AllExpts(nx).start+timeoffset;
    Expt.Header.End = AllExpts(nx).end+timeoffset;
 
    if isfield(Expt.Trials,'cL')
 %       vals = [Expt.Trials.cL]; %unsafe syntax = varying Fr
        if median(nframes) > 1 %% subspace map for cL
            Expt.Header.expname = [Expt.Header.expname 'C'];
        end
    end
    if Expt.Header.rc
        Expt.Header.expname = [Expt.Header.expname 'RC'];
    end
    if min(durs) < 0
        id = find(durs < 0);
        err = sprintf('Expt %d at %.2f(id%d) has (%d) negative durations',nexpts,Expt.Trials(id(1)).Start(1),Expt.Trials(id(1)).id,length(id));
        Idx = AddError(err, Idx);
        id = find(durs > 0);
        Expt.Trials = Trials(id);
    end
    if ~isfield(Expt,'e1vals') & isfield(Expt.Trials,et)
        Expt.e1vals = unique([Expt.Trials.(et)]);
    end

    if strcmp(Expt.Stimvals.et,'fx') && strcmp(Expt.Stimvals.e2,'fy')
        emname = strrep(Expt.Header.Name,'.mat','.eyecal.mat');
        if ~exist(emname,'file')
            MakeEyeCalfile(Expt,emname);
        else
            load(emname);
            fprintf('Eye gains RH %.2f LH %.2f RV %.2f LV %.2f\n',gains(1),gains(2),gains(3),gains(4));
        end
    end
%     if strmatch(et,'Dc')
%         ors = cat(2,Expt.Trials.or);
%     end
        Expts{nexpts} = Expt;
    nexpts = nexpts+1;
    end
%  max(spkids);
end
    
function  gains = MakeEyeCalfile(cExpt,emname);
 
emfile = strrep(cExpt.Header.Name,'.mat','.em.mat');
if ~exist(emfile,'file')
    gains = NaN;
    return;
end
load(emfile);
if ~isfield(Expt.Trials,'fx')
    for j = 1:length(cExpt.Trials)
        id = find([Expt.Trials.id] == cExpt.Trials(j).id);
        Expt.Trials(id).fx = cExpt.Trials(j).fx;
        Expt.Trials(id).fy = cExpt.Trials(j).fy;
    end
end
[gains, positions] = EyeCal(Expt,'ids',[cExpt.Trials.id]);
save(emname,'gains','positions');


function d = CreationDate(Text)

did = strmatch('uf',Text.text);
if isempty(did) %online file
    did = strmatch('bt',Text.text);
end
if length(did)
    ds = Text.text(did(1),:);
    did = strfind(ds,'Creat');
    if length(did)
        d = datenum(ds(did(1)+8:end));
    else
        d = 0;
    end
else 
    d = 0;
end

function Idx = AddError(err,Idx)

if ~isfield(Idx,'errs') | isempty(strmatch(err,Idx.errs))
    msgbox(err);
    Idx.errs = {Idx.errs{:} err};
    Idx.newerrs = Idx.newerrs+1;
else
    fprintf([err '\n']);
end


function name = BuildName(name)
if isempty([strfind(name,'/') strfind(name,'\')])
    name = [pwd '/' name];
end
id = strfind(name,':');
if id
    name = name(id(1)+1:end);
else
    name = name;
end
if isunix
    name = strrep(name,'\','/');
end


