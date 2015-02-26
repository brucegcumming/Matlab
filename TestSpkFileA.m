function [Expt, Expts, AllData, Raw] = APlaySpkFile(name, varargin)

% [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%
% Builds/retnes lists of trials/expts from a matlab file
% APlaySpkFile(name, 'relist')  Re-builds the index filt
% APlaySpkFile(name, 'setprobe', n)  sets which probes data is loaded
% initially
% APlaySpkFile(name, 'setprobe', -1) rebuilds probe list
% APlaySpkFile(name, 'usealltrials') includes bad fix trials too
%
%Error Fixing A file nameAdd.txt (i.e. jbeG044.mat -> jbeG044Add.txt)
%If it exists is read in and adds lines to the text stream.

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
stimchange = [];
logfid = 0;
fixup = 1;
UstimV = 0;
setprobe = 1;
state.nospikes = 0;
nerr = 0;
preperiod=2000;
postperiod = 2000;
savedvdt = 0;
dvfile = [];
Raw = [];
timeoffset = 0;
ignoreSpikeO = 0;  %% for files with full voltage records elsewhere, ignore these
quickload = 0;
%method for matching up StimChan starts with text events. With many channel
%recordings, can get big delays between these, so the original method,
%finding the nearest stimchan event, is unreliable. New method (1) finds
%all events firts, and if the numbers are the same, matches in order.
%dfeault set to method 1 Jan 6 2012 by bgc
state.method = 1; 
state.showerrs = 0;
saveexpts = 0;
state.needframes = 1;
state.tt = [];
state.nospikes = 0;
state.alltrials = 0;
state.profiling = 0;
errs = {};

mkidx =1;  %%need this up here so taht relist works

j = 1;
while j <= length(varargin)
    vg = varargin{j};
    if strncmpi(vg,'alltrials',4)
        onlinedata = 1;
    elseif strncmpi(vg,'profile',6)
        state.profiling = 1;
    end
    j = j+1;
end
    
clusterdate = now;
thecluster = 1;

if ischar(name)
    if ~exist(name,'file')
        fprintf('No file %s\n',name);
        return;
    end
    
    

   load(name);
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
    elseif strncmpi(vg,'alltrials',8)
        argon = {argon{:} varargin{j}};
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
    elseif strncmpi(vg,'mkufl',4)
        MkUfl(name, Ch30);
    elseif strncmpi(vg,'name',4)
        j = j+1;
        name = varargin{j};
   elseif strncmpi(vg,'noidx',5)
       mkidx = 0;
    elseif strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif strncmpi(vg,'cluster',2)
        j = j+1;
        thecluster = varargin{j};
    elseif strncmpi(vg,'play',4)
        playspikes = 1;
    elseif strncmpi(vg,'rfs',3)
        MkUfl(name, Text,'overwrite');
    elseif strncmpi(vg,'relistonly',8)
        mkidx = 2;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    end
    j = j+1;
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
fprintf('Reading Text took %.2f\n',toc);
end

for j = 1:length(Expts)
    Expts(j).midtrial = (Expts(j).lasttrial+Expts(j).firsttrial)/2;
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
lastfix.fx = 0;
lastfix.fy = 0;
Stimulus.Flag = '';


   vars = who('Ch*');
   for j = 1:length(vars)
       if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
           
            if strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
                bstimes = stimlvl.times(1:2:end);
                estimes = stimlvl.times(2:2:end);
            end
       end
   end
 
   settrials =1;
   Trials.Start = bstimes;
   Trials.End = estimes;
   ntrials = length(Trials.Start);
fprintf('Text->stims .....');
if state.profiling == 1
    profile on;
end
inexpt = 0;
Stimulus.OptionCode = '+se';
nfpj=0;
bsctr = 0;

tic;
id = [];
codes = zeros(1,length(aText.text));
findstrs = {'puA' 'puF' 'USd' 'USp' 'USf' 'nph' 'ijump' 'mixac' 'baddir' ...
        'e1max' 'backMov' 'FakeSig' 'pBlack' 'aOp' 'aPp' 'seof' 'serange' ...
        'nimplaces' 'usenewdirs' 'choicedur' 'cha' 'imi' 'choicedur' 'ePr' ...
        'coarsemm' 'psyv' 'imi'};    
    charstrs = {'Covariate'  'hxtype' 'cx' 'adapter' 'exp'};
    extrastrs = {'StartDepth' 'Electrode'};
    findstrs = [findstrs charstrs extrastrs];
    for j = 1:length(findstrs)
        f = findstrs{j};
        slens(j) = length(f);
        id = find(strncmp(f,aText.text,length(f)));
        codes(id) = j;
        if sum(strcmp(f, charstrs))
            vartypes(j) = 'C';
        elseif sum(strcmp(f, extrastrs))
            vartypes(j) = 'X';
        else
            vartypes(j) = 'N';
        end
    end
    id = find(strncmp('#',aText.text,1));
    j = j+1;
    codes(id) = j;
    vartypes(j) = '#';
    slens(j) = 0; 
    readmethod = 0;



    if readmethod == 1 && isfield(Trials,'EndTxt')
        endtimes = Trials.EndTxt;
    else
        readmethod = 0;
        endtimes = Trials.End;
    end
waitloop = round(length(aText.text)/20); %update waitbar 20 times
%waitbar(0,sprintf('parsing %d txt lines',length(aText.text)));
toc
nc = 0;

tendid = aText.codes(:,4) == 2 & ismember(aText.codes(:,1),[WURTZOKW WURTZOK BADFIX]);
tendid = tendid';
tdelay = 0;
flipdir = 1;
trialendtime = Trials.End(1);


acodes = aText.codes(:,1);
acodes(find(tendid)) = 100;
codes(find(tendid)) = 100;
slens(100) = 0;
vartypes(100) = '!';
newstarts = acodes == FRAMESIGNAL;

tic;
for j = 1:length(aText.text)
    txt = aText.text{j};
    t = aText.times(j);
    acode = acodes(j);
    dcode = aText.codes(j,4);
    if codes(j) > 0
        c = codes(j);
        slen = slens(c);
        vartype = vartypes(c);
        ss = txt(1:slen);
        if txt(slen+1) == '='
            val= txt(slen+2:end);
        else
            val= txt(slen+1:end);
        end
    else
        id = findstr(txt,'=');
        if ~isempty(id)
            slen = id(1)-1;
            val= txt(slen+2:end);
        elseif length(txt) > 2
            slen = 2;
            val= txt(slen+1:end);
        else
            slen =0;
            val = '';
        end
        if ~isempty(txt)
            ss = txt(1:slen);
        end
        vartype = 'N';
    end
    if slen > 0
        if vartype == 'X'
        elseif dcode == 2  % this was FROM spike 3
            if acode == 3 && instim == 1 %end stim
                instim = 2;
            end
            
        elseif acode == 5 %stim start
            instim = 1;
            Stimulus.Seedseq = {};  %% these must be set for each stim
            Stimulus.Stimseq = {};
            Stimulus.xoseq = {};
            Stimulus.yoseq = {};
            Stimulus.Phaseseq = [];
            Stimulus.cLseq = [];
            Stimulus.cRseq = [];
            gotend = 0;
            nfpj=0;
            bsctr = bsctr+1;
        elseif acode == 3 %end stim
            instim = 2;
        elseif acode == STARTEXPT %end stim
            inexpt = 1;
            if readmethod == 1
                id = find([Expts.end] < aText.times(j));
                ix = length(id)+1;
            end
        elseif acode == ENDEXPT %end stim
            if inexpt
                ix = ix+1;
            end
            inexpt = 0;
        elseif vartype == 'C'
            Stimulus.(ss) = val;
        else
            %            if strmatch(ss,{'0' '1' '2' '3' '4'})
            if regexp(ss,'^[0-9]')
                ss = ['x' ss];
            end
            
            try
                [Stimulus.(ss), ok] = sscanf(val,'%f');
                if ~ok
                    Stimulus.(ss) = val;
                end
            catch
                ok = 0;
            end
            if strncmp(ss,'et',2)
                Stimulus.(ss) = val;
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                if(ix) Expts(ix).et = val; end
            end
            if strncmp(ss,'e2',2) & length(ix)==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e2 = val;
            end
            if strncmp(ss,'e3',2) & length(ix) ==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e3 = val;
            end
        end
    end

    %if  (aText.codes(j,4) == 2 || (isempty(txt) && acode) > 0)  %empty text = code from spike2 -> binoc
    if codes(j) == 100
        if acodes(j) == 100
        correctdir = 0;
        nc = nc+1;

%Real ON/Off times are set from the events above. But need to know that
%text following WURTZOK applies to the next stimulus. So instim = 2 means
%text has been received ending trial, but not officially over yet. 
%trial gets incremented at end stim. So response applies to trial -1
%but check that the time is sensible. Can get one of these events when
%storage is off, resetting the last stored trials
        if acode == WURTZOKW & trial > 1 & tdelay < 10000
            Trials.RespDir(trial-1) = -1 * correctdir.*flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
  %          Trials.score(trial-1) = 0;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif acode == WURTZOK & trial > 1  & tdelay < 10000
            Trials.RespDir(trial-1) = 1 * correctdir.* flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
 %           Trials.score(trial-1) = 1;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif acode == BADFIX
            if trial <= length(Trials.End)
            if aText.times(j) > Trials.End(trial) && settrials== 0
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
        else %shouldb't happen
            trial = trial;
       end
        end
    end

    if t >= trialendtime && newstarts(j) == 0 && instim 
%need to read past the end of the last trial a litle way to get things like
%Stimseq which come afterwards

        if isfield(Stimulus,'st')
            Stimulus.inexpt = inexpt;
            %Trials = SetTrial(Stimulus, Trials, trial, ntrials);
            Stimulus.CorLoop = 0;
            Stimulus.uf = '';
            Stimulus.fx = lastfix.fx;
            Stimulus.fy= lastfix.fy;
            Stimulus.rptframes = [];
            Stimulus.endevent = acode;
            if 0 %~isempty(ustimmarkch)
                marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1 ...
                    & bitand(ustimmarkch.codes(:,1),1));
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
            if length(endtimes) >= trial
                trialendtime = endtimes(trial);
            end
        end
        instim = 0;
    elseif trial > length(Trials.Start)
        if gotend
          break;
        end
% if we have passed stim Start time, text must refer to next stimulus. 
%But don't increment based on codes coming back from Spike
    elseif t > Trials.Start(trial) && instim == 0 && dcode ~= 2
        instim = 1;
    end

end

if state.profiling == 1
    profile viewer;
end
fprintf('Text->stims %.2f\n',toc);

function Trials = SetTrial(Stimulus, Trials, trial, ntrials)

if isfield(Stimulus,'op')
Trials.op(trial) = Stimulus.op;
end
if isfield(Stimulus,'optionb')
Trials.optionb(trial) = Stimulus.optionb;
end
return;
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
    if strcmp(F,'St')
        Trials.St(trial) = Stimulus.(F);
    elseif strncmp(F,'trode',5)
    elseif ischar(Stimulus.(F)) 
        if ~isfield(Trials,F) || ischar(Trials.(F)) || iscell(Trials.(F))
        if ~isempty(Stimulus.(F))
            Trials.(F){trial} = Stimulus.(F);
        end
        else
            fprintf('%s is Char in Stimulus, not in Trials\n',fn{k});
        end
    elseif sum(strcmp(F,{'Seedseq' 'Stimseq' 'Phaseseq' 'cLseq' 'cRseq' 'xoseq' 'yoseq' 'rptframes' 'rwtimes'}))    
           Trials.(fn{k}){trial} = Stimulus.(fn{k});
    else
        if ~isfield(Trials,F)
            Trials.(F)(1:ntrials,1) = NaN; % pre-allocate memory
        end
        if isempty(Stimulus.(F))
            Trials.(F)(trial,1) = NaN;
        else
%            Trials.(F)(trial,1:length(Stimulus.(F))) = Stimulus.(F);
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
%    Trials.FlipDir(length(Trials.Start)) = 1;
end


function ix = FindExptn(Expts, readmethod, bsctr, trial, ix)
if readmethod == 1
    trial = bsctr;
end
    
if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
    lastix = ix;
    ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
    if ix
        ix = ix(end);
    else
        ix = lastix;
    end
end
