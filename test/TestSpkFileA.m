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
state.showerrs = 1;
saveexpts = 0;
state.needframes = 1;
state.tt = [];
state.nospikes = 0;
state.alltrials = 0;
state.profiling = 0;
errs = {};

mkidx =0;  %%need this up here so taht relist works


if isdir(name)
    [Expt, Expts] = ReadExptDir(name, varargin{:});
    return; 
end

j = 1;
while j <= length(varargin)
    vg = varargin{j};
    if strncmpi(vg,'alltrials',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'setprobe',4)
        j = j+1;
        setprobe = varargin{j};
    elseif strncmpi(vg,'findprobes',6)
       mkidx = 1;
       setprobe = -1; %force relisting of probes
    elseif strncmpi(vg,'noerrs',5)
        state.showerrs = 0;
    elseif strncmpi(vg,'nospikes',5)
        state.nospikes = 1;
    elseif strncmpi(vg,'rfs',3)
       state.nospikes = 2;
    elseif strncmpi(vg,'bysuffix',7)
       ignoreSpikeO = 2;
    elseif strncmpi(vg,'method',6)
        j = j+1;
        state.method = varargin{j};
    elseif strncmpi(vg,'noframes',6)
        state.needframes = 0;
    elseif strncmpi(vg,'profile',6)
        state.profiling = 1;
    elseif sum(strncmpi(vg,{'quicksuffix' 'quickload'},9))
       ignoreSpikeO = 2;
       quickload = 1;
       if strcmp(vg,'quicksuffix')
           quickload = 2;
       end
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    elseif strncmpi(vg,'saveexpts',6)
        saveexpts = 1;
    elseif strncmpi(vg,'sortexpts',6)
        idxfile = strrep(name,'.mat','idx.mat');
        load(idxfile);
        [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header,1, Expt, state);
        return;
    elseif strncmpi(vg,'usealltrials',8)
        state.alltrials = 1;
    elseif strncmpi(vg,'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
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
tic; 

inexpt = 0;
Stimulus.OptionCode = '+se';
nfpj=0;
bsctr = 0;


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
    
    readmethod = 0;




waitloop = round(length(aText.text)/20); %update waitbar 20 times
%waitbar(0,sprintf('parsing %d txt lines',length(aText.text)));
for j = 1:length(aText.text)
%    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);

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
        elseif aText.codes(j,4) == 2  % this was FROM spike 3
            if aText.codes(j,1) == 3 && instim == 1 %end stim
                instim = 2;
            end

        elseif aText.codes(j,1) == 5 %stim start
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
        elseif aText.codes(j,1) == 3 %end stim
            instim = 2;
        elseif aText.codes(j,1) == STARTEXPT %end stim
            inexpt = 1;
            if readmethod == 1
                id = find([Expts.end] < aText.times(j));
                ix = length(id)+1;
            end
        elseif aText.codes(j,1) == ENDEXPT %end stim
            if inexpt
                ix = ix+1;
            end
            inexpt = 0;
        elseif regexp(txt,'sb[+,-,0]') %ignore these lines
        elseif strncmp(txt,'{}',2) %bug!!
        txt = aText.text{j};
            
        elseif strncmp(txt,'EndStim',7) %finished reading all text related to last stim
            gotend = 1;
        elseif strncmp(txt,'manexpt=',8)
            Stimulus.manexpt = txt(9:end);
        elseif strncmp(txt,'manexvals',8)
            stimid = sscanf(txt,'manexvals%d');
            id = strfind(txt,' ');
            if ~isempty(id)
                a = sscanf(txt(id(1)+1:end),'%f');
            end
            Stimulus.exvals = a;
        elseif strncmp(txt,'exvals',6)
            a = sscanf(txt,'exvals %f %f %f %d');
            if isfield(Stimulus,'et')
                Stimulus.(Stimulus.et) = a(1);
            end
            if isfield(Stimulus,'e2')
                if sum(strcmp(Stimulus.e2,{'backMov' 'Dc'}))
                    Stimulus.(Stimulus.e2) = a(2);
                end
            end
            Stimulus.ex3val = a(3);
            if isfield(Stimulus,'e3')
                if strcmp(Stimulus.e3,'mixac')
                    Stimulus.(Stimulus.e3) = a(3);
                end
            end
        elseif strncmp(txt,'mixac',5)
            Stimulus.mixac = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'Off at',6) %Storage turned off - should be outside trial
            if intrial
                fprintf('Storage Off in Trial at %.1f',atext.times(j));
            end
        elseif sum(strncmp(txt,{'Write On'},8)) %lines to ignore
            
        elseif sum(strncmp(txt,{'RightHemi' 'Electrode' ' Electrode' 'Experime'},8)) %lines to include in Comments
            txtid = [txtid j];
        elseif strncmp(ss,'cx',2)
            Stimulus.cx = txt(3:end);
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
            if length(Stimulus.Phaseseq) > length(Stimulus.Stimseq)
                Stimulus.Phaseseq = Stimulus.Phaseseq(1:length(Stimulus.Stimseq));
            end
            if trial > length(Trials.Start)
                Trials.Phaseseq{trial} = Stimulus.Phaseseq;
            elseif instim ~= 1  && trial > 1 && t < Trials.Start(trial)
                Trials.Phaseseq{trial-1} = Stimulus.Phaseseq;
            end
            if Stimulus.id == 6136 || Stimulus.id > 570
                trial;
            end
        elseif strncmp(txt,'mtFl=',5) 
            if length(txt) > 5
                Stimulus.mtFl = sscanf(txt(6:end),'%d');
            end                
        elseif strncmp(txt,'mtFn=',5) && length(txt) > 50
            framets = sscanf(txt(6:end),'%f');
            id = find(diff(framets(1:end-1)) > 1.5);
            if diff(framets(end-1:end)) > 2.8
                id = cat(1,id ,length(framets));
            end
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.rptframes{trial-1} = id;
                    Trials.ndrop(trial-1) = length(id);
                    Trials.framet{trial-1} = framets;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.framet{trial} = framets;
                    Trials.ndrop(trial) = length(id);
                    Trials.rptframes{trial} = id;
                end
                
            end
        elseif strncmp(txt,'mtrX=',5)
            Stimulus.xoseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.xoseq{trial-1} = Stimulus.xoseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.xoseq{trial} = Stimulus.xoseq;
                end
                
            end
        elseif strncmp(txt,'mtrY=',5)
            Stimulus.yoseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.yoseq{trial-1} = Stimulus.yoseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.yoseq{trial} = Stimulus.yoseq;
                end
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
        elseif strncmp(txt,'mtxo=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.midtrial] > trial);
            if ix
                ix = ix(1);
                Expts(ix).xovals = sscanf(txt(6:end),'%f');            
            end
            end
        elseif strncmp(txt,'dx:',3)
            Trials.dxseq{trial} = sscanf(txt(4:end),'%f');
        elseif strncmp(txt,'mtei=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.midtrial] > trial);
            if ix
                ix = ix(1);
                Expts(ix).e1vals = sscanf(txt(6:end),'%f');            
            else
                err = sprintf('No Expt for mtei at trial %d',trial);
                Expt = AddError(Expt, err);
                ix = lastix;
            end
            end
        elseif strncmp(txt,'mte3=',5)
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
            ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
            Expts(ix).e2vals = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'Off at',5)
        elseif strncmp(txt,'EndExpt',5)
            if inexpt
                ix = ix+1;
            end
            inexpt = 0;
        elseif strncmp(txt,'sonull',5)
        elseif strncmp(txt,'NewConnect',7)
        elseif strncmp(txt,'BGCS Version',7)
        elseif strncmp(txt,'testflag',7)
        elseif strncmp(txt,'rptframes ',3)
            Stimulus.rptframes = sscanf(txt(10:end),'%d');
        elseif strncmp(txt,'ijump',5)
            Stimulus.ijump = sscanf(txt(6:end),'%d');
        elseif strncmp(txt,'nph',5)
            Stimulus.nph = sscanf(txt(6:end),'%d');
        elseif strncmp(txt,'seof',4)
            Stimulus.seof = sscanf(txt(5:end),'%d');
        elseif strncmp(txt,'/local',6)
            Stimulus.imprefix = txt;
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
            if txt(3) == '='  %char description
                Stimulus.OptionCode = txt(4:end);
            elseif ~isempty(a)
                Stimulus.op = a(1);
            end
            if length(a) > 1
                Stimulus.optionb = a(2);
            end
            if isempty(Stimulus.op)
                fprintf('Missing op stim %d\n',trial);
                Stimulus.op = 0;
            end
        elseif strncmp(txt,'annTyp',6)
            Stimulus.annTyp = sscanf(txt(7:end),'%f');
        elseif ss(1) == '#'
            comment = ss;
        elseif ~isstrprop(ss(1),'alphanum')
            fprintf('Non-Printing Name %s\n',txt);
        elseif strncmp(txt,'vs',2) || strncmp(txt,'sq',2)
            try
                [a, ok] = sscanf(val,'%f');
                if ~ok
                    Stimulus.(ss) = val;
                else
                    Stimulus.(ss) = a(1);
                    if length(a) > 2
                        Stimulus.FlipDir = a(3);
                    end
                end
            catch
                ok = 0;
                Stimulus.FlipDir = 0;
            end
        elseif strncmp(txt,'fp',2)
            a = sscanf(val,'%f');
            if length(a) > 1
            if instim == 1
                nfpj = nfpj+1;
                Stimulus.dfx(nfpj) = a(1);
                Stimulus.dfy(nfpj) = a(2);
            else
                Stimulus.fx = a(1);
                Stimulus.dfx = a(1);                
                Stimulus.fy = a(2);
                Stimulus.dfy = a(2);
                nfpj = 0;
            end
            end
        elseif strncmp(txt,'fx',2) || strncmp(txt,'fy',2)
%if really in a stimulus, make note of new fx but keep original also
%if instim ==2, don't change anything.
            a = sscanf(val,'%f');
            if instim == 1
                Stimulus.(['d' ss]) = a;
            elseif instim == 0
                Stimulus.(ss) = a;
                Stimulus.(['d' ss]) = a;
                lastfix.(ss) = a;
            elseif instim == 2 %Post stim, but trial counter not yet incremented. Store value
                Stimulus.(ss) = a;
                Stimulus.(['d' ss]) = a;
                lastfix.(ss) = a;
            end
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
    if aText.codes(j,1) == FRAMESIGNAL
            newstart = 1;
    else
        newstart = 0;
    end
    if (aText.codes(j,4) == 2 || isempty(txt))  && exist('Stimulus','var') %empty text = code from spike2 -> binoc
        correctdir = 0;
        if isfield(Stimulus,'OptionCode') && ...
                (~isempty(strfind(Stimulus.OptionCode,'+2a')) || ~isempty(strfind(Stimulus.OptionCode,'+afc'))) ...
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
%            Stimulus.FlipDir = 1;
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
       if isfield(Trials,'FlipDir')
           flipdir = Trials.FlipDir(trial-1);
       else
           flipdir = 1;
       end
        if aText.codes(j,1) == WURTZOKW & trial > 1 & tdelay < 10000
            Trials.RespDir(trial-1) = -1 * correctdir.*flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
  %          Trials.score(trial-1) = 0;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == WURTZOK & trial > 1  & tdelay < 10000
            Trials.RespDir(trial-1) = 1 * correctdir.* flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
 %           Trials.score(trial-1) = 1;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == BADFIX
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
    if readmethod == 1 && isfield(Trials,'EndTxt')
        endtimes = Trials.EndTxt;
    else
        readmethod = 0;
        endtimes = Trials.End;
    end
    if length(endtimes) >= trial && t >= endtimes(trial) & instim & trial <= length(Trials.Start) && newstart == 0
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
            Stimulus.inexpt = inexpt;
        Trials = SetTrial(Stimulus, Trials, trial, ntrials);
        Stimulus.CorLoop = 0;
        Stimulus.uf = '';
        Stimulus.fx = lastfix.fx;
        Stimulus.fy= lastfix.fy;
        Stimulus.rptframes = [];
        Stimulus.endevent = aText.codes(j,1);
        if ~isempty(ustimmarkch)
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
        end
        instim = 0;
    elseif trial > length(Trials.Start)
        if gotend
          break;
        end
% if we have passed stim Start time, text must refer to next stimulus. 
%But don't increment based on codes coming back from Spike
    elseif t > Trials.Start(trial) && instim == 0 && aText.codes(j,4) ~= 2
        instim = 1;
    end
end

fprintf('Actual Run Time %.2f\n',toc);
if state.profiling == 1
    profile off;
    profsave(profile('info'),'SpkFile_profile');
    profile viewer;
end

if isfield(Trials,'ve') && iscellstr(Trials.ve)
    for j = 1:length(Trials.ve)
        if strncmp(Trials.ve{j},'binoclean',8)
            x = sscanf(Trials.ve{j}(11:end),'%f');
            ve(j) = 10+x(1);
            if length(x) > 1
                ve(j) = ve(j) + x(2)./100;
            end
        else
            ve(j) = 0;
        end
    end
    Trials.ve = ve;
    if length(Trials.ve) < length(Trials.Start)
        Trials.ve(length(Trials.Start)) = median(ve);
    end
end

if isempty(Expts(end).et) && length(Expts) > 1
    Expts(end).et = Expts(end-1).et;
    Expts(end).e2 = Expts(end-1).e2;
    Expts(end).e3 = Expts(end-1).e3;
end

if ~isfield(Expts,'result') || isempty(Expts(end).result)
    Expts(end).result = 0;
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
fprintf('Clearing %d cells took %.2f\n',length(cellids),toc);


function Trials = SetTrial(Stimulus, Trials, trial, ntrials)

if isfield(Stimulus,'op')
Trials.op(trial) = Stimulus.op;
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
