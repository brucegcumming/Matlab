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
state.tt = TimeMark(state.tt,'Start');
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

[aText, Text] = AddText(regexprep(name,'\.[0-9]*.mat','Add.txt'), aText, Text);
[aText, Text] = AddText(strrep(name,'.mat','Add.txt'), aText, Text);

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
fprintf('Text->stims .....');
waitbar(0,sprintf('parsing %d txt lines',length(aText.text)));
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
        'coarsemm' 'psyv' };
    charstrs = {'Covariate'  'hxtype' 'cx' 'adapter' 'exp'};
    findstrs = [findstrs charstrs];
    for j = 1:length(findstrs)
        f = findstrs{j};
        slens(j) = length(f);
        id = find(strncmp(f,aText.text,length(f)));
        codes(id) = j;
        if sum(strcmp(f, charstrs))
            vartypes(j) = 'C';
        else
            vartypes(j) = 'N';
        end
    end

for j = 1:length(aText.text)
%    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);

    if codes(j) > 0
        c = codes(j);
        slen = slens(c);
        vartype = vartypes(c);
        ss = txt(1:slen);
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
    end
end

toc
waitloop = round(length(aText.text)/20); %update waitbar 20 times
for j = 1:length(aText.text)
    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);

%    id = find(txt == 0);
%    txt = txt(1:id(1)-1);
    if length(txt)>1;
        ss = txt(1:2);
        slen = 2;
        id = findstr(txt,'=');
        if ~isempty(id)
            slen = id(1)-1;
            ss = txt(1:slen);
        elseif 0
            for f = {'puA' 'puF' 'USd' 'USp' 'USf' 'nph' 'ijump' 'mixac' 'baddir' ...
                    'e1max' 'backMov' 'FakeSig' 'pBlack' 'aOp' 'aPp' 'seof' 'serange' ...
                    'nimplaces' 'usenewdirs' 'choicedur' 'cha' 'imi' 'choicedur' 'ePr' ...
                    'coarsemm' 'psyv' }
                if sum(strncmp(txt,f,length(f{1})))
                    ss = txt(1:length(f{1}));
                    slen = length(f{1});
                end
            end
        end
        vartype = 'N';
        for f = {'Covariate'  'hxtype' 'cx' 'adapter' 'exp'}
            if sum(strncmp(txt,f,length(f{1})))
                ss = txt(1:length(f{1}));
                slen = length(f{1});
                vartype = 'C';
            end
        end
        if length(txt) > 2 & txt(slen+1) == '='
            val= txt(slen+2:end);
        else
            val= txt(slen+1:end);
        end
        if aText.codes(j,4) == 2  % this was FROM spike 3
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
            a = InterpretLine(txt);
            if (isfield(a,'electrode') && ~strcmp(a.electrode,'default')) || ~isfield(Peninfo,'electrode')
                Peninfo = CopyFields(Peninfo, a);
                Peninfo.trode = txt;
            end
            id = strfind(txt,'Contact');
            if length(id)
                x = id(1);
                id = strfind(txt(id:end),' ');
                Peninfo.probesep = sscanf(txt(id(1)+x:end),'%d');
            end
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
    if mod(j,waitloop) == 0
        waitbar(j/length(aText.text));
    end
end
delete(waitbar(1));
drawnow;
if state.profiling == 1
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
Expt.Comments.Peninfo.trode = BuildProbeString(Peninfo);

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
fprintf('Text->stims %.2f\n',toc);
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
for trial = 1:length(Trials.Start)
    if isnan(Trials.End(trial)) && trial > length(Trials.End)
        Trials.End(trial) = Trials.Start(trial) + nomdur;
    end
    spkids = find(Spks.times > Trials.Start(trial)-preperiod & Spks.times < Trials.End(trial)+postperiod);
    Trials.Spikes{trial} = (Spks.times(spkids) - Trials.Start(trial));
    if spkids
        spkid(trial,:) = [spkids(1) spkids(end)];
        Trials.Cluster{trial} = Spks.codes(spkids,1);
        frames = [];
 %this must be just for caliration of timing. Can't see need for this
 %if we recording real data
        if frametimes & testframe
            for k = spkids'
                id = find(frametimes < Spks.times(k));
                if id
                    frames = [frames frametimes(id(end))-Spks.times(k)];
                end
            end
        end
        if isempty(frames)
            Trials.Frames(trial) = NaN;
        else
            Trials.Frames(trial) = frames(1);
        end
    else
        spkid(trial,:) = [0 0];
        Trials.Frames(trial) = NaN;
        Trials.Cluster{trial} = [];
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
trial = 1;
for spk = 1:size(Spks.times)
    t = Spks.times(spk);
    if t > (Trials.End(trial)+postdur) & instim
        Trials.Spikes{trial} = (Spks.times(lastspk:maxspk) - Trials.Start(trial));
        spkid(trial,:) = [lastspk, maxspk];%
%        Trials(trial).Frames = frames;
        Trials.Cluster{trial} = Spks.codes(lastspk:maxspk,1);
        trial = trial+1;
        instim = 0;
        lastspk = spk;
        elseif t > Trials.Start(trial)
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
        if trial > length(Trials.Start)
            break;
        end
end
fprintf('Spikes Took %.3f\n',toc);
fixfile = strrep(idxfile,'idx','fix'); 
if exist(fixfile)
    load(fixfile);
    f = fields(fixTrials);
    [tid, fid] = ismember(fixTrials.id,Trials.id);
    tid = find(tid);
    for j = 1:length(f)
        if ~strcmp(f{j},'id')
            if ~isfield(Trials,f{j})
                Trials.(f{j}) = ones(size(Trials.Start));
            end
            Trials.(f{j})(fid(tid)) = fixTrials.(f{j})(tid);
        end
    end
end
%Trials = Trials(2:end);
AllData.datenum = clusterdate;
AllData.quickload = quickload;
AllData.Events = Events;
AllData.Text = aText;
AllData.Spikes = Spks;
if exist('Ch32','var') && strcmp(Ch32.title,'DigMark')
    Expt.DigMark = Ch32; %So its saved
end


stored = find(bitand([Trials.op],STOREBIT));
% bitand STOREBIT (16) identifies if storage was on;
%Trials = Trials(stored);

%Expts #20 seems to be empty, but storage was on...
if isfield(Trials,'TrueEnd') & length(Trials.TrueEnd) < length(Trials.End)
    Trials.TrueEnd(length(Trials.End)) = NaN;
end
if isfield(probes,'var')
for j = 1:length(probes)
   probes(j).var = strrep(probes(j).var,'A','');
end
end
Expt.Trials = Trials;
Expt.Spkid = spkid;
Expt.Probes = probes;
Expt.bstimes= bstimes;
Expt.estimes= estimes;
Expt.ExptList = Expts;
if exist('mainsch','var') && isfield(mainsch,'times')
    Expt.mainstimes = mainsch.times * 10000;
end

Header.Name = BuildName(name);
Header.Spike2Version = version;
Header.unstored = nonstore;
Header.CreationDate = CreationDate(Text);
Header.ReadMethod = readmethod;
if isfield(Expt,'DataType')
    Header.DataType = Expt.DataType;
end
clear Text;
if isempty(Expts) %nothing in this file
    fclose(logfid);
    return;
end

if mkidx == 1
        fprintf('Saving index %s\n',idxfile);
        if logfid > 0
        fprintf(logfid, '%s Saving index %s\n',datestr(now),idxfile);
        end
        Expt.bstimes = bstimes;
        Expt.estimes = estimes;
        Expt.Header = Header;
    save(idxfile,'Expt','Expts');
    iExpts = Expts;
end

Trials = AddStimsFromFile(strrep(name,'.mat','.SetStim'),Trials);


if isfield(Expts,'firsttrial')
    if timeoffset
        args = {args{:} 'timeoffset' timeoffset};
    end
[Expts, Expt] = SortExpts(Expts, Trials, Header, thecluster, Expt, state);
Expts = AddComments(Expts,Expt);
end

if mkidx == 1
    tExpts = Expts;
    ExptList = MkExList(Expts);
    Expts = iExpts;
    Expt.starttrial = starttrial;
    save(idxfile,'Expt','Expts','ExptList');
    Expts = tExpts;
end
if saveexpts
    SaveExpts(name, Expts);
end
if mkidx
%    save(idxfile,'Expt','Expts');
end
fclose(logfid);
WriteErrors(idxfile, Expt);


function  StimCh = AddStimsFromFile(AddTxtFile, StimCh)
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%f %f','delimiter','\n'); %ON, Off
    for j = 1:length(a{1})
        id = find(StimCh.times < a{1}(j));
        id = id(end);
        lvl = cat(1,StimCh.level(1:id), 1, 0 ,StimCh.level(1+id:end));
        t = cat(1,StimCh.times(1:id), a{1}(j), a{2}(j), StimCh.times(1+id:end));
        StimCh.level = lvl;
        StimCh.times = t;
    end
end

function Trials = AddStimValsFromFile(name, Trials)
fid = fopen(name,'r');
if fid > 0
    mycprintf('blue','Reading Stimulus Properties from %s\n',name);
    a = textscan(fid,'id%d %s','delimiter','\n');
    tid = a{1};
    s = a{2};
    for j = 1:length(a{1})
        if strncmp('badexpt',s{j},7) %error in file - ignore this expt
            Trials.baddexpts(tid(j)) = 1;
        else
        id = find(Trials.id  == a{1}(j));
        if length(id) ==1
            x = sscanf(s{j}(4:end),'%f');
            if strncmp(s{j},'dx:',3)
                Trials.dxvals{id} = x;
            elseif strncmp(s{j},'ce:',3)
                Trials.cevals{id} = x;
            end
        end
        end
    end
end

function SaveExpts(name, Expts)
%make separate files for each expt on disk so that can access in parallel
%when needed
for j = 1:length(Expts)
    outfile = strrep(name,'.mat',['Expt' num2str(j) '.mat']);
    Expt = Expts{j};
    save(outfile,'Expt');
end

function str = BuildProbeString(P)

str = '';
if isfield(P,'electrode')
    str = [str 'Electrode ' P.electrode ' '];
end
if isfield(P,'tube')
    str = [str 'Tube ' P.tube ' '];
end
if isfield(P,'hemishpere')
    str = [str P.hemisphere 'HemiSphere '];
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

function WriteErrors(idxfile, Idx)
if isfield(Idx,'errs') & length(Idx.errs)
    ename = strrep(idxfile,'idx.mat','err.txt');
    fid = fopen(ename,'a');
    fprintf(fid,'%s\n',Idx.errs{:});
    fclose(fid);
end

function [aText, Text] = AddText(AddTxtFile, aText, Text)
%Add text reads lines of the form
%t txt
%where t is in seconds, not timestamps
%Addts these to the txt record.
%
%can also add Stimulus Values by id with
%idxxxx str

[a,b,c] = fileparts(AddTxtFile);
if ~strcmp(c,'.txt')
    return;
end
SpkDefs;
fid = fopen(AddTxtFile,'r');
ts = now;

if fid > 0
    cprintf('blue','Adding Text from %s\n',AddTxtFile)
    a = textscan(fid,'%d %s','delimiter','\n');
    if isempty(a{1})
        if isempty(aText.text)
            return;
        end
        b = textscan(fid,'%2s%d %s','delimiter','\n');
        frewind(fid);
        idid = find(strncmp('id',aText.text,2));
        for j = 1:length(idid)
            tid(j) = sscanf(aText.text{idid(j)},'id%d');
        end
        ti = b{2};
        for j = 1:length(ti)
            id = find(tid == ti(j));
            if ~isempty(id)
                t(j) = aText.times(idid(id(end)));
            else
                t(j) = NaN;
            end
        end
        s = b{3};        
    else
    t = a{1};
    s = a{2};
    end
    fclose(fid);
    
    newlines = 0;
    did = find(strncmp(s,'delete',6));
    if ~isempty(did) && ~isempty(aText.times)
        for j = 1:length(did)
            cprintf('blue','%d deleting %s\n',t(j),aText.text{t(did(j))});
        end
        id = setdiff(1:length(aText.times), t(did));
        aText.text = aText.text(id);
        aText.times = aText.times(id);
        aText.codes = aText.codes(id,:);
    end
    for j = 1:length(s)
        %?why do we look for whitespace?? removed Aug 2010.
        id = findstr(s{j},' ');
        id = [];
        if length(id)
            txt = s{j}(id(1)+1:end);
        else
            txt = s{j};
        end
        if j < 500
            cprintf('blue','%d Adding Text %s\n',t(j),s{j})
        end
        if ~isempty(aText.times)
        if t(j) < 0  && t(j) > -1000%special case for fixing lines
            if strncmp(s{j},'cm=rf',5)
                id = strmatch('cm=rf',aText.text);
                for k = 1:length(id)
                    aText.text(id,:) = s(j);
                end
            end
        elseif t(j) == -1000  %put these at the end
            aText.text = {aText.text{:} txt};
            aText.times = [aText.times; max(aText.times)+1];
            aText.codes = [aText.codes; 0 0 0 0];
        elseif t(j) == 0 || t(j) <= aText.times(1)
            aText.text = {txt aText.text{:}};
            aText.times = [0; aText.times];
            aText.codes = [0 0 0 0; aText.codes];
        elseif sum(strcmp(s(j),{'delete' 'badexpt'})) %special lines not going into text
            txt = '';
        else
            aText.text{end+1} = txt;
            aText.times(end+1)=t(j);
            aText.codes(end+1,:)=0;
            if strcmp(txt,'EndExpt')
                aText.codes(end,1) = ENDEXPT;
            end              
            newlines = newlines+1;
        end
        end
        if ~isempty(txt)
            Text.text(end+1,1:length(txt)) = txt;
        end
    end
    if newlines
        [t, tid] = sort(aText.times);
        aText.text = aText.text(tid);
        aText.times = aText.times(tid);
        aText.codes = aText.codes(tid,:);
    end
    fprintf('Took %.2f\n',mytoc(ts));
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
                    Expts{j} = ParseExptComment(Expts{j}, Expt.Comments.text{cid(bid(k))});
                end
            end
        end
        Expts{j}.Comments.text = {Expt.Comments.text{cid}};
        Expts{j}.Comments.times = {Expt.Comments.times(cid)};
        cid = strmatch('cm=VisualArea',Expt.Comments.text);
        for k = 1:length(cid)
            Expts{j}.Header.Area{k} = Expt.Comments.text{cid(k)}(15:end);
        end
        if isfield(Expt.Comments,'Peninfo')
            Expts{j}.Comments.Peninfo = Expt.Comments.Peninfo;
        end
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
            Trials.(F)(trial,1) = Stimulus.(F)(1);
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

function [Expts, Idx, state] = SortExpts(AllExpts, AllTrials, Header, thecluster, Idx,state,  varargin)
SpkDefs;
timeoffset = 0;
%stimnames = {'None', 'Gabor', 'RDS', 'Grating', 'bar', 'circle', 'rectangle', 'test', 'square', 'probe', '2grating', 'Cylinder', 'twobar', 'rls', 'annulus', 'rdssine', 'nsines'};
Expts = [];
state.tt = TimeMark(state.tt,'Start Sorting');
spikid = Idx.Spkid;
if ~isfield(Idx,'newerrs')
    Idx.newerrs = 0;
end
if isfield(Header,'frameperiod')
frameperiod = Header.frameperiod;
else
frameperiod = 167;
end

if isfield(Header,'Name') && ~isempty(regexp(Header.Name,'[0-9]\.[0-9]*\.mat'))
    Header.bysuffix = 1;
else
    Header.bysuffix = 0;
end

findtrial = 0;
if state.alltrials
    usebadtrials = 1;
else
usebadtrials = 0;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'alltrials',5)
        usebadtrials = 1;
    elseif strncmpi(varargin{j},'findtrial',5)
        j = j+1;
        findtrial = varargin{j};
    elseif strncmpi(varargin{j},'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    j = j+1;
end

if isfield(AllTrials,'ve') && length(AllTrials.ve) < length(AllTrials.Start)
        AllTrials.ve(length(AllTrials.Start)) = median(AllTrials.ve);
end


for j = 1:length(AllExpts)
    if isempty(AllExpts(j).result)
        AllExpts(j).result = -1;
    end
end

AllExpts = AllExpts([AllExpts.result] >= 0);
STIM_BAR = 4;
STIM_GRATING = 3;
for j = 1:length(stimnames)
    if strmatch(stimnames{j},'bar')
        STIM_BAR = j-1; %starts with 0
    end
end

fn = fieldnames(AllTrials);
% make a list of fileds that are NOT automatically set 
state.tt = TimeMark(state.tt,'Fixed Names');
ids = find(strcmp('Spikes',fn));
ids = [ids find(strcmp('OptionCode',fn))];
ids = [ids find(strcmp('Start',fn))];
ids = [ids find(strcmp('End',fn))];
ids = [ids find(strcmp('EndTxt',fn))];
ids = [ids find(strcmp('Result',fn))];
ids = [ids find(strcmp('op',fn))]; % has op and optionb
ids = [ids find(strcmp('Stimseq',fn))];
ids = [ids find(strcmp('Seedseq',fn))];
ids = [ids find(strcmp('rwtimes',fn))];
ids = [ids find(strcmp('xoseq',fn))];
ids = [ids find(strcmp('yoseq',fn))];
ids = [ids find(strcmp('rptframes',fn))];
ids = [ids find(strcmp('framet',fn))];
ids = [ids find(strcmp('cLseq',fn))];
ids = [ids find(strcmp('cRseq',fn))];
ids = [ids find(strcmp('Phaseseq',fn))];
ids = [ids find(strcmp('Cluster',fn))];
ids = [ids find(strcmp('Frames',fn))];
ids = [ids find(strcmp('StartEv',fn))];
ids = [ids find(strcmp('Events',fn))];
ids = [ids find(strcmp('rws',fn))];
ids = [ids find(strcmp('endelay',fn))];
ids = [ids find(strcmp('rwset',fn))];
ids = [ids find(strcmp('estimes',fn))];
ids = [ids find(strcmp('uStimt',fn))];
ids = [ids find(strcmp('Flag',fn))];
ids = [ids find(strcmp('bsstimes',fn))];
ids = [ids find(strcmp('esstimes',fn))];
ids = [ids find(strcmp('ex3val',fn))];
ids = [ids find(strcmp('bsdelay',fn))];
%ids = [ids strmatch('imver',fn)];
%ids = [ids strmatch('imse',fn)];
state.tt = TimeMark(state.tt,'Set Fields');


%do not include PhaseSeq here - has special cases
seqstrs = {'dxvals' 'cevals'};
cellfields = {'rwtimes'};
seqvars = {};
for j = 1:length(seqstrs)
    id = find(strcmp(seqstrs{j},fn));
    if ~isempty(id)
        ids = [ids id];
        seqvars = {seqvars{:} fn{id}};
    end
end


fn = fn(setdiff([1:length(fn)],ids));
for nf = 1:length(fn)
    if iscell(AllTrials.(fn{nf}))
    else
        if length(AllTrials.(fn{nf})) < length(AllTrials.Start)
            fprintf('Forcing end values for %s\n',fn{nf});
            AllTrials.(fn{nf})(length(AllTrials.Start)) = 0;
        end
    end
end
state.tt = TimeMark(state.tt,'Checked Ends');
%non-zero values of falsestart should indicate the time gap
if isfield(AllTrials,'TrueStart') & ~isfield(AllTrials,'FalseStart')
    AllTrials.FalseStart = AllTrials.TrueStart;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 0)) = 1;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 1)) = 0;
end
if ~isfield(AllTrials,'Nf') || length(AllTrials.Nf) < length(AllTrials.Start)
    AllTrials.Nf(length(AllTrials.Start)) = NaN;
end

if isfield(AllTrials,'Phaseseq')
    for j = 1:length(AllTrials.Phaseseq)
        id = find(AllTrials.Phaseseq{j} == 0);
        AllTrials.Phaseseq{j}  = 1;
    end
end
if ~isfield(AllTrials,'Events') || length(AllTrials.Events) < length(AllTrials.Start)
    AllTrials.Events{length(AllTrials.Start)} = [];
end
%If there are a minority of trials with False Starts (usually missing the
%StimON channel, or a long delay to the matching FRAMESIGNAL on the serial line
%set the result to 2, so that these can be exluded if necessary.
% only do this if the result is already not 0 (otherwise Add in Badfix
% trials
if isfield(AllTrials,'FalseStart') && sum(AllTrials.FalseStart > 0) < length(AllTrials.FalseStart)/5
    AllTrials.Result(find(AllTrials.FalseStart > 0 & AllTrials.Result > 0)) = 2;
end
    %fn is now a list of fields that are set for each trial automatically.
nexpts = 1;
phasevals  = [0:360];
phasevals(2:4) = [pi pi/4 3*pi/4];
lphasevals = [0 pi pi 0]; %%see binoc.c SetRandomPhase();
rphasevals = [0 pi 0 pi];
for nx = 1:length(AllExpts)
    Idx.exptno = nexpts;
    state.tt = TimeMark(state.tt,sprintf('Expt %d',nx));
    clear Trials;
    
    frpt = 1;
    a = AllExpts(nx).firsttrial:AllExpts(nx).lasttrial;
    if isfield(AllTrials,'sM') && median(AllTrials.sM(a)) == 26
        badtr = 1;
    else
        badtr = usebadtrials;
    end
    
    if badtr
        Header.usebadtrials = 1;
        igood = 1:length(a);
        fn = {fn{:} 'Result'};
    else
        igood = find(AllTrials.Result(a) > 0);
    end
    nt = length(igood);
    igood = a(igood);
    nu = 0; %number of ustim pulses
    if nt <= 3 && AllExpts(nx).result ~= CANCELEXPT && state.showerrs
        err = sprintf('%s Expt %d only %d good trials',Header.Name,nx,nt);
        Idx = AddError(err, Idx, 0);
    end
    if ~isfield(AllExpts,'result')
        AllExpts(1).result = 1;
    end
    if ~isfield(AllExpts,'e3')
        AllExpts(1).e3 = 'e0';
    end
    if nt> 3 & igood(1) < length(AllTrials.Start) & ismember(AllExpts(nx).result,[2 0])
       spkids = [];
       needfields = {};
       if isempty(AllExpts(nx).e3)
           fprintf('Empty Type Expt %d, trials %d - %d\n',nx, AllExpts(nx).firsttrial, AllExpts(nx).lasttrial)
       elseif sum(strcmp(AllExpts(nx).e3,'ar'))
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
        fastseq = 0;
        for nf = 1:length(fn)
            fsz= size(AllTrials.(fn{nf}));
            if iscell(AllTrials.(fn{nf}))
                if strcmp(fn{nf},'uStimt')
                else
                nv = unique({AllTrials.(fn{nf}){igood}});
                Expt.Stimvals.(fn{nf}) = AllTrials.(fn{nf}){a(1)};
                end
            else
                if fsz(1) > 1 && fsz(2) > 1
%for now only check unique of element 1. If 1s all same, and 2s all same, but 1 ~= 2, 
%don't need it
                    nv = unique([AllTrials.(fn{nf})(igood)]);
                else
                    nv = unique([AllTrials.(fn{nf})(igood)]);
                end
            end
            
            if isnumeric(nv) & sum(~isnan(nv))> 1
                if fsz(1) > 1 && fsz(2) > 1
                    for j = 1:nt
                        Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j),:);
                    end
                else
                    for j = 1:nt
                        Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j));
                    end
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
        if isfield(AllTrials,'dfx') && size(AllTrials.dfx,2) > 1
            for j = 1:nt
                id = find(AllTrials.dfx(igood(j),:) ~= 0);
                Trials(j).dfx = AllTrials.dfx(igood(j),id);
                Trials(j).dfy = AllTrials.dfy(igood(j),id);
            end
            
        end
        if isfield(AllTrials,'rf');
            Expt.Stimvals.rf = median(AllTrials.rf(igood,:));
            Expt.Stimvals.st = mode(AllTrials.st(igood));
        end
    duration = mean([Trials.End] - [Trials.Start]);
    et = Expt.Stimvals.et;
    e2 = Expt.Stimvals.e2;
    if ~isfield(Expt.Stimvals,'e3')
        Expt.Stimvals.e3 = 'e0';
    end
    if isfield(Expt.Stimvals,'Fr') && Expt.Stimvals.Fr > 0
        frpt = Expt.Stimvals.Fr;
    end
    if isfield(AllExpts(nx),'e1vals') & ~isempty(AllExpts(nx).e1vals)
        Expt.e1vals = AllExpts(nx).e1vals;
    end
    if isfield(AllExpts(nx),'xovals') & ~isempty(AllExpts(nx).xovals)
        Expt.xovals = AllExpts(nx).xovals;
    elseif isfield(Expt,'xovals')
        Expt = rmfield(Expt,'xovals');
    end
    if sum(strncmp('ce',{et e2},2))
        Expt.Stimvals.ce = median(abs(AllTrials.ce(igood)));
    end
     if sum(strcmp(et,{'Op' 'Pp'})) & isfield(Expt.Stimvals,'rf')
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
    if isfield(AllExpts(nx),'e2vals') && ~isempty(AllExpts(nx).e2vals)
        Expt.e2vals = AllExpts(nx).e2vals;
        if strcmp(Expt.Stimvals.e2,'ce') && max(Expt.e2vals) > 1
            bid = find(Expt.e2vals > 1);
            err = sprintf('Fixing ce values (%s->1)',sprintf('%.1f ',Expt.e2vals(bid)));
            Idx = AddError(Idx, err);
            Expt.e2vals(bid) = 1;
        end

    end
                if isfield(Expt.Stimvals,'Fs')
        xovals = [-16:16] .* Expt.Stimvals.Fs;
        yovals = [-16:16] .* Expt.Stimvals.Fs;
        end
    allev = [];
    serrid = [];
    nframes = [];
    psychtrial = 0;
    seqtrial = 0;
    crtrial = 0; %count # with contrast reversal
    timesexpt = 0;
    nu = 0; %count trials with uStimt;
    for k = 1:nt
        Idx.t = AllTrials.Start(igood(k));
        if ~isempty(AllTrials.Cluster{igood(k)})
        sid = find(ismember(AllTrials.Cluster{igood(k)}, thecluster));
        Trials(k).Spikes = round(AllTrials.Spikes{(igood(k))}(sid));
        Trials(k).count = sum(find(Trials(k).Spikes > 500 & Trials(k).Spikes < duration+500));
        else
            Trials(k).Spikes = [];
            Trials(k).count = 0;
        end
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
        if ~isempty(strfind(Trials(k).OptionCode,'+2a')) || ~isempty(strfind(Trials(k).OptionCode,'+afc'))
            psychtrial = psychtrial+1;
        end
        if isfield(AllTrials,'Dc')
            Dcval = AllTrials.Dc(igood(k));
        else
            Dcval = 0;
        end
        if strfind(Trials(k).OptionCode,'+fS') & Dcval < 1
            seqtrial = seqtrial+1;
            fastseq = 1;
        else
            fastseq = 0;
        end
        if strfind(Trials(k).OptionCode,'+cr')
            crtrial = crtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+x2')
            timesexpt = timesexpt+1;
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
        %If just changing im seed slowly, don't mess wiht end/start
        if frpt > 1 && fastseq == 0
            frpt = 1;
        end
% for image seqs, Stimseq records the order of seeds, and there will not
% be a conversion to stimulus type
%NB. Seedseq is NOT in seqvars, and is handled separately. Don't want to
%build substpace maps, to don't make Start/End into vectors
       needseq = 0;
       for ns = 1:length(seqvars)
           if length(AllTrials.(seqvars{ns}){igood(k)}) > 1
               Trials(k).(seqvars{ns}) = AllTrials.(seqvars{ns}){igood(k)};
               needseq = length(Trials(k).(seqvars{ns}));
           end
       end
        if isfield(AllTrials,'Seedseq') && length(AllTrials.Seedseq{igood(k)}) > 1
            Trials(k).Seedseq = AllTrials.Seedseq{igood(k)};
        end
        for nf = 1:length(cellfields) %fileds that are cell arrays
            f = cellfields{nf};
            if isfield(AllTrials,f) && length(AllTrials.(f)) >= k
                Trials(k).(f) = AllTrials.(f){k};
            end
        end
        for f = 1:length(needfields)
            Trials(k).(needfields{f}) = AllTrials.(needfields{f})(igood(k));
        end
        if isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 1 && ...
                ((isfield(Expt,'e1vals') && ~isempty(Expt.e1vals)) || isfield(Expt,'xovals')) 
            if AllTrials.Stimseq{igood(k)}(end) > length(Expt.e1vals) %last value sometimes junk               
               AllTrials.Stimseq{igood(k)} = AllTrials.Stimseq{igood(k)}(1:end-1);
            end
            if length(Expt.e1vals) && ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ((max(AllTrials.Stimseq{igood(k)}) < length(Expt.e1vals) && seqtrial > k/2) || strcmp(et,'backMov'))
            evid = AllTrials.Stimseq{igood(k)}+1;
            if frpt > 1
                evid = evid(1:frpt:end);
            end
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
            if sum(strcmp(et,'Dc'))
                ev = Expt.e1vals(evid);
                Trials(k).(et) = AllTrials.Dc(igood(k));
                Trials(k).(Expt.Stimvals.e2)= ev;
                Trials(k).(e2)(find(ev == ISIGNALFRAME)) = AllTrials.(e2)(igood(k)); %% blanks
                if strcmp(Expt.Stimvals.e2,'or')
                    Trials(k).ori = AllTrials.(e2)(igood(k));
                end
            elseif sum(strcmp(et,'backMov'))
                ev = AllTrials.Stimseq{igood(k)}+1;
                Trials(k).(et) = ev;
             %   eb = Expt.e2vals(AllTrials.Stimseq{igood(k)}+1);
              %  Trials(k).(Expt.Stimvals.e2)= eb;
            elseif sum(strcmp(et,{'ic' 'pR'})) && Dcval > 0 && Dcval < 1
                ev = Expt.e1vals(evid);
                Trials(k).ori = AllTrials.or(igood(k));
                Trials(k).or = ev;
            else
                ev = Expt.e1vals(evid);
                Trials(k).(et) = ev;
                 if isfield(Expt,'e2vals') && length(Expt.e2vals) > max(AllTrials.Stimseq{igood(k)})
                    eb = Expt.e2vals(evid);
                    Trials(k).(Expt.Stimvals.e2)= eb;
                 else
                     eb = zeros(size(ev));
                end
            end
            if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
            end
            if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
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
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 1
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
                fastseq && ...
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
            if strcmp(Expt.Stimvals.et,'serange')
                if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                    Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
                end
                if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                    Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
                end
            end
            serrid = [serrid k];
            end
        elseif needseq
            Trials(k).Start = Trials(k).Start + [0:needseq-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;            
        end
        if isfield(AllTrials,'rptframes') && ~isempty(AllTrials.rptframes{igood(k)})
            Trials(k).rptframes = AllTrials.rptframes{igood(k)};
        end
        
        durs(k) = Trials(k).End(end) - Trials(k).Start(1);
        nframes(k) = length(Trials(k).Start);
    end
    if isfield(Trials,'ob') && sum([Trials.ob] < 0) && ~isfield(Trials,'or')
        aid = find([Trials.ob] >= 0);
        bid = find([Trials.ob] < 0);
        [Trials(aid).or] = deal(Expt.Stimvals.or);
        if Expt.Stimvals.or > 25
            [Trials(bid).or] = deal(Expt.Stimvals.or-90);
        else
            [Trials(bid).or] = deal(Expt.Stimvals.or+90);
        end
        for j = 1:length(Trials)
            Trials(j).ob = abs(Trials(j).ob);
        end
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

%some expts may have nfr > 2 on some trials, nut not other
%e.g. image tuning exps (seedseq ->nfr > 2) with interleaved blanks.
%only trim for missing sequence if its really an fastseq expt
    if seqtrial > nt/2 && mean(nframes) > 2
        fastseq = 1;
    else
        fastseq = 0;
    end
    if fastseq && state.alltrials == 0
       id = find(nframes > 1);
       Trials = Trials(id);
       bid = find(nframes  ==1);
       if length(bid) > 2
           Idx = AddError(Idx, 'Ex %d Removing %d Trials becuase no RC sequence\n',nexpts,length(bid));
       end
    end
    if nu == 0  && isfield(Trials,'uStimt')
        Trials = rmfield(Trials,'uStimt');
    end
    if isfield(Trials,'uStim') && sum([Trials.uStim]) == 0
            Trials = rmfield(Trials,'uStim');
    end

    if isfield(Trials,'inexpt') && sum([Trials.inexpt] ==0 > 0)
       bid = find([Trials.inexpt] == 0);
       id = find([Trials.inexpt] > 0);
       ids = unique([Trials(bid).id]);
       if ~isempty(bid)
           Idx.t = Trials(bid(1)).Start(1);
           if length(ids) <= 1
               Trials = Trials(id);
               Idx = AddError(Idx, 'Ex %d Removing %d Trials becuase not in Expt\n',nexpts,length(bid));
           else
               Idx = AddError(Idx, 'Ex %d Has  %d Trials not in Expt\n,',nexpts,length(bid));
           end
       end
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
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
% If reaal End of expt is long time after last trial, probabaly means user did
% EndExpt. Still want to record real time of endexpt
        if length(AllExpts) < nexpts || Expt.Trials(end).End(end)+10000 > AllExpts(nx).end
            AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
        else
            AllExpts(nx).end = AllExpts(nx).end+5000;
        end
    else %%online, unfinished, edpt
        if nx > 1  
            fprintf('Expt %d No end. Adding.\n');
        end
        AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
    end
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
        Expt.Header.trange(2) = AllExpts(nx).end;
    else
        Expt.Header.trange(2) = 0;
    end



    if isfield(Idx,'DigMark')
        t = Expt.Header.trange./10000;
        id = find(Idx.DigMark.times > t(1) & Idx.DigMark.times < t(2));
        if length(id)
            Expt.DigMark.times = Idx.DigMark.times(id);
            Expt.DigMark.codes = Idx.DigMark.codes(id,1);
        end
    end
    if psychtrial > nt/2
        Expt.Header.psych = 1;
    else
        Expt.Header.psych = 0;
    end
    %if Expt.Stimvals.x2 is set, it means it was manually set via AddTxt
    if ~isfield(Expt.Stimvals,'x2') || Expt.Stimvals.x2 == 0
        if timesexpt > nt/2
            Expt.Stimvals.x2 = 1;
        else
            Expt.Stimvals.x2 = 0;
        end
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
        Idx = AddError(Idx, err);
        id = find(durs > 0);
        Expt.Trials = Trials(id);
    end
    if ~isfield(Expt,'e1vals') & isfield(Expt.Trials,et)
        Expt.e1vals = unique(cat(1,Expt.Trials.(et)));
    end
 
    if ~strcmp('e0',Expt.Stimvals.e2) && strcmp('dx',Expt.Stimvals.et) && Expt.Stimvals.ve < 10.1027
        Expt.Stimvals.x2 = 1;
    end
    
%replace off screen positiosn with 0, but only for 
%regular tuning curves
    if sum(strncmp(Expt.Stimvals.et,{'xo' 'yo' 'Op' 'Pp'},2)) && Expt.Header.rc == 0
        exv = [Expt.Trials.(Expt.Stimvals.et)];
        if size(exv,1) == 1
        id = find(abs(exv) > 35); %off screen
        for j = 1:length(id)
            Expt.Trials(id(j)).xo = Expt.Stimvals.rf(1);
            Expt.Trials(id(j)).yo = Expt.Stimvals.rf(2);
            Expt.Trials(id(j)).Op = 0;
            Expt.Trials(id(j)).Pp = 0;
            Expt.Trials(id(j)).st = 0;
        end
        end
    end
    
    if strcmp(Expt.Stimvals.et,'fx') && strcmp(Expt.Stimvals.e2,'fy')
        emname = strrep(Expt.Header.Name,'.mat','.eyecal.mat');
        if ~exist(emname,'file')
            MakeEyeCalfile(Expt,emname);
        else
            a = load(emname);
            if isfield(a,'gains')
                fprintf('Eye gains RH %.2f LH %.2f RV %.2f LV %.2f\n',a.gains(1),a.gains(2),a.gains(3),a.gains(4));
            end
        end
    end
%     if strmatch(et,'Dc')
%         ors = cat(2,Expt.Trials.or);
%     end
    if isfield(Idx,'errexpt')
        id = find(Idx.errexpt == nexpts);
        if ~isempty(id)
            Expt.errs.msg = Idx.errs(id);
            if isfield(Idx,'errimtes')
                Expt.errs.t = Idx.errtimes(id);
            else
                Expt.errs.t = zeros(size(id));
            end
        end
    end
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
d = 0;
if length(did)
    for j = 1:length(did)
        ds = Text.text(did(j),:);
        dsid = strfind(ds,'Creat');
        if length(dsid)
            d = datenum(ds(dsid(1)+8:end));
            break;
        end
    end
end
if d == 0 %still not found
    did = strmatch('vebinoc',Text.text);
    if ~isempty(did)
        ds = Text.text(did(1),:);
        did = strfind(ds,' ');
         if length(did)
            d = datenum(ds(did(1)+1:end));
        end
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
ufl = strrep(name,'.mat','.ufl');
rid = strmatch('cm=rf',Text.text);


AddTxtFile = strrep(name,'.mat','Add.txt');
if ~exist(AddTxtFile)
    AddTxtFile= regexprep(name,'\.[0-9]*.mat','Add.txt');
end
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
    a = sscanf(Text.text(rid(j),:)','cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    rfs(j,1:length(a)) = a;
end
else
    for j = 1:length(rfstrs)
        a = sscanf(rfstrs{j},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(j,1:length(a)) = a;
    end
end
% find lines suggesting RF was changed after a quantitative measure
oid = strmatch('RO',Text.text);
pid = strmatch('RP',Text.text);
sid = [oid pid strmatch('RO',Text.text)];
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

function Idx = AddError(err, varargin)
%Idx = AddError(Idx, varargin)
% or Idx = AddError(err, Idx, show) old style
if ischar(err) %% old style
    Idx = varargin{1};
    show = varargin{2};
else  
    Idx = err;
    err = sprintf(varargin{:});
    
    if isfield(Idx,'state')
        show = Idx.state.showerrs;
    else
        show = 0;
    end
end
        
if ~isfield(Idx,'errs') | isempty(strmatch(err,Idx.errs)) 
    if show
        msgbox(err,'APlaySpkFile Error!!','modal');
    end
    mycprintf('errors',[err '\n']);
    Idx.errs = {Idx.errs{:} err};
    nerr = length(Idx.errs);
    if isfield(Idx,'t')
        Idx.errtimes(nerr) = Idx.t;
    end
    if isfield(Idx,'exptno')
        Idx.errexpt(nerr) = Idx.exptno;
    end
    Idx.newerrs = Idx.newerrs+1;
else
    mycprintf('errors',[err '\n']);
end

function FindMissingTimes(Events, Text, bstimes, estimes,bsid,fsid,ts,te)

bid = strmatch('bss',Text.text);
bsstimes = Text.times(bid);
tid = find(bsstimes >  ts & bsstimes < te);
eid = find(bstimes > ts & bstimes < te);
x = bsstimes-bstimes(1:length(bsstimes));
hold off;
plot(bsstimes,x,'o');
hold on;
x = bsstimes(tid)-bstimes(eid(1:length(tid)));
plot(bsstimes(tid),x,'ro');
plot([ts ts],get(gca,'ylim'));
plot([te te],get(gca,'ylim'));


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


