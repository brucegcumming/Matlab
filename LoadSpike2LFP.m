function [Expt, details] = LoadSpike2LFP(Expt, varargin)
%[Expt, details] = LoadSpike2LFP(Expt, varargin)
%Load LFP data from spike2 .lfp.mat file
%Also calculates power spectrum

reload = 0;
loaded = 0;
setid = 0;
fixcrit = 5;
drive = '';
details = [];
lfpfile = [];
autoscale = 0;
verbose = 0;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j},'Trials') & isfield(varargin{j}.Trials,'LFP') %% its an lfp.mat file
            LFP = varargin{j};
            loaded = 1;
        end
    elseif strncmpi(varargin{j},'autoscale',6)
        autoscale = 1;
    elseif strncmpi(varargin{j},'lfpfile',6)
        j = j+1;
       lfpfile = varargin{j};
    elseif strncmpi(varargin{j},'drive',5)
        j = j+1;
       drive = varargin{j};
    elseif strncmpi(varargin{j},'fixshort',4)
        j = j+1;
       fixcrit = varargin{j};
    elseif strncmpi(varargin{j},'reload',3)
       reload = 1;
    elseif strncmpi(varargin{j},'verbose',4)
        verbose =1;
   end
   j = j+1;
end

if isfield(Expt.Trials,'LFP') & ~reload & ~loaded
    return;
end
expname = GetEval(Expt,'name');
if ~isfield(Expt.Header,'bysuffix')
    Expt.Header.bysuffix = 0;
end

if iscell(lfpfile) 
    if exist(lfpfile{1},'file')
    load(lfpfile{1});

    nLFP = LFP;
    for j = 2:length(lfpfile)
        if exist(lfpfile{j},'file')
        load(lfpfile{j});
        nLFP.Trials = [nLFP.Trials LFP.Trials];
        end
    end
    LFP = nLFP;
    LFP.Header.Name = lfpfile{1};
    lfpfile = lfpfile{1};
    else
        return;
    end
elseif ~loaded
    if Expt.Header.bysuffix  && (isfield(Expt.Header,'fileprefix') || isfield(Expt.Header,'loadname'))
        if ~isfield(Expt.Header,'suffix')
            Expt.Header.suffix = GetExptNumber(Expt);
        end
        if isfield(Expt.Header,'Combineids')
        firsttrial = 1;
        details.gotlfp = 0;
        details.missinglfp = {};
        
        for j = 1:length(Expt.Header.Combineids)
            if isfield(Expt.Header,'suffixes') && length(Expt.Header.suffixes) == length(Expt.Header.Combineids)
                eid = Expt.Header.suffixes(j);
            else
                eid = Expt.Header.Combineids(j);
            end
            if isfield(Expt.Header,'loadname')
                [a,b] = fileparts(Expt.Header.fileprefix);
                lfpfile = [fileparts(Expt.Header.loadname) '/' b 'A.' num2str(eid) '.lfp.mat'];
            else
                lfpfile = [Expt.Header.fileprefix 'A.' num2str(eid) '.lfp.mat'];
            end
            if exist(lfpfile)
                load(lfpfile);
                ids = [LFP.Trials.id];
                breaks = find(diff(ids) <0);
                LFP.program = 'LoadSpike2LFP';
                if ~isempty(breaks);
                    if ~isfield(LFP.Header,'s2ver')
                        LFP = AddError(LFP,'%s Missing some ids. Error at id=%d/%d ReRun MakeMat in Spike2',LFP.Header.Name,ids(breaks(1)+1));
                    else
                        LFP = AddError(LFP,'%s Ids go backwards at %d->%d',LFP.Header.Name,ids(breaks(1)),ids(breaks(1)+1));
                    end
                end
                if isfield(LFP.Trials,'err')
                    for t = 1:length(LFP.Trials)
                        if ~isempty(LFP.Trials(t).err)
                            if ismember(LFP.Trials(t).id,[Expt.Trials.id])
                                str = 'Data Needed';
                            else
                                str = 'Data Not Needed';
                            end
                            LFP = AddError(LFP,'%s At id = %d: %s%s',LFP.Header.Name,LFP.Trials(t).id,LFP.Trials(t).err,str);
                            
                        end
                    end
                    LFP.Trials = rmfield(LFP.Trials,'err');
                end
                if ~isfield(LFP,'errs')
                    LFP.errs = [];
                    LFP.errdata = [];
                end
                TrialRange(j,1)  = firsttrial;
                LFPS(j) = LFP;
                Names{j} = lfpfile;
                idrange(j,:) = minmax([LFP.Trials.id]);
                details.gotlfp = details.gotlfp+1;
            else
                fprintf('Missing LFP File %s\n',lfpfile);
                details.missinglfp = {details.missinglfp{:} lfpfile};
            end
        end
        LFP.Trials = cat(1,[LFPS.Trials]);
        LFP.errs = cat(2,LFPS.errs);
        LFP.Header.Names = Names;
        LFP.Header.idrange = idrange;
        elseif isfield(Expt.Header,'suffix')
            if isfield(Expt.Header,'fileprefix')
                [a,b] = fileparts(Expt.Header.fileprefix);
            else
                [a,b] = fileparts(Expt.Header.loadname);
                b = regexprep(b,'\.[0-9]*','');
            end
            lfpfile = [fileparts(Expt.Header.loadname) '/' b 'A.' num2str(Expt.Header.suffix) '.lfp.mat'];
            if ~exist(lfpfile)
                fprintf('Missing LFP File %s\n',lfpfile);
                details.gotlfp = 0;
            else
                details.gotlfp = 1;
                load(lfpfile);
                LFP.Header.Names{1} = lfpfile;
            end
        end
    else
        if Expt.Header.bysuffix
            lfpfile = regexprep(expname,'([0-9])(\.[0-9]*)\.mat','$1A$2.lfp.mat');
        end
        if isempty(lfpfile)
            if strncmp(drive,expname,length(drive)) %already has drive in name
                lfpfile = [strrep(expname,'.mat','.lfp.mat')];
            else
                lfpfile = [drive strrep(expname,'.mat','.lfp.mat')];
            end
            if ~exist(lfpfile,'file')
                lfpfile =  strrep(lfpfile,'.lfp.mat','A.lfp.mat');
            end
            
        end
        if ~exist(lfpfile,'file')
            a = lfpfile;
            lfpfile = name2path(splitpath(expname),'smr');
            lfpfile = strrep(lfpfile,'.mat','.lfp.mat');
            if ~exist(lfpfile,'file')
                details.gotlfp = 0;
                fprintf('No LFP datafile: %s or %s\n',lfpfile,a);
                return;
            end
        end
        load(lfpfile);
        details.gotlfp = 1;
        
    end
end

if ~exist('LFP','var') %no data
    return;
end
if ~isfield(LFP.Trials,'LFP')
    fprintf('No LFP Data in %s\n',lfpfile);
    return;
end
nmissing = 0;
clipped = [];
missing = [];
if isempty(LFP.errs)
    LFP.errs = {};
end
    if ~isfield(LFP.Header,'amps')
        LFP.Header.amps = LFPGains(LFP);
    end
for j = 1:length(LFP.Trials)
    if isempty(LFP.Trials(j).Start)
        LFP.Trials(j).Start = NaN;
    end
    if ~isfield(LFP.Trials,'id') | isempty(LFP.Trials(j).id)
        LFP.Trials(j).id = NaN;
    end
    if isempty(LFP.Trials(j).ftime)
        if ismember(LFP.Trials(j).id,[Expt.Trials.id])
            str = 'Data Needed';
            nmissing = nmissing+1;
        else
            str = 'Data Not Needed';
        end
        AddError(LFP,'Missing ftime in Trials%d,Id%d %s',j,LFP.Trials(j).id,str);
        endtimes(j) = NaN;
        LFP.Trials(j).ftime = NaN;
    else
        endtimes(j) = LFP.Trials(j).ftime + size(LFP.Trials(j).LFP,1) .* LFP.Header.CRsamplerate .* 10000;
    end
end

lfpid = [LFP.Trials.id];
Tid = [Expt.Trials.id];
eid = GetExptNumber(Expt);
for j = 1:length(Expt.Trials)
    T = Expt.Trials(j);
    id = find([LFP.Trials.ftime] < Expt.Trials(j).Start(1) & endtimes > Expt.Trials(j).End(end));
    idid = find(lfpid == Expt.Trials(j).id);
    if length(id) == 1
        trial(j) = id;
        lT =  LFP.Trials(id);
        starts(j) = Expt.Trials(j).Start(1);
        diffs(j) = 0; %if one and only one  LFP trial spans the Expt Trial, it must be the one
    elseif length(idid) ==1 % ? && ismember(idid,id)?
        LT =  LFP.Trials(idid);
        d = (Expt.Trials(j).Start(1)-LFP.Trials(idid).Start);
        tmiss(1) = (T.Start(1) - LT.ftime)/10; %ms 
        tmiss(2) = (endtimes(idid)-T.End(end))/10; %ms 
        if isempty(id)
            if verbose || tmiss(1) < 0 || tmiss(2) < 0
                LFP = AddError(LFP,'E%d LFP Trial %d %.2f -> %.2f But Trial %d id%d %.2f->%.2f. have %.1fms before StimOn,%.1fms after Off',eid,idid,...
                LT.ftime./10000,endtimes(idid)./10000,j,T.id,T.Start(1)./10000,T.End(end)./10000,tmiss(1),tmiss(2));
                clipped = [clipped j];
            end
        elseif ismember(lfpid(idid),Tid(j)) && abs(d) < 700
            if abs(d) > 100 %< 10ms is just because exact frame time not in LFP Files
                fprintf('Using LFP Trial %d for Trials %d at %.2f. Start diff is %.1fms\n',idid,j,T.Start(1)./10000,d/10);
            end
        else
           LFP = AddError(LFP,'No LFP Data for Trial %d Id%d at %.2f: mismatched Id',j,T.id,T.Start(1)./10000);           
        end
        diffs(j) = abs(d);
        trial(j) = idid;
        starts(j) = Expt.Trials(j).Start(1);        
    elseif isfield(Expt.Trials,'Start')
        [diffs(j), trial(j)] = min(abs(Expt.Trials(j).Start(1)-[LFP.Trials.Start]));
        starts(j) = Expt.Trials(j).Start(1);
    elseif isfield(Expt.Trials,'TrialStart')
            [diffs(j), trial(j)] = min(abs(Expt.Trials(j).TrialStart-[LFP.Trials.Start]));
            starts(j) = Expt.Trials(j).TrialStart;
    end
    xdiffs(j) = Expt.Trials(j).Start(1)-LFP.Trials(trial(j)).Start;
    lfplen(j) = length(LFP.Trials(trial(j)).LFP);
    lfpos(j) = round((Expt.Trials(j).Start(1) - LFP.Trials(trial(j)).ftime)/(LFP.Header.CRsamplerate.*10000));

    if isfield(Expt.Trials,'Result')
        goodtrial(j) = Expt.Trials(j).Result;
    else
        goodtrial(j) = 1;
    end
end
details.ltrials = trial;
lfpchans = size(LFP.Trials(1).LFP,2);

%Select Trials with teh same data lenth
%If minimum is close to median, use minimum
%If fixcit is set, throw away fixcrit % of short trials
%Otherwise set them all to the shortest
lfplens = lfplen;
gid = find(lfplen > 0 & diffs < 700 & goodtrial);
if min(lfplens)+5 > prctile(lfplen(gid),50)
    lfplen = min(lfplen(gid));
elseif fixcrit > 0
 lfplen = prctile(lfplen(gid),fixcrit);
 lfplen = floor(lfplen-1);
else
 lfplen = min(lfplen(gid));
end

if isempty(lfplen) || isnan(lfplen)
    fprintf('No matching Trials in %s for %s\n',Expt.Header.Name,LFP.Header.Name);
    return;
end
Expt.Header.lfplen = lfplen;
%CRsamplerate, is from BinSiz() in Spike2, so is in fact the time distance
%between sample points
Expt.Header.LFPsamplerate = LFP.Header.CRsamplerate;
Expt.Header.LFPamps = LFP.Header.amps;

powerspec = zeros(lfplen,lfpchans);
if isfield(LFP.Header,'needscale') && LFP.Header.needscale > 0
    autoscale = 1;
end

if autoscale
        Expt.Header.LFPscaled = 1;
end

Expt = CopyErrs(Expt,LFP);

presamples = floor(prctile(lfpos,90));
Expt.Header.LFPpreperiod = presamples .* LFP.Header.CRsamplerate .* 10000;
for j = 1:length(Expt.Trials)
    k = trial(j);
    lfpos(j) = round((Expt.Trials(j).Start(1) - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000));
%don't worry about size of LFP records. This will be fixed later in FixLFPTrials    
    if abs(diffs(j)) < 700 & size(LFP.Trials(k).LFP,1) >= lfplen
        Expt.Trials(j).LFP = LFP.Trials(k).LFP;
        Expt.Trials(j).lfpo = round(lfpos(j)); 
        Expt.Trials(j).lfptime = LFP.Trials(k).ftime;
        if setid
            LFP.Trials(k).id = Expt.Trials(j).id;
        end
        if autoscale 
            for c = 1:size(Expt.Trials(j).LFP,2)
                Expt.Trials(j).LFP(:,c) = Expt.Trials(j).LFP(:,c)./ LFP.Header.amps(c);  
            end
        end
        goodtrial(j) = 1;
    elseif abs(diffs(j)) < 700 && endtimes(k) > Expt.Trials(j).End(end) %Missing start
        Expt.Trials(j).LFP = LFP.Trials(k).LFP;
        Expt.Trials(j).lfpo = round(lfpos(j)); 
        Expt.Trials(j).lfptime = LFP.Trials(k).ftime;
        goodtrial(j) = 3;
    elseif goodtrial(j) == 0
        lfpos(j) = (Expt.Trials(j).Start(1) - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000);
        Expt.Trials(j).lfpo = round(lfpos(j)); 
        Expt.Trials(j).lfptime = LFP.Trials(k).ftime;
        Expt.Trials(j).LFP = LFP.Trials(k).LFP;
        if autoscale 
            for c = 1:size(Expt.Trials(j).LFP,2)
                Expt.Trials(j).LFP(:,c) = Expt.Trials(j).LFP(:,c)./ LFP.Header.amps(c);  
            end
        end
        goodtrial(j) = 2;
    else
        if isfield(LFP.Header,'Names') && length(LFP.Header.Names) > 1 %mulitple files
            b = FindBlock(LFP, Expt.Trials(j).id);
            LFP = AddError(LFP,'%s No LFP data for Trial %d id %d:%.2f',LFP.Header.Names{b},j,Expt.Trials(j).id,Expt.Trials(j).Start(1)./10000);
        else
            LFP = AddError(LFP,'%s No LFP data for Trial %d id %d:%.2f',lfpfile,j,Expt.Trials(j).id,Expt.Trials(j).Start(1)./10000);
        end
        Expt.Trials(j).LFP = zeros(size(powerspec));
        lfpos(j) = (LFP.Trials(k).Start - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000);
        Expt.Trials(j).lfpo = 0;
        Expt.Trials(j).lfptime = 0;
        goodtrial(j) = 0;
        missing = [missing j];
        nmissing= nmissing+1;
    end
    tlfp = Expt.Trials(j).LFP;
%make sure FTlpf has same size in all trials    
    lfpft = fft(tlfp);
    if size(lfpft,1) >= lfplen
        Expt.Trials(j).FTlfp = lfpft(1:lfplen,:);
    else
        Expt.Trials(j).FTlfp = lfpft;
        Expt.Trials(j).FTlfp(end+1:lfplen,:) = NaN;
    end
    if size(lfpft,2) == lfpchans && goodtrial(j) && size(lfpft,1) == length(powerspec)
    powerspec = powerspec + abs(lfpft);
    end
end
Expt.Header.LFPsnr = 0;
Expt.Header.LFPtimes = ([1:lfplen] - mean(lfpos)) .* (LFP.Header.CRsamplerate);
Expt.Header.LFPgoodtrials = goodtrial;
Expt.Header.LFPmissing = length(missing);
Expt.Header.LFPclipped = length(setdiff(clipped,missing));
Expt = CopyErrs(Expt,LFP);
if isfield(Expt.Header,'trange')
    trange = Expt.Header.trange;
    trange(2) = trange(2) + 1000;
else
    trange(1) = Expt.Trials(1).Start(1) - 10000;
    trange(2) = Expt.Trials(end).End(end) + 10000;
end
    
if isfield(LFP,'Pulses')
    cid = find([LFP.Pulses.ftime] > trange(1) & [LFP.Pulses.ftime] < trange(2));
    if ~isempty(cid)
    Expt.Pulses = LFP.Pulses(cid);
    end
end

function Expt = CopyErrs(Expt,LFP)

if isfield(LFP,'errs') && ~isempty(LFP.errs)
    Expt.LFPerrs = LFP.errs;
end

if isfield(LFP,'errs')
    Expt.LFPerrs = LFP.errs;
    if isfield(LFP,'errdata')
        Expt.LFPerrdata = LFP.errdata;
    else
        for j = 1:length(LFP.errs)
            Expt.LFPerrdata(j).time = LFP.Header.ConvDate;
        end
    end
end

function block = FindBlock(LFP,id)

a = find(id >= LFP.Header.idrange(:,1) & id <= LFP.Header.idrange(:,2)) ;
if isempty(a)
    a = find(id >= LFP.Header.idrange(:,1));
    if isempty(a)
        a = 1;
    else
        a = a(1);
    end
end
block = a;

