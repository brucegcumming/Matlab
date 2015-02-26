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
    else
        return;
    end
elseif ~loaded
    if Expt.Header.bysuffix  && (isfield(Expt.Header,'fileprefix') || isfield(Expt.Header,'loadname'))
        if isfield(Expt.Header,'Combineids')
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
                if ~isempty(breaks);
                    LFP = AddError(LFP,'%s Ids go backwards at %d->%d',LFP.Header.Name,ids(breaks(1)),ids(breaks(1)+1));
                end
                if ~isfield(LFP,'errs')
                    LFP.errs = [];
                end
                LFPS(j) = LFP;
                Names{j} = lfpfile;
            else
                fprintf('Missing LFP File %s\n',lfpfile);
            end
        end
        LFP.Trials = cat(1,[LFPS.Trials]);
        LFP.errs = cat(1,LFPS.errs);
        LFP.Header.Names = Names;
        elseif isfield(Expt.Header,'suffix')
            if isfield(Expt.Header,'fileprefix')
                [a,b] = fileparts(Expt.Header.fileprefix);
            else
                [a,b] = fileparts(Expt.Header.loadname);
                b = regexprep(b,'\.[0-9]*','');
            end
            lfpfile = [fileparts(Expt.Header.loadname) '/' b 'A.' num2str(Expt.Header.suffix) '.lfp.mat'];
            load(lfpfile);
            LFP.Header.Names{1} = lfpfile;
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
                fprintf('No LFP datafile: %s or %s\n',lfpfile,a);
                return;
            end
        end
        load(lfpfile);
    end
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
    endtimes(j) = LFP.Trials(j).ftime + size(LFP.Trials(j).LFP,1) .* LFP.Header.CRsamplerate .* 10000;
end

for j = 1:length(Expt.Trials)
    id = find([LFP.Trials.ftime] < Expt.Trials(j).Start(1) & endtimes > Expt.Trials(j).End(end));
    if length(id) == 1
        trial(j) = id;
        starts(j) = Expt.Trials(j).Start(1);
        diffs(j) = 0; %if one and only one  LFP trial spans the Expt Trial, it must be the one
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

if isfield(LFP,'errs') && ~isempty(LFP.errs)
    Expt.LFPerrs = LFP.errs;
end

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
    elseif abs(diffs(j)) < 700 && endtimes(k) > Expt.Trials(j).End(end)
        Expt.Trials(j).LFP = LFP.Trials(k).LFP;
        Expt.Trials(j).lfpo = round(lfpos(j)); 
        Expt.Trials(j).lfptime = LFP.Trials(k).ftime;
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
        fprintf('%s No LFP data for Trial id %d:%.2f\n',lfpfile,Expt.Trials(j).id,Expt.Trials(j).Start(1)./10000);
        Expt.Trials(j).LFP = zeros(size(powerspec));
        lfpos(j) = (LFP.Trials(k).Start - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000);
        Expt.Trials(j).lfpo = 0;
        Expt.Trials(j).lfptime = 0;
        goodtrial(j) = 0;
    end
    tlfp = Expt.Trials(j).LFP;
    lfpft = fft(tlfp);
    Expt.Trials(j).FTlfp = lfpft;
    if size(lfpft,2) == lfpchans && goodtrial(j) && size(lfpft,1) == length(powerspec)
    powerspec = powerspec + abs(lfpft);
    end
end
Expt.Header.LFPsnr = 0;
Expt.Header.LFPtimes = ([1:lfplen] - mean(lfpos)) .* (LFP.Header.CRsamplerate);
Expt.Header.LFPgoodtrials = goodtrial;
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
