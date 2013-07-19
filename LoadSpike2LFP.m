function [Expt, details] = LoadSpike2LFP(Expt, varargin)

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
    if Expt.Header.bysuffix
        for j = 1:length(Expt.Header.Combineids)
            lfpfile = [Expt.Header.fileprefix 'A.' num2str(Expt.Header.Combineids(j)) '.lfp.mat'];
            if exist(lfpfile)
                load(lfpfile);
                LFPS(j) = LFP;
                Names{j} = lfpfile;
            else
                fprintf('Missing LFP File %s\n',lfpfile);
            end
        end
        LFP.Trials = cat(1,[LFPS.Trials]);
        LFP.Header.Names = Names;
    else
        if isempty(lfpfile)
            lfpfile = [drive strrep(Expt.Header.Name,'.mat','.lfp.mat')];
        end
        if ~exist(lfpfile,'file')
            lfpfile = name2path(splitpath(Expt.Header.Name),'smr');
            lfpfile = strrep(lfpfile,'.mat','.lfp.mat');
            if ~exist(lfpfile,'file')
                fprintf('No LFP datafile: %s\n',lfpfile);
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
end
for j = 1:length(Expt.Trials)
    if isfield(Expt.Trials,'TrialStart')
            [diffs(j), trial(j)] = min(abs(Expt.Trials(j).TrialStart-[LFP.Trials.Start]));
            starts(j) = Expt.Trials(j).TrialStart;
    else
            [diffs(j), trial(j)] = min(abs(Expt.Trials(j).Start(1)-[LFP.Trials.Start]));
            starts(j) = Expt.Trials(j).Start(1);
    end
    lfplen(j) = length(LFP.Trials(trial(j)).LFP);
end
details.ltrials = trial;
lfpchans = size(LFP.Trials(1).LFP,2);

%Select Trials with teh same data lenth
%If minimum is close to median, use minimum
%If fixcit is set, throw away fixcrit % of short trials
%Otherwise set them all to the shortest
lfplens = lfplen;
gid = find(lfplen > 0 & diffs < 700);
if min(lfplens)+5 > prctile(lfplen(gid),50)
    lfplen = min(lfplen(gid));
elseif fixcrit > 0
 lfplen = prctile(lfplen(gid),fixcrit);
 lfplen = floor(lfplen);
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


for j = 1:length(Expt.Trials)
    k = trial(j);
    if abs(LFP.Trials(k).Start - starts(j)) < 700 & size(LFP.Trials(k).LFP,1) >= lfplen
        Expt.Trials(j).LFP = LFP.Trials(k).LFP(1:lfplen,:);
        lfpos(j) = (Expt.Trials(j).Start(1) - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000);
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
    else
        Expt.Trials(j).LFP = zeros(size(powerspec));
        lfpos(j) = (LFP.Trials(k).Start - LFP.Trials(k).ftime)/(LFP.Header.CRsamplerate.*10000);
        Expt.Trials(j).lfpo = 0;
        Expt.Trials(j).lfptime = 0;
    end
    lfpft = fft(Expt.Trials(j).LFP);
    Expt.Trials(j).FTlfp = lfpft;
    if size(lfpft,2) == lfpchans
    powerspec = powerspec + abs(lfpft);
    end
end
Expt.Header.LFPsnr = 0;
Expt.Header.LFPtimes = ([1:lfplen] - mean(lfpos)) .* (LFP.Header.CRsamplerate);
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
