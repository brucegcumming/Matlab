function [avg, details] = SpTrigLFP(Trials, duration, samplerate, w, varargin)
%
%SpTrigLFP(Trials, duration, samplerate, w)
%Builds a spike triggered LFP for a set of trials
% By tdefault also builds a voltage histogram for the LFP, which can be
% slow
%SpTrigLFP(...,'nohist')  turns  this off;
%
%sample rate is samples per tic
% for spike2 LFP this is 1/(Expt.Header.LFPsamplerate * 10000)
%
% details.times = time values for the average (in tics)
latency = 500;
nspk = 0;
lfpch = 1;
histx = [];
shuffle = 0;
frameperiod = 166.7;
mkhist = 1;
nch = size(Trials(1).LFP,2);
j =1;
while j <= length(varargin)
    if isfield(varargin{j},'loadname')
        details = varargin{j};
    elseif strncmpi(varargin{j},'framer',5)
        j = j+1;
        frameperiod = varargin{j};
    elseif strncmpi(varargin{j},'histvals',5)
        j = j+1;
        histx = varargin{j};
    elseif strncmpi(varargin{j},'shuffle',5)
        shuffle = 1;
    elseif strncmpi(varargin{j},'nohist',5)
        mkhist = 0;
    end
    j = j+1;
end
slfpavg = zeros(2*w+1,nch);
lfpavg = zeros(2*w+1,nch);
lfpavn = zeros(2*w+1,1);
lfprange = round([latency  latency+duration] .* samplerate);
len = length(lfpavg);
splfpv = [];
allfpv = [];
lfpv = [];
toff = 0;
dtimes = {};
for tr = 1:length(Trials);
    lfp = Trials(tr).LFP;
    if shuffle
        if tr > 1
            slfp = Trials(tr-1).LFP;
        else
            slfp = Trials(end).LFP;
        end
        slfp(isnan(slfp)) = 0;
    end
    lfp(isnan(lfp)) = 0;
    lfppre = Trials(tr).Start(1)-Trials(tr).lfptime;
    if ~isempty(lfp)
        if isfield(Trials,'spkuse')
            spk = Trials(tr).spkuse;
            dtimes{tr} = 1+round((Trials(tr).stdelays+lfppre)* samplerate);
        else
            spk = find(Trials(tr).Spikes > latency & Trials(tr).Spikes < latency+duration);
        end
    nspk = nspk + length(spk);
%spiketimes is in tics. tics/(tics per sample) = samples
    
    lfptimes = 1+round(Trials(tr).Spikes(spk) * samplerate);
    alltimes = [lfprange(1):lfprange(2)];
    if isfield(Trials,'lfpo')
        lfptimes = lfptimes+Trials(tr).lfpo;
        alltimes = alltimes+Trials(tr).lfpo;
        alltimes = alltimes(find(alltimes < length(Trials(tr).LFP)));
    end
% use start time of LFP in tics if available, to  reduce rounding errors
% N.B if the time of a spike was  exactly at lfptime, this would mean 
% suing lfp sample 1, not 0, so we have to add 1 to the calculation. 
    if isfield(Trials,'lfptime') & Trials(tr).lfptime
        lfptimes = 1+round((Trials(tr).Spikes(spk) + Trials(tr).Start(1)-Trials(tr).lfptime)* samplerate);
    end
    sptimes{tr} = lfptimes;
    
    if mkhist
    allfpv = [allfpv Trials(tr).LFP(alltimes,:)'];
    end
%    lfptimes = round(Trials(tr).Spikes(spk) ./(samplerate * 10000));
    for lfpt = lfptimes'
        if lfpt <= w
            start = 1;
            last = lfpt+w;
            if last > size(lfp,1)
                last = size(lfp,1);
            end
            istart = length(lfpavg)- last+1;
            iend = length(lfpavg);
        elseif lfpt >= 1+length(lfp)-w
            last = length(lfp);
            istart = 1;
            start = lfpt -w;
            iend = last - start+1;
        else
            start = lfpt-w;
            last = lfpt+w;
            istart = 1;
            iend = len;
        end
        if istart > w | iend < w
            err = iend-istart;
        end
        lfpavg(istart:iend,:) = lfpavg(istart:iend,:) + lfp(start:last,:);
        if shuffle
            slfpavg(istart:iend,:) = slfpavg(istart:iend,:) + slfp(start:last,:);
        end
        lfpavn(istart:iend) = lfpavn(istart:iend)+1;
    end
%        lfpv = [lfpv lfp(lfptimes)'];
    end
end
idx = find(lfpavn == 0);
lfpavg(idx,:) = 0;
lfpavn(idx) = -1;
if shuffle
    details.shuffle = slfpavg ./ repmat(lfpavn,1,nch);;
end
avg = lfpavg ./ repmat(lfpavn,1,nch);
lfp = cat(3,Trials.LFP);
stimlfp = mean(lfp,3);
lfpid = [lfprange(1):lfprange(end)];

if isfield(Trials,'lfpo')
lfpo = round(mean([Trials.lfpo]));
else
    lfpo = 0;
end
id = find(lfpid+lfpo > 0 & lfpid++lfpo <= size([Trials.LFP],1));
lfpid = lfpid(id);
spktimes = round(lfpid./samplerate);
if length(spktimes)
spksdf = trigsdfa(Trials,100,spktimes);
stimcorr = stimlfp(lfpid+lfpo).*spksdf';
else
    spksdf = zeros(size(spktimes));
end
if mkhist
if isempty(histx)
[y,x] = smhist(allfpv,'sd',0.1);
else
[y,x] = smhist(allfpv,'sd',0.1,'xval',histx);
end    
for toff = -4:20
    splfpv = [];
    for tr = 1:length(Trials)
        splfpv = [splfpv Trials(tr).LFP(sptimes{tr}+toff)'];
    end
    meanv(toff+5) = mean(splfpv);
end
end

% if just looking at certain stims in a subspace map
% then compare lfp at peak delay with spike-triggered lfp
% at +- 0.5 frames of that time
if ~isempty(dtimes) 
    splfpv = [];
    tplfpv = []; %lfp V at frametime
    stplfpv = []; %lfp V at frametime, if spike
    nstplfpv = []; %lfp V at frametime, if no spike
    prespv = [];
    sdtimes = {};
    for tr = 1:length(Trials)
        tplfpv = [tplfpv Trials(tr).LFP(dtimes{tr})'];
        sdtimes{tr} = [];
        for k = 1:length(dtimes{tr})
            dt = sptimes{tr} - dtimes{tr}(k);
            if min(abs(dt)) < (frameperiod/2) .* samplerate
                sdtimes{tr} = [sdtimes{tr} dtimes{tr}(k)];
            end
            counts{tr}(k) = sum(abs(dt) < (frameperiod/2) .*samplerate);
        end
        nsdtimes{tr} = setdiff(dtimes{tr},sdtimes{tr});
        ssdtimes{tr} = [];
        for k = 1:length(sptimes{tr})
            dt = sptimes{tr}(k) - dtimes{tr};
                 [a, b] = min(abs(dt));
                 ssdtimes{tr}(k) = dtimes{tr}(b);
                 ssdiffs{tr}(k) = dt(b);
        end
        splfpv = [splfpv Trials(tr).LFP(sptimes{tr})']; %V at time of spikes
        prespv = [prespv Trials(tr).LFP(ssdtimes{tr})']; %V at stim t nearest to spike
        stplfpv = [stplfpv Trials(tr).LFP(sdtimes{tr})']; %V for frames followed by spike
        nstplfpv = [nstplfpv Trials(tr).LFP(nsdtimes{tr})']; % V for frames followed by no spike
    end
    hold off;
    cmphist(splfpv,tplfpv);
    legend('spike','stim')
    allcounts = cell2mat(counts);
    for j = 0:max(allcounts)
        mv(j+1) = mean(tplfpv(find(allcounts == j)));
    end
    [a,b] = sort(tplfpv);
    v = smooth(a,100);
    p = smooth(allcounts(b),100);
else
    mv = 0;
    allcounts = [];
    tplfpv = [];
end



if mkhist
[maxm, toff]  = max(meanv);
splfpv = [];
for tr = 1:length(Trials)
    splfpv = [splfpv Trials(tr).LFP(sptimes{tr}+toff)'];
end
toff = toff -5;
if length(splfpv)
[ay, ax] = smhist(splfpv,'sd',0.1,'xval',x);
else
    ay = 0;
    ax = 0;
end
details.trigvy = ay;
details.vx = x;
details.vy = y;
details.toff = toff;
details.meanv = meanv(5);
end
details.nspk = nspk;
details.lfpavn = lfpavn;
details.times = [-w:w] ./ samplerate; %in tics
details.mv = mv;
details.counts = allcounts;
details.framev = tplfpv;
