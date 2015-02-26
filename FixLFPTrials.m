function Expt = FixLFPTrials(Expt, varargin)
%Expt = FixLFPTrials(Expt)
%Makes sure that All Trials.LFP fields have the same length, padding
%front and back with NaNs;
%Expt = FixLFPTrials(Expt,'double') converts int to doubles
%Expt = FixLFPTrials(Expt,'int') converts double back to int (do this
%before saving

MAXPROBES = 96;
pad = NaN;
needft = 0;
gotft = 0;
fixdouble = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',6)
        pad = 0;
    elseif strncmpi(varargin{j},'double',5)
        fixdouble = 1;
    elseif strncmpi(varargin{j},'int',3)
        fixdouble = 2;
    elseif strncmpi(varargin{j},'ft',2)
        needft = 1;
    elseif strncmpi(varargin{j},'fixspike',6)
        Expt = FixLFPSpike(Expt,varargin{:});
        return;
    end
    j = j+1;
end

if isfield(Expt,'Spikes') && iscell(Expt.Spikes) && isfield(Expt,'Expt') %AllExpt
    Expt.Expt = FixLFPTrials(Expt.Expt,varargin{:});
    return;
end


if ~isfield(Expt.Trials,'LFP') %load faileed
    return;
end

if fixdouble ==2 %convert to int
    if isfield(Expt.Header,'LFPclass') && strcmp(Expt.Header.LFPclass,'int')
        return;
    end    
    
    if ~isfield(Expt.Header,'lfpscale')
        for j = 1:length(Expt.Trials)
            lfprange(j,:) = minmax(Expt.Trials(j).LFP(~isnan(Expt.Trials(j).LFP)));
        end
        lfprange = minmax(lfprange(:));
        lfpscale =  max(abs(lfprange));
        Expt.Header.lfpscale = lfpscale;
        maxint = 32000;
        Expt.Header.lfpscale(2) = maxint;
    end
    maxint = Expt.Header.lfpscale(2);

    for j = 1:length(Expt.Trials)
        x = Expt.Trials(j).LFP;
        x = round(x .*maxint./lfpscale);
        x(isnan(x)) = maxint+2;
        iLFP{j} = int16(x);
    end
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).LFP = iLFP{j};
    end
    Expt.Header.LFPclass = 'int';
    return;
elseif fixdouble && isfield(Expt.Header,'lfpscale')
    maxint = Expt.Header.lfpscale(2);
    for j = 1:length(Expt.Trials)
        x = Expt.Trials(j).LFP;
        dLFP{j} = (double(x).*Expt.Header.lfpscale(1)./maxint);
        dLFP{j}(x == maxint+2) = pad;
    end
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).LFP = dLFP{j};
    end
    Expt.Header.LFPclass = 'double';
    return;
end

for j = 1:length(Expt.Trials)
    sizes(j,:) = size(Expt.Trials(j).LFP);
    if isempty(Expt.Trials(j).lfptime)
        starts(j) = 0;
    else
    starts(j) = Expt.Trials(j).Start(1) - Expt.Trials(j).lfptime;
    end
end

if max(sizes(:,2)) > MAXPROBES
    flip = 1;
    sizes = fliplr(sizes);
else
    flip = 0;
end
maxl = floor(prctile(sizes(:,1),90));
ts = median(starts);

if isfield(Expt.Header,'CRsamplerate')
    samplerate = Expt.Header.CRsamplerate .* 1000;
    if samplerate < 1 %must have been in 1/10 of ms
        samplerate = samplerate * 10;
    end
elseif isfield(Expt.Header,'LFPsamplerate')
    samplerate = Expt.Header.LFPsamplerate .* 10000;
else
    samplerate = 10;
end

%samplerate is samplerate in ms.

if isfield(Expt.Header,'preperiod')
    preperiod = Expt.Header.preperiod;
else
    preperiod = ts;
end
    

for j = 1:length(Expt.Trials)
    if flip
        Expt.Trials(j).LFP = Expt.Trials(j).LFP';
    end
    np(j) = floor((ts - starts(j)) ./(samplerate));
    lfprange(j,:) = minmax(Expt.Trials(j).LFP(:));
    if np(j) > 0
        pre = ones(np(j),sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,pre, Expt.Trials(j).LFP);
        Expt.Trials(j).lfpvalid(1) = np(j)+1;
    elseif np(j) < 0
        t = -np(j);
        Expt.Trials(j).LFP =  Expt.Trials(j).LFP(1+t:end,:);
        Expt.Trials(j).lfpvalid(1) = 1;
    end
    npost = maxl-size(Expt.Trials(j).LFP,1);
    if npost > 0
        post = ones(npost, sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,Expt.Trials(j).LFP,post);
        Expt.Trials(j).lfpvalid(2) = maxl-npost;
    elseif npost <= 0
        Expt.Trials(j).LFP = Expt.Trials(j).LFP(1:maxl,:);
        Expt.Trials(j).lfpvalid(2) = maxl;
    end
    nposts(j) = npost;
    newsizes(j,:) = size(Expt.Trials(j).LFP);
    if needft
       Expt.Trials(j).FTlfp = fft(Expt.Trials(j).LFP);
    elseif isfield(Expt.Trials,'FTlfp')
        if size(Expt.Trials(j).FTlfp,1) < maxl
            Expt.Trials(j).FTlfp(end+1:maxl,:) = NaN;
        end
    end
end
    Expt.Header.lfplen = min(newsizes(:,1));

Expt.Header.LFPtimes = (([1:Expt.Header.lfplen] .* samplerate) - (preperiod)) * 10;
if needft
    Expt.Header.LFPfreq = [1:Expt.Header.lfplen] .* 1./(Expt.Header.lfplen .* Expt.Header.LFPsamplerate);
end
if fixdouble ~= 1
    Expt = FixLFPTrials(Expt,varargin{:},'int');
else
    Expt.Header.LFPclass = 'double';
end
max(np);