function Expt = FixLFPTrials(Expt, varargin)
%Expt = FixLFPTrials(Expt)
%Makes sure that All Trials.LFP fields have the same length, padding
%front and back with NaNs;

MAXPROBES = 96;
pad = NaN;
needft = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',6)
        pad = 0;
    elseif strncmpi(varargin{j},'ft',2)
        needft = 1;
    elseif strncmpi(varargin{j},'fixspike',6)
        Expt = FixLFPSpike(Expt,varargin{:});
        return;
    end
    j = j+1;
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
elseif isfield(Expt.Header,'LFPsamplerate')
    samplerate = Expt.Header.LFPsamplerate .* 10000;
else
    samplerate = 10;
end

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
    end
end
Expt.Header.lfplen = min(newsizes(:,1));
Expt.Header.LFPtimes = ([1:Expt.Header.lfplen] .* samplerate) - (preperiod);
if needft
    Expt.Header.LFPfreq = [1:Expt.Header.lfplen] .* 1./(Expt.Header.lfplen .* Expt.Header.LFPsamplerate);
end
max(np);