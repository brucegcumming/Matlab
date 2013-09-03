function Expt = FixLFPTrials(Expt, varargin)
%Expt = FixLFPTrials(Expt)
%Makes sure that All Trials.LFP fields have the same length, padding
%front and back with NaNs;

MAXPROBES = 96;
pad = NaN;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',6)
        pad = 0;
    elseif strncmpi(varargin{j},'fixspike',6)
        Expt = FixLFPSpike(Expt,varargin{:});
        return;
    end
    j = j+1;
end

for j = 1:length(Expt.Trials)
    sizes(j,:) = size(Expt.Trials(j).LFP);
    starts(j) = Expt.Trials(j).Start(1) - Expt.Trials(j).lfptime;
end

if max(sizes(:,2)) > MAXPROBES
    flip = 1;
    sizes = fliplr(sizes);
else
    flip = 0;
end
maxl = max(sizes(:,1));
ts = min(starts);
for j = 1:length(Expt.Trials)
    if flip
        Expt.Trials(j).LFP = Expt.Trials(j).LFP';
    end
    np(j) = floor((starts(j)-ts) ./(Expt.Header.CRsamplerate * 10000));
    if np(j) > 0
        pre = ones(np(j),sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,pre, Expt.Trials(j).LFP);
    end
    npost = maxl-size(Expt.Trials(j).LFP,1);
    if npost > 0
        post = ones(npost, sizes(j,2)).* pad;
        Expt.Trials(j).LFP = cat(1,Expt.Trials(j).LFP,post);
    elseif npost < 0
        Expt.Trials(j).LFP = Expt.Trials(j).LFP(1:maxl,:);
    end
    newsizes(j,:) = size(Expt.Trials(j).LFP);
end
Expt.Header.lfplen = min(newsizes(:,1));
Expt.Header.lfptimes = (10000 .* [1:Expt.Header.lfplen] .* Expt.Header.LFPsamplerate) - (Expt.Header.preperiod);
max(np);