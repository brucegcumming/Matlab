function Expt = FixLFPSpike(Expt, varargin)
%Expt = FixLFPSpike(Expt)
%Removed mean spike from LFP record
%Uses MeanSpike from ClusterTimes, convolved with filter used to build LFP
%To estimate mean spike in the LFP. Then removes this mean from LFP record
%for each Spike 


MAXPROBES = 96;
pad = NaN;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'zeropad',6)
        pad = 0;
    end
    j = j+1;
end


for j = 1:length(Expt.Header.Clusters)
    C = Expt.Header.Clusters{j};
    if isfield(C,'MeanSpike')
        ms(j,:) = C.MeanSpike.ms;
        trigt(j) = C.trigt;
    end
end

ms = mean(ms,1);
subsample = Expt.Header.LFPdecimate;
G = Gauss(Expt.Header.LFPsigma,[-100:100]);
lfps = conv(ms,G);
spts = downsample([1:length(lfps)]-trigt(1),subsample);
spts = round(spts./subsample);
lfps = downsample(lfps,subsample);
p = Expt.Header.probe;
Expt.Header.LFPspike = lfps;

for j = 1:length(Expt.Trials)
    t = round((Expt.Trials(j).Spikes - Expt.Trials(j).lfpstart)/(Expt.Header.LFPsamplerate .* 10000));
    for s = 1:length(t);
        idx = t(s) + spts-1;
        idx = idx(idx > 0 & idx <= size(Expt.Trials(j).LFP,1));
        Expt.Trials(j).LFP(idx,p) = Expt.Trials(j).LFP(idx,p) - lfps(1:length(idx))'; 
    end
end
