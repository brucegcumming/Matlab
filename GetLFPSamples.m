function LFP = GetLFPSamples(Expt, times)
%LFP = GetLFPSamples(Expt, times) Extract a matrix of LFP values from
%Trials
LFP = [];

if isfield(Expt.Header,'LFPsamplerate')
    samplerate = Expt.Header.LFPsamplerate .* 10000;
elseif isfield(Expt.Header,'CRsamplerate')
    samplerate = Expt.Header.CRsamplerate .* 10000;
else
    samplerate = 10;
end
    
for j = 1:length(Expt.Trials)
    if isfield(Expt.Trials,'lfptime')
        start = Expt.Trials(j).lfptime;
    else
        start = Expt.Trials(j).ftime;
    end
    ts = [1:size(Expt.Trials(j).LFP,1)] .* samplerate + start - Expt.Trials(j).Start(1);
    id = find(ts > times(1) & ts < times(end));
    toff = (ts(1) - times(1)) ./ (samplerate .* 10000);
    if toff < 1
        toff = 1;
    end
    LFP(toff:toff+length(id)-1,:,j) = Expt.Trials(j).LFP(id,:);
end