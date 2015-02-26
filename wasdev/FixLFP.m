function [Expt, spksub] = FixLFP(Expt, AllData)
% FixLPF removes the effect os a spike on the LFP
%
submode = 3;
calcfakelfp = 0;
gain = 1;
if submode ~= 3
spkid = find(Spikes.codes(:,2) ==1 & Spikes.times * 10000> Expt.Trials(1).Start(1)...
    & Spikes.times * 10000 < Expt.Trials(end).End(end));
allspkid = find(Spikes.times * 10000> Expt.Trials(1).Start(1)...
    & Spikes.times * 10000 < Expt.Trials(end).End(end));


spkwv = mean(Spikes.values(spkid,:)) .* gain;
end
if submode == 1
for j = 1:length(Expt.Trials)
    starts = round(Expt.Trials(j).lfpo+(Expt.Trials(j).Spikes./(Expt.Header.LFPsamplerate * 10000)));
    for k = starts'
        Expt.Trials(j).LFP(k:k+length(spksub)-1) = Expt.Trials(j).LFP(k:k+length(spksub)-1) - spksub'; 
        lfpavg(k:k+length(spksub)-1) = lfpavg(k:k+length(spksub)-1)+spksub'; 
        lfpn(k:k+length(spksub)-1) = lfpn(k:k+length(spksub)-1)+1; 
    end
end
elseif submode == 2
    spksub = mean(spkwv) .* [0.2 0.6 0.2];
for j = 1:length(Expt.Trials)
    starts = toff+round(Expt.Trials(j).lfpo+(Expt.Trials(j).Spikes./(Expt.Header.LFPsamplerate * 10000)));
    for k = starts'
        Expt.Trials(j).LFP(k:k+length(spksub)-1) = Expt.Trials(j).LFP(k:k+length(spksub)-1) - spksub'; 
    end
    if calcfakelfp
    id = find(spktimes>Expt.Trials(j).Start(1) & spktimes < Expt.Trials(j).End(end));
    starts = toff+Expt.Trials(j).lfpo+round((spktimes(id)-Expt.Trials(j).Start(1))./...
        (Expt.Header.LFPsamplerate * 10000));
    lfpavg(starts) = lfpavg(starts)+spkmeans(id)';
        lfpn(starts) = lfpn(starts)+1; 
    end
end
elseif submode == 3
    [a,n, spksub] = CalcLFPPulse(Expt, AllData);
    Expt.Header.LFPspksub = spksub;
    idm = find(spksub > max(spksub)/20);
    spksub = spksub(1:idm(end));
    for j = 1:length(Expt.Trials)
        starts = round(Expt.Trials(j).lfpo+(Expt.Trials(j).Spikes./(Expt.Header.LFPsamplerate * 10000)));
        starts = starts(find(starts < length(Expt.Trials(j).LFP) - (length(spksub))))-1;
        for k = starts'
            Expt.Trials(j).LFP(k:k+length(spksub)-1) = Expt.Trials(j).LFP(k:k+length(spksub)-1) - gain * spksub';
        end
    end

end

if calcfakelfp
    lfpavg = lfpavg./lfpn;
end