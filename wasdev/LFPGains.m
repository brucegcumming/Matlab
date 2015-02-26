function gains = LFPGains(LFP, varargin)

for j = 1:length(LFP.Trials)
    nch = size(LFP.Trials(j).LFP,2);
pwr(j,1:nch) = std(LFP.Trials(j).LFP);
end
gains = mean(pwr);