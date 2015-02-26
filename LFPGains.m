function gains = LFPGains(LFP, varargin)

if ~isfield(LFP.Trials,'LFP') %load failed
    gains = 0;
    return;
end

for j = 1:length(LFP.Trials)
    nch = size(LFP.Trials(j).LFP,2);
pwr(j,1:nch) = std(LFP.Trials(j).LFP);
end
gains = mean(pwr);