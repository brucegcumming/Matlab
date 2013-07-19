function LFP = ScaleLFP(LFP, scale)


for j = 1:length(LFP.Trials)
   LFP.Trials(j).LFP(:,1:16) = LFP.Trials(j).LFP(:,1:16) .* scale;
end
    