function LFP = AddLFPFT(LFP, varargin)

nch = 23;
for j = 1:length(LFP.Trials)
    chs(j) = size(LFP.Trials(j).LFP,2);
end

for j = find(chs >= nch)
   LFP.Trials(j).FTlfp = fft(LFP.Trials(j).LFP);
end

