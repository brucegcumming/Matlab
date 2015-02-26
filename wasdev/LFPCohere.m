function [cx, details] = LFPCohere(LFP)

nt = 0;
LFP = CheckLFP(LFP,'fix'); %make sure sizes are constant
for t =length(LFP.Trials):-1:1
    if sum(var(LFP.Trials(t).LFP)) > 0
    for j = size(LFP.Trials(t).LFP,2):-1:1
        for k = size(LFP.Trials(t).LFP,2):-1:1
        cx(:,j,k) = mscohere(LFP.Trials(t).LFP(:,j),LFP.Trials(t).LFP(:,k));
        end
    end
    if nt == 0
        csum = cx;
    else
        csum = csum+cx;
    end
    nt = nt+1;
    lens(nt) = size(LFP.Trials(t).LFP,1);
    end
end

len = min(lens);
cx = csum./nt;
for j = 1:size(cx,1)
    details.absdiff(j,:) = mean(squeeze(abs(diff(cx(j,:,:)))),2);
    details.freqs = (0:len-1)/(len * LFP.Header.LFPsamplerate);
end