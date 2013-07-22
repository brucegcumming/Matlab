function LFP = TrackLFPMains(LFP, tics)
%TrackLFPMains  adjusts an  LFP file for flucutating mains noise.
%requires that the file already be fixed. Then takes inner product of the
%average signal with  the record on each trial, and subtracts out the
%averagle proportion to that.
rate = LFP.Header.CRsamplerate .* 10000; %tics/sample
len = ceil(333.3./rate);

isi = diff(tics);
period = median(isi(isi < 200));
avg = LFP.Header.MainsNoise;
ugain = mean(avg.*avg);

for j = 1:length(LFP.Trials)
    kid = 1;
    id = [];

    start = LFP.Trials(j).ftime;
    last = start+ length(LFP.Trials(j).LFP) .* rate;
    t = start+ [1:length(LFP.Trials(j).LFP)] .* rate;
    
    mid = find(tics > LFP.Trials(j).ftime & tics < last+170);
        
    for k = 1:size(LFP.Trials(j).LFP,1)
        while t(k) > tics(mid(kid)) & kid < length(mid)
            kid = kid+1;
        end
        if kid > length(mid)
            fprintf('Trial longer than Mains pulses at %.3f\n', tics(mid(end))./10000);
            id(k) =  NaN;
        else
            p = len .* (t(k)-tics(mid(kid)))./period;
            id(k) = len + ceil(p);
            rem(k) = mod(p,1);
        end
    end
nid = find(id <=0);
id(nid) = id(nid)+len;
nid = find(id <=0);
if ~isempty(nid)
    fprintf('Trial %.3f precedes Mains pulses at %.3f\n', t(k)./10000,tics(1)./10000);
    id = 1+mod(id-1,len);
end
nid = find(id > len);
if ~isempty(nid)
    fprintf('Trial %d (%.3f - %.3f) %.3f long gap before %.3f\n', j,start./10000,last./10000,t(k)./10000,tics(mid(end))./10000);
    id = 1+mod(id-1,len);
end
ids{j} = id(find(~isnan(id)));
if isempty(ids{j})
    fprintf('Trial %d No ids\n', j);
end

pwr = LFP.Trials(j).LFP .* avg(ids{j},:);
pwrs(j,:) = mean(pwr)./ugain;
end
spwr = smooth(pwrs',10,'gauss')';
imagesc(cov(pwrs));
for j = 1:length(LFP.Trials)
for k = 1:size(LFP.Trials(j).LFP,2)
LFP.Trials(j).LFP(:,k) = LFP.Trials(j).LFP(:,k) - (avg(ids{j},k) .* spwr(j,k));
end
end


LFP.MainsPwrs = spwr;