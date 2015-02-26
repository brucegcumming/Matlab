function Expt = FilterLFPTrials(Expt)

lfplen = Expt.Header.lfplen;
for j = 1:length(Expt.Trials)
    [Expt.Trials(j).LFP Expt.Trials(j).FTlfp removed(:,j)] = FilterLFPa(Expt.Trials(j).LFP(1:lfplen),Expt.Header.CRsamplerate,Expt.Trials(j).Start(1));
end
return;



function [lfp, ft, removed] = FilterLFPa(lfp, samplerate, start)

mainsperiod = 166.6155;
%mainsperiod = 166;
%mainsperiod = 332;
mainsperiod = 166.71;
id = 1;

%
% first caculate the mean periodic variation with period mainsperiod. 
% To avoid aliasing, build two averages with value added to bins either side
% of correct one. 
%
lfppad = [lfp; zeros(length(lfp) * 4,1)];
ftfrq = (0:length(lfppad)-1)/(length(lfppad)/(10000 * samplerate));
pwr = abs(fft(lfppad))';
[a,b] = min(abs(ftfrq - 60));
ids = b-3:b+3;
mf = sum(pwr(ids).* ftfrq(ids))./sum(pwr(ids));

t = start+((0:length(lfp)-1))./samplerate;
amean = zeros(1,1+floor(mainsperiod * samplerate));
abins = 1+round(mod(t,mainsperiod).*samplerate);
avgl = max(abins);
for bin = 1:max(abins);
        amean(bin) = mean(lfp(find(abins == bin)));
end
nm = 1;
mps = [166.0:0.01:168.0];

for mp = mps
    amean = zeros(1,1+floor(mp * samplerate));
    abins = 1+round(mod(t,mp).*samplerate);
    for bin = 1:max(abins);
        amean(bin) = mean(lfp(find(abins == bin)));
    end
    ameans(:,nm) = amean;
    pwrs(nm) = var(lfp - amean(abins)');
    nm = nm+1;
end
    
[mins(1) mint(1)] = min(lfp(1:avgl));
tm = mint(1);
cbins(tm:tm+avgl-1) = 1:avgl;
if tm > 1
cbins(1:tm) = avgl-tm:tm;
end
j = 2;
while tm < (length(lfp)-avgl)
    [mins(j) step] = min(lfp((tm+2):tm+avgl));
    tm= tm+step+1;
    mint(j) = tm;
    cbins(tm:tm+avgl-1) = 1:avgl;
    j = j+1;
end
cbins = cbins(1:length(abins));
dbins = 1+ round(mod(t, ceil(mainsperiod * samplerate)));
for bin = 1:max(abins);
    amean(bin) = mean(lfp(find(abins == bin)));
    cmean(bin) = mean(lfp(find(cbins == bin)));
end
for bin = 1:max(dbins);
    dmean(bin) = mean(lfp(find(dbins == bin)));
end
lfp = lfp - amean(abins)';
removed = amean;
%lfp = lfp - amean(abins)'/2;
%lfp = lfp - dmean(dbins)'/2;

ft = fft(lfp);
