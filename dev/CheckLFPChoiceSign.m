function CheckLFPChoiceSign(LFP, probe)

fid =1:250;
uid = find([LFP.Trials.Dc] == 0 & [LFP.Trials.RespDir] > 0);
did = find([LFP.Trials.Dc] == 0 & [LFP.Trials.RespDir] < 0);
a = abs(cat(3,LFP.Trials(uid(1:end)).FTlfp));
upchoicepwr = mean(a(fid,7,:),3);
a = abs(cat(3,LFP.Trials(did(1:end)).FTlfp));
dnchoicepwr = mean(a(fid,7,:),3);

mo = mean([LFP.Trials.ori]);
did = find([LFP.Trials.Dc] > 0.1 & [LFP.Trials.ori] > mo);
uid = find([LFP.Trials.Dc] > 0.1 & [LFP.Trials.ori] < mo);
a = abs(cat(3,LFP.Trials(uid(1:end)).FTlfp));
upsigpwr = mean(a(fid,7,:),3);
uprespdir = mean([LFP.Trials(uid).RespDir]);
dnrespdir = mean([LFP.Trials(did).RespDir]);
upstim = mean([LFP.Trials(uid).ori]);
dnstim = mean([LFP.Trials(did).ori]);
a = abs(cat(3,LFP.Trials(did(1:end)).FTlfp));
dnsigpwr = mean(a(fid,7,:),3);
len = length(LFP.Header.LFPtimes);
freqs = (0:len-1)/(len * LFP.Header.LFPsamplerate);

hold off;
plot(freqs(fid),upsigpwr-dnsigpwr);
hold on;
plot(freqs(fid),upchoicepwr-dnchoicepwr,'r');
legend('signal','choice');
title(sprintf('%.0f-%.0f',upstim,dnstim));
