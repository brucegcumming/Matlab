function [DATA, counts] = CountSpikesB(DATA, expid, pid, spikelist, varargin)
%use histc to find spkies/trials alignment.
%
%expid is hte indext to DATA.expts. May not match Expt Suffix number
Spks = DATA.AllClusters{expid}(pid(1));
codes = zeros(size(Spks.codes(:,1)));
if size(Spks.codes,2) == 1
ctype = 1;
else
ctype = 2;
end
for j = 1:length(spikelist)
codes(Spks.codes(:,ctype)==spikelist(j)) = 1;
end
nt = length(DATA.Expts{expid}.Trials);
for t = 1:nt
x(t*2-1)  = DATA.Expts{expid}.Trials(t).Start(1) - DATA.state.preperiod;
x(t*2)  = DATA.Expts{expid}.Trials(t).End(end) + DATA.state.postperiod;
end
%since histc needs monotonicaly increasing x, if trial(n)-preperiod < trial(n-1)+postperiod
%this fails.  So do even and odd trials separately.
aid = 1:2:nt;
bid(1:2:length(aid)*2) = aid*2-1;
bid(2:2:length(aid)*2) = 2*aid;
aid = 2:2:nt;
cid(1:2:length(aid)*2) = aid*2-1;
cid(2:2:length(aid)*2) = 2*aid;
[a,b] = histc(Spks.times,x(bid));
[a,c] = histc(Spks.times,x(cid));
for t = 1:length(DATA.Expts{expid}.Trials)
trial = DATA.Expts{expid}.Trials(t);
if bitand(t,1)
mid = find(b == t);
else
mid = find(c == t-1);
end
id = find(codes(mid,1) ==1);
oid = find(codes(mid,1) ==0);
DATA.Expts{expid}.Trials(t).Spikes = double(round(Spks.times(mid(id)) - trial.Start(1)));
spks = DATA.Expts{expid}.Trials(t).Spikes;
counts(t) = sum(spks > 500 &  spks < trial.dur + 500);
DATA.Expts{expid}.Trials(t).count = counts(t);
DATA.Expts{expid}.Trials(t).OSpikes = round(Spks.times(mid(oid)) - trial.Start(1));
DATA.Expts{expid}.Trials(t).Ocodes = codes(mid(oid));
end

