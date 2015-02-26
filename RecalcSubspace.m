function [sdf, details] = RecalcSubspace(E, tx, xval,trange, varargin)
%[sdf, count] = RecalcSubspace(E, tx, xvals,trange) %recalculate an sdf.
gap = 0;
if nargin < 4
    trange = [];
end
smoothw = 200;
ctx = tx; %second copy of tx for condtioning gap

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'conditional',5)
        j = j+1;
        ctx = varargin{j};
    elseif strncmpi(varargin{j},'gap',3)
        j = j+1;
        gap = varargin{j};
    end
    j = j+1;
end

details.times = trange(1)-smoothw:10:trange(2)+smoothw;
[a,b] = find(tx == xval);
Triggers = {[]};
Triggers{length(E.Trials)} = [];
tid = 1:length(E.Trials); %use all trials for now

useflash = zeros(1,length(a));
for j = 1:length(a)
    trial = b(j);
    if(a(j) > gap) && sum(ctx(a(j)-gap:a(j)-1,b(j))) ==0
        useflash(j) = 1;
        Triggers{trial}(end+1) = E.Trials(tid(b(j))).Start(a(j))-E.Trials(tid(b(j))).TrialStart;
    end
end
for j = 1:length(tid)
    E.Trials(tid(j)).Trigger = Triggers{j};
end
[sdf, n, nspks, spikes, t, details] = trigsdfa(E.Trials,smoothw, details.times, 'box');
spkt = find(details.spikebins >= trange(1) & details.spikebins < trange(2));
details.count = sum(spikes(spkt))./n;
details.rate = details.count * 10000./diff(trange);
details.usetrials = a(find(useflash==1));
details.useframes = b(find(useflash==1));