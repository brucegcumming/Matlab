function [lfp,n] = CalcLFPPulse(Expt, All, varargin)

lfp = []
n = 0;

plottype = 0;
trange(1) = Expt.Trials(uset(1)).Start(1);
trange(2) = Expt.Trials(uset(end)).End(end);

cid = find(All.Events.codes(:,1) == 'c' & ...
    All.Events.times > trange(1) & ...
    All.Events.times < trange(2));

if isempty(cid)
    return;
end

Expt = LoadSpike2LFP(Expt);

spks = [];
for j = 1:length(cid)
    spkid = find(All.Spikes.times > All.Events.times(cid(j)) & ...
        All.Spikes.times < All.Events.times(cid(j))+ 4000);
    id = find(max(All.Spikes.values(spkid,:)') > 3.5);
    spkid = spkid(id);
    spks = [spks spkid'];
    if plottype == 1
        plot(All.Spikes.values(spkid,:)');
        hold on;
    end
end

if ~isempty(spks)
    [lfp, n] = PulseTrigLFP(All.Spikes.times(spks),Expt.Pulses,20,Expt.Header);
    if plottype = 1
    plot(lfp.*10,'linewidth',2);
    title(sprintf('%d pulses',n));
    end
end



function [lfpavg,n] = PulseTrigLFP(times, pulses, w, Header)
% lfp = PulseTrigLFP(times, pulses, w, Header)

samplerate = 1./(Header.LFPsamplerate * 10000);
lfpavg = zeros(w+1,1);
for j = 1:length(times)
    [d, ii] = min(abs(times(j)-[pulses.ftime]));
    lfptimes = 1+round((times(j)-pulses(ii).ftime)* samplerate);
    lfp(:,j) = pulses(ii).LFP(lfptimes:lfptimes+w);
end
lfpavg = mean(lfp,2);
n = size(lfp,2);