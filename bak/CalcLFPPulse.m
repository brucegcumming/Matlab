function [lfp,n, adj] = CalcLFPPulse(Expt, All, varargin)
%[lfp,n] = CalcLFPPulse(Expt, All, ...)

lfp = [];
adj = [];
n = 0;
plottype = 0;


j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',3)
        plottype = 1;
    end
    j = j+1;
end

if isfield(Expt.Header,'trange');
trange  = Expt.Header.trange;
trange(2) = trange(2)+1000; %allow 100ms for the pulses after receiving endexpt
else
trange(1) = Expt.Trials(1).Start(1);
trange(2) = Expt.Trials(end).End(end);
end
cid = find((All.Events.codes(:,1) == 'c' | All.Events.codes(:,3) > 5)& ...
    All.Events.times > trange(1) & ...
    All.Events.times < trange(2));

if isempty(cid)
    return;
end

Expt = LoadSpike2LFP(Expt);
if ~isfield(Expt,'pulses')
    return;
end

spks = [];
for j = 1:length(cid)
    spkid = find(All.Spikes.times > All.Events.times(cid(j)) & ...
        All.Spikes.times < All.Events.times(cid(j))+ 4000);
    pks = max(abs(All.Spikes.values(spkid,1:20)'));
    id = find(pks > median(pks) * 0.9);
    spkid = spkid(id);
    spks = [spks spkid'];
end

sui = mean(All.Spikes.values(spks,:));
sui = median(All.Spikes.values(spks,:));
[a,b] = max(abs(sui));
if sui(b) < 0  %% Cal pulse looks to be negative
    ncal = 1;
    spkid = find(All.Spikes.values(spks,b) < sui(b)/2);
    bid = find(All.Spikes.values(spks,b) > sui(b)/2);
else
    ncal = 0;
    spkid = find(All.Spikes.values(spks,b) > sui(b)/2);
    bid = find(All.Spikes.values(spks,b) < sui(b)/2);
end
bspks = spks(bid);
spks = spks(spkid);

if plottype == 1 & ~isempty(spks)
    plot(All.Spikes.values(spks,:)','r');
    hold on;
    if ~isempty(bspks)
        plot(All.Spikes.values(bspks,:)','b');
    end
end

suv = GetMeanSpike(Expt.Trials, All);
DATA.meanspk = suv;
DATA.impulse = sui;

% subtract of mean endpts, to estimate DC. i.e. if ongoing voltage is
% correlated with p(spike), do NOT want to remove this from the average
suv = suv - mean(suv([1 end]));
if ncal
spkp = FitSpike(suv,-sui);
[err,spk] = TrySpike(spkp, suv,-sui);
else
spkp = FitSpike(suv,sui);
[err,spk] = TrySpike(spkp, suv,sui);
end


if ~isempty(spks) 
    sratio = Expt.Header.LFPsamplerate./All.Spikes.interval;
    [lfp, n] = PulseTrigLFP(All.Spikes.times(spks),Expt.Pulses,20,Expt.Header);
    G = FitGauss([1:10],lfp(1:10)','freebase');
    if ncal
        G.params(3) = -G.params(3);
    end
    sfit = FitGauss([1:1/sratio:10],G.params,'eval')-G.params(4);
    c = conv(sfit,spk);
    adj = c(10:sratio:end);
    if plottype == 1
    plot(lfp.*10,'linewidth',2);

    title(sprintf('%d pulses',n));
    end
end

function suv = GetMeanSpike(Trials, All)

spks = [];
for j = 1:length(Trials)
    spkid = find(All.Spikes.times > Trials(j).Start(1) & ...
        All.Spikes.times < Trials(j).End(end)+500);
    spks = [spks spkid'];
end
suv = mean(All.Spikes.values(spks,:));


function [lfpavg,n] = PulseTrigLFP(times, pulses, w, Header)
% lfp = PulseTrigLFP(times, pulses, w, Header)

samplerate = 1./(Header.LFPsamplerate * 10000);
lfpavg = zeros(w+1,1);
for j = 1:length(times)
    [d, ii] = min(abs(times(j)-[pulses.ftime]));
    lfptimes = 1+round((times(j)-pulses(ii).ftime)* samplerate);
    if lfptimes > 0 & lfptimes < length(pulses(ii).LFP) - w
        lfp(:,j) = pulses(ii).LFP(lfptimes:lfptimes+w);
    end
end
lfpavg = mean(lfp,2);
n = size(lfp,2);


function [err, insp, diffs] = TrySpike(x, spike, impulse, varargin)

insp = zeros(size(spike));
imp = insp;
imp(10) = 1;
mode = 2;

if mode == 1
%insp(round(x(1))) = x(2);
npts = round(x(4));
npost = round(x(5));
if npts > 1
    insp(round(x(1))+[1:npts]) = interp1([1 npts],[x(2) x(3)],[1:npts]);
else
    insp(round(x(1))+1) = x(2);
    insp(round(x(1))+2) = x(3);
end
    
insp(round(x(1))+npts-1+[1:npost]) = interp1([1 npost],[x(3) 0],[1:npost]);
elseif mode == 2 % exponential, starts at x(1), voltage x(2), decay constant x(3)
    a = round(x(1)); % start time;
    if a < 1
        a = 1;
    end
    npts = length(insp) - a;
    insp(a:end) = x(2) .* exp(-[0:npts]/x(3)) + x(2) .* exp(-[0:npts]/x(5));
end
p = conv(insp,impulse) + x(4);
p = p(10:55);
%p = p - mean(p-spike);
diffs = (p-spike);
err = sum(diffs.^2);

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',3)
        plot(p);
        hold on;
        plot(spike,'r');
        plot(insp,'g');
    end
    j = j+1;
end


function spkp = FitSpike(spike, impulse)
%fit a simple descriptive spike shape such than when convolved with the 
%impulse function it produces the observed mean spike

x(1) = 8;
x(2) = 0.1;
x(3) = 8;
x(4) = 0.1;
x(5) = 16;
options = optimset('TolFun',1e-6);
spkp = fminsearch(@TrySpike,x,options, spike, impulse);


