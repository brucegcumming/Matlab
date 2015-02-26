function [dprime, details] = CalcDprime(DATA,cluster, varargin)

trackcl = 0;
yr = DATA.Spikes.cy(DATA.spklist);
xr = DATA.Spikes.cx(DATA.spklist);
ispk = DATA.spklist;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'track',4)
trackcl = 1;
elseif strncmpi(varargin{j},'smooth',4)
j = j+1;
smw = varargin{j};
elseif strncmpi(varargin{j},'xparam',4)
j = j+1;
Spikes = DATA.AllData.Spikes;
xr = GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),varargin{j}, 1,[]);
elseif strncmpi(varargin{j},'yparam',4)
j = j+1;
Spikes = DATA.AllData.Spikes;
yr = GetSpikeVals(DATA, ispk, Spikes.values(ispk,:), Spikes.dVdt(ispk,:),varargin{j}, 1,[]);
end
j = j+1;
end
spkr = minmax(DATA.spklist);

id = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == cluster);
nid = find(DATA.AllData.Spikes.codes(DATA.spklist,2) == 0);
sy = ((yr-mean(yr(id)))./std(yr));
sx = ((xr-mean(xr(id)))./std(xr));
sd = abs(sx+i*sy);
o = [0:pi/40:pi];
for j = 1:length(o)
d = sx .* cos(o(j)) + sy.* sin(o(j));
ds{j} = d;
%
%        dprimes(j) = abs(mean(d(nid(tid)))-mean(d(id)))./sqrt(mean([var(d(nid(tid))) var(d(id))]));
dprimes(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
%        dprimes(j) = abs(mean(d(nid))-mean(d(id)))./std(d(nid));
end
[dprime, maxi] = max(dprimes);
detials.distances = ds{maxi};
details.id = id;
details.angle = o(maxi);
if trackcl
smw = 0;
d = ds{maxi};
T = DATA.Expts{1}.Trials;
for j = 1:length(T);
ispk = [];
codes = [];
for k = max([1 j-smw]):min([length(T) j+smw])
[a, ts, b] = FindSpikes(DATA,[T(k).Start(1) T(k).End(end)+500], DATA.probe, spkr);
ispk = [ispk; a];
codes = [codes; b];
end
nid = ispk(find(codes == 0));
did = ispk(find(codes == cluster));
dp(j) = abs(mean(d(nid))-mean(d(id)))./sqrt(mean([var(d(nid)) var(d(id))]));
end
details.dps = dp;
end


