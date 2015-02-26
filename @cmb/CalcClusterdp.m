function dp = CalcClusterdp(DATA, cl)  %old method using ROC of radial distances.
expspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
if strncmp(DATA.filetype,'Grid', 4)
ctype = 1;
elseif isfield(DATA,'AllClusters')
ctype = 1;
else
ctype = 2;
end
if length(expspks) < 10
dp = 0;
return;
end
if isfield(DATA,'AllSpikes')
id = find(DATA.AllSpikes{DATA.probe}.codes(expspks,ctype) == cl-1);
nid = find(DATA.AllSpikes{DATA.probe}.codes(expspks,ctype) ~= cl-1);
elseif isfield(DATA,'AllClusters')
probe = GetProbe(DATA,DATA.currentexpt(1), DATA.probe);
id = find(DATA.AllClusters{DATA.currentexpt(1)}(probe).codes(expspks,ctype) == cl-1);
nid = find(DATA.AllClusters{DATA.currentexpt(1)}(probe).codes(expspks,ctype) ~= cl-1);
else
id = find(DATA.AllData.Spikes.codes(expspks,ctype) == cl-1);
nid = find(DATA.AllData.Spikes.codes(expspks,ctype) ~= cl-1);
end
rx = std(DATA.Spikes.cx(expspks(id)));
cx = mean(DATA.Spikes.cx(expspks(id)));
ry = std(DATA.Spikes.cy(expspks(id)));
cy = mean(DATA.Spikes.cy(expspks(id)));
dx = DATA.Spikes.cx(expspks(id))- mean(DATA.Spikes.cx(expspks(id)));
dy = DATA.Spikes.cy(expspks(id))- mean(DATA.Spikes.cy(expspks(id)));
dc = abs(dx+i*dy);
dx = DATA.Spikes.cx(expspks)- mean(DATA.Spikes.cx(expspks(id)));
dy = DATA.Spikes.cy(expspks)- mean(DATA.Spikes.cy(expspks(id)));
d = abs(dx+i*dy);
[s, did] = sort(d);
dc = d;
dc(id) = 0;
d(nid) = 0;
% need to order these somehow first before doing fpos/drate
% ? rank order d before setting id, nid to zero.

%detection rate
drate = cumsum(dc(did)) ./ sum(dc);
%false positive rate
fpos = cumsum(d(did)) ./ sum(d);
%area gives ROC
dp = trapz(fpos,drate);


