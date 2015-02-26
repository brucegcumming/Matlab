function DATA = GetNS5Spikes(DATA, eid, probe);

if isfield(DATA.AllData.Spikes,'expt') && DATA.AllData.Spikes.expt == eid && DATA.AllData.Spikes.probe == probe
return;
end
highpass = 100;
if isempty(strfind(path,'BlackRock'))
addpath([GetFilePath('bgcmatlab') '/BlackRock']);
end
set(DATA.toplevel,'name','Loading Spikes from ns5');
drawnow;
spts = [-16:20];
gid = find(DATA.grididx.expt == eid);
if DATA.state.usensx ==2 %use NEV files with spikes in
nsname = [DATA.grididx.datdir '/' DATA.grididx.names{gid}];
fprintf('Loading %s\n',nsname);
nev = openNEV(nsname, 'read', 'compact','nomat','nosave');
id = find(nev.Data.Spikes.Electrode == probe);
Spikes.values = double(nev.Data.Spikes.Waveform(id,:));
Spikes.times = DATA.grididx.toff(gid)+double(nev.Data.Spikes.TimeStamp(id))/3.0000237;
Spikes.nevtimes = nev.Data.Spikes.TimeStamp(id);
else
nsname = strrep([DATA.grididx.datdir '/' DATA.grididx.names{gid}],'.nev','.ns5');
nsx = openNSx(nsname, 'read', ['e:' num2str(probe)]);
V = double(nsx.Data);
DATA.ArrayConfig = GetArrayConfig(nsx);
for j = 1:length(DATA.ArrayConfig.id)
DATA.probenames{j} = sprintf('%d(E%d %d,%d)',j,DATA.ArrayConfig.id(j),DATA.ArrayConfig.X(j),DATA.ArrayConfig.Y(j));
end
clear nsx;
if highpass
smv = smooth(V,highpass);
V = V-smv;
end
sgn = diff(sign(diff(V)));
id = find(sgn > 0);
nspk = length(V) .* 50./30000;
prc = 100 * nspk./length(id);
th = prctile(V(id+1),prc);
tid = find(V(id+1) < th);
id = id(tid)+1;
id = id(id < length(V)-max(spts));
allid = repmat(id,length(spts),1) + repmat(spts',1,length(id));
Spikes.values = reshape(V(allid),[size(allid)])';
Spikes.times = id'/3 + DATA.grididx.toff(gid);
end
scale = prctile(abs(Spikes.values(:)),99.99);
Spikes.values = Spikes.values .*5./scale;
Spikes.spkscale = scale/5;
Spikes.codes = zeros(length(id),1);
Spikes.dVdt = diff(Spikes.values,1,2);
Spikes.probe = probe;
Spikes.expt = eid;
DATA.AllData.Spikes = Spikes;
DATA = CalcClusterVars(DATA,  1:length(DATA.AllData.Spikes.codes),'force');
DATA.spklist = [];
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).cx = DATA.Spikes.cx;
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).cy = DATA.Spikes.cy;
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).times = DATA.AllData.Spikes.times;
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).codes = DATA.AllData.Spikes.codes;
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).suffix = DATA.currentexpt(1);
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).dips = 0;
DATA.AllClusters{DATA.currentexpt(1)}(DATA.probe).dropi = 0;
if probe > size(DATA.cluster,2)
DATA.cluster{1,probe}.touched = 0;
end
cmb.NotBusy(DATA);

