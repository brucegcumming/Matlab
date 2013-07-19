function result = PlotSpikeFile(spkfile, varargin)
% PlotSpikeFile(file, varargin)
% Plots the spike waveform data from a single file on disk
% useful for looking at spikes made by autocutting in AllVPcs
Spks = [];
AllSpikes = {};
preperiod = 500;
postperiod = 500;
result = [];
Expt = [];
C = [];
excludetrials = [];
probe =1;
files{1} = spkfile;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        if isfield(varargin{j}, 'values')
            Spks = varargin{j};
        elseif isfield(varargin{j},'Trials')
            Expt = varargin{j};
        elseif isfield(varargin{j},'probe')
            C = varargin{j};
            probe = C.probe;
        end
    elseif strncmp(varargin{j},'exclude',6)
        j = j+1;
        excludetrials = varargin{j};
    elseif strncmp(varargin{j},'spkfile',6)
        j = j+1;
        files{length(files)+1} = varargin{j};
    end
    j = j+1;
end


    
[probes(1), ex, xs ] = Name2Probe(spkfile);
if length(files) > 1
    for j = 1:length(files)
        probes(j) = Name2Probe(files{j});
        AllSpikes{probes(j)} = ReadSpikeFile(files{j});
    end
else
    Spikes = ReadSpikeFile(spkfile);
end
probe = probes(1);

cfile = strrep(spkfile,'/Spikes','');
cfile = sprintf('%s/Expt%d%sClusterTimesDetails.mat',fileparts(cfile),ex,xs);
afile = sprintf('%s/Expt%d%sAutoClusterTimesDetails.mat',fileparts(cfile),ex,xs);
if exist(cfile)
    load(cfile);
end
if exist(afile)
    a = load(afile);
    if exist('ClusterDetails')  %have both
        for j = 1:length(a.ClusterDetails)
            if j > length(ClusterDetails) || isempty(ClusterDetails{j})
                ClusterDetails{j} = a.ClusterDetails{j};
            end
        end
    else
        ClusterDetails = a.ClusterDetails;
    end
end
if exist('ClusterDetails') && isfield(ClusterDetails{probe},'xy')
    Spikes.xy = ClusterDetails{probe}.xy;
    DATA.Spikes.cx = Spikes.xy(:,1);
    DATA.Spikes.cy = Spikes.xy(:,2);
end

if length(AllSpikes)
    AllSpikes{probes(1)}.xy = Spikes.xy;
    Spks = AllSpikes{probes(1)};
else
    Spks = Spikes;
end
if ~isfield(Spks,'values') %failed to load - 
    return;
end
if size(Spks.xy,1) < size(Spks.values,1)
    fprintf('Size mismatch %d XY vals, %d Spike Vals\n',size(Spks.xy,1) ,size(Spks.values,1));
end
if isfield(C,'clst') && size(C.clst,1) == size(Spks.codes,1)
    Spks.codes(:,1) = C.clst -1;
    Spks.codes(:,2) = C.clst -1;
    if length(AllSpikes)
     AllSpikes{probes(1)}.codes(:,1) = C.clst-1;        
     AllSpikes{probes(1)}.codes(:,2) = C.clst-1;        
    end
end
    
t = Spks.times(1);
nt = 1;
tw = 10000;
if isempty(Expt)
while t < max(Spks.times)
    Expt.Trials(nt).Start = t;
    Expt.Trials(nt).End = t + tw;
    Expt.Trials(nt).Trial = nt;
    Expt.Trials(nt).id = nt;
    t = t+tw+preperiod+postperiod;
    id = find(Spks.times > t);
    if length(id)
    t = Spks.times(id(1))-preperiod;
    end
    nt = nt+1;
end
end
Expt.Header.trange = minmax(Spks.times);
Expt.Header.Name = spkfile;
Expt.Header.cellnumber = 0;
DATA.plot.syncoverlay = 0;
if length(AllSpikes)
    DATA.AllSpikes = AllSpikes;
    DATA.plot.syncoverlay = 1;
else
DATA.AllData.Spikes = Spks;
end
DATA.Expts{1} = Expt;
DATA.spklist = 1:length(Spks.times);
DATA.AllData.pcs = [];
DATA.tag.clusterxy = 'Spike XY Plot';
DATA.tag.spikev = 'Spike Waveforms';
DATA.tag.top = DATA.tag.clusterxy;
DATA.state.uselfp = 0;
DATA.plot.clusterX = 1;
DATA.plot.clusterY = 2;
DATA.plot.clusterXrange = [0 40];
DATA.plot.clusterYrange = [0 5];
DATA.plot.clusterZrange = [0 5];
DATA.plot.clusterZ = 3;
DATA.plot.SpikeMaxV = max(abs(Spks.values(:)));
DATA.plot.SpikeVsep = 2;
DATA.plot.dvdt = 0;
DATA.plot.showwave = 0;
DATA.plot.prettyfigs = 1;
DATA.plot.synccluster = 0;
DATA.probe = probe;
DATA.currenttrial = 1;
DATA.plot.showsync = 0;
DATA.state.fixrange = 0;
DATA.state.preperiod = preperiod;
DATA.state.postperiod = postperiod;
DATA.syncsign = 6;
DATA.state.recut = 1;
DATA.cluster = [];
DATA.currentcluster = 1;
DATA.firsttrial = 1;
DATA.clusterArange = [7:11];
DATA.clusterBrange = [12:18];
DATA.clusterErange = [1:40];
DATA.plot.nodc = 1;
DATA.plot.voffsets(1) = 0; 
if length(probes) > 1
DATA.xprobes = probes(2:end);
np = (length(files)-1);
for j = 1:length(probes)
    range(j,1) = prctile(DATA.AllSpikes{probes(j)}.values(:),1);
    range(j,2) = prctile(DATA.AllSpikes{probes(j)}.values(:),99);
end
id = find(probes > probes(1));
id = sort(probes(id));
for j = 1:length(id);
    p = find(probes == id(j));
    q = find(probes == id(j)-1);
    DATA.plot.voffsets(p) = DATA.plot.voffsets(q)+range(p-1,2)-range(p,1);
end
id = sort(find(probes < probes(1)));
id = sort(probes(id));
for j = 1:length(id);
    p = find(probes == id(j)); %this probe
    q = find(probes == id(j)+1); %the probe above
    DATA.plot.voffsets(p) = DATA.plot.voffsets(q)+range(q,1)-range(q,2);
end
DATA.plot.voffsets = DATA.plot.voffsets .* 1.2;
end
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.ptsize = 1;
    DATA.plot.xcorr = 0;
DATA.probelist = probes;
DATA.plot.autoscale = 1;
%DATA.AllData.Trialids = [Expt.Trials.Trial];
%DATA.Expts{1} = Expt;
DATA.explabels{1} = spkfile;
DATA.currentexpt = 1;

DATA.plot.showwave = 0;
DATA.vstep = 1;
DATA.densityplot = 0;
DATA.plot.setptsize = 0;
DATA.spikelist = [0:4];
DATA.sids{1} = 0;
DATA.sids{2} = 0;
for j = 1:length(excludetrials)
    DATA.Expts{1}.Trials(excludetrials(j)).Trial  =  abs(DATA.Expts{1}.Trials(excludetrials(j)).Trial)  .* -1; 
end
DATA = SpoolSpikes(DATA, varargin{:});
if isfield(DATA,'dprime')
    result.dprime = DATA.dprime;
end
if isfield(DATA,'dprimes')
    result.dprimes = DATA.dprimes;
end
if isfield(DATA,'dprimet')
    result.dprimet = DATA.dprimet;
end
result.toplevel = DATA.toplevel;
result.svfig = DATA.svfig;


function [probe, ex, xs] = Name2Probe(name)

cfile = strrep(name,'/Spikes','');
id = regexp(cfile,'p[0-9]*t[0-9]*');
if length(id)
probe = sscanf(cfile(id+1:end),'%d');
j = id+2;
while cfile(j) ~= 't'
    j = j+1;
end
ex = sscanf(cfile(j+1:end),'%d');
id = regexp(cfile,'p[0-9]*t[0-9]*a');

if length(id)
    xs = 'a';
else
    xs = [];
end
end


function Spikes = ReadSpikeFile(spkfile)
if ~exist(spkfile)
    fprintf('Cannot read %s\n',spkfile);
    Spikes = [];
    return;
end
load(spkfile);
if size(Spikes.values,2) > 100
    Spikes.values = Spikes.values';
end
    Spikes.times = Spikes.times .* 10000;
    Spikes.times = reshape(Spikes.times,length(Spikes.times),1);
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end

