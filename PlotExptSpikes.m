function result = PlotExptSpikes(Expt, varargin)
% PlotExptSpikes(Expt, varargin)
% Plots the spike waveforms belongint to an Expt.  Finds them
% on the disk first, includeing for cells build with dirrent 
% probes etc.
Spks = [];
AllSpikes = {};
preperiod = 500;
postperiod = 500;
result = [];
if isempty(Expt)
    result.dprime = NaN;
    return;
end
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) & isfield(varargin{j}, 'values')
        Spks = varargin{j};
    end
    j = j+1;
end
t(1) = Expt.Trials(1).Start(1);
t(2) = Expt.Trials(end).End(end);
if ~isfield(Expt.Header,'trange')
Expt.Header.trange = t;
end

spkfile = Expt.Header.Name;
if isfield(Expt,'probes')
    ps = unique(Expt.probes);
    spkfile = strrep(Expt.Header.Name,'.mat','probes.mat');
    spkfile = CheckPath(spkfile); 
    load(spkfile);
    for j = 1:length(Expt.Trials)
        a = find([probes.first] < Expt.Trials(j).Start(1)./10000 & [probes.last] > Expt.Trials(j).End(1)./10000 ...
            & [probes.probe] == Expt.probes(j));
        probeid(j) = a;
        subprobes(j) = mod(a,10);
    end

    for j = 1:length(ps)
        id = find(Expt.probes == ps(j));
        subs{j} = unique(subprobes(id));
        [a,b] = fileparts(strrep(Expt.Header.Name,'\','/'));
        clusterfile = CheckPath(sprintf('%s/%s.p%dcl.mat',a,b,ps(j)))
        clear clid;
        load(clusterfile);

        for k = 1:length(subs{j})
            pid = unique(probeid(find(Expt.probes == ps(j) & subprobes == subs{j}(k))));
            [a,b] = fileparts(strrep(Expt.Header.Name,'\','/'));
            spkfile = CheckPath(sprintf('%s/Spikes/%s',a,probes(pid).filename))
            a = load(spkfile);
            spkfile = strrep(spkfile,'\','/'); 
            
            if k == 1
            AllSpikes{ps(j)} = a.(probes(pid).var);
            AllSpikes{ps(j)}.times = AllSpikes{ps(j)}.times .*10000;
            AllSpikes{ps(j)}.codes(:,2) = clid(probes(pid).firsti:(probes(pid).firsti+probes(pid).nspk-1));
            else
                AllSpikes{ps(j)}.times = [AllSpikes{ps(j)}.times; a.(probes(pid).var).times.*10000];
                AllSpikes{ps(j)}.values = [AllSpikes{ps(j)}.values; a.(probes(pid).var).values];
                a.(probes(pid).var).codes(:,2) = clid(probes(pid).firsti:(probes(pid).firsti+probes(pid).nspk-1));
                AllSpikes{ps(j)}.codes = [AllSpikes{ps(j)}.codes; a.(probes(pid).var).codes];
            end
        end
    end
    Spks.times = [];
    Spks.values = [];
    Spks.codes = [];
    allid = [];
    for j = 1:length(Expt.probes)
        id = find(AllSpikes{Expt.probes(j)}.times > Expt.Trials(j).Start(1)-preperiod ...
            & AllSpikes{Expt.probes(j)}.times < Expt.Trials(j).End(end)+postperiod);
        allid = [allid; id];
        Spks.values = [Spks.values; AllSpikes{Expt.probes(j)}.values(id,:)];
        Spks.times = [Spks.times; AllSpikes{Expt.probes(j)}.times(id,:)];
        Spks.codes = [Spks.codes; AllSpikes{Expt.probes(j)}.codes(id,:)];
    end
    if isfield(Expt.Trials,'SpikeGain') && length(unique([Expt.Trials.SpikeGain])) > 1
        refgain = mean([Expt.Trials.SpikeGain]);
        did = find(abs(diff([Expt.Trials.SpikeGain])) > 0);
        did = [did length(Expt.Trials)];
        last = 1;
        for j = 1:length(did)
            t(1) = Expt.Trials(last).Start(1)-preperiod;
            t(2) = Expt.Trials(did(j)).End(end)+preperiod;
            id = find(Spks.times > t(1) & Spks.times < t(2));
            gain = median([Expt.Trials(last:did(j)).SpikeGain]);
            Spks.values(id,:) = Spks.values(id,:) .* refgain/gain;
            last = did(j)+1;
        end
    end
    Spks.dVdt = diff(Spks.values,1,2);
    DATA.AllData.Spikes = Spks;
elseif isfield(Expt.Header,'probe') % a multi channel file
    spkfile = [];
end


if isempty(Spks) & isempty(AllSpikes)
load(spkfile);
Spks = CleanSpikes(Ch3);
Spks.times = Spks.times .* 10000;
end
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
DATA.plot.SpikeMaxV = 5;
DATA.plot.SpikeVsep = 2;
DATA.plot.dvdt = 0;
DATA.plot.showwave = 0;
DATA.probe = 1;
DATA.currenttrial = 1;
DATA.plot.showsync = 0;
DATA.state.fixrange = 0;
DATA.state.preperiod = preperiod;
DATA.state.postperiod = postperiod;
DATA.syncsign = 3;
DATA.state.recut = 1;
DATA.cluster = [];
DATA.currentcluster = 1;
DATA.firsttrial = 1;
DATA.clusterArange = [7:11];
DATA.clusterBrange = [12:18];
DATA.clusterErange = [1:40];
DATA.plot.nodc = 1;
DATA.plot.voffsets = 0;
    colors = mycolors;
    DATA.spkcolor{1} = [0.5 0.5 0.5];
    DATA.spkcolor(2:20) = colors(1:19);
    DATA.ptsize = 1;
DATA.probelist = DATA.probe;
DATA.plot.autoscale = 1;
DATA.AllData.Trialids = [Expt.Trials.Trial];
DATA.Expts{1} = Expt;
DATA.explabels{1} = Expt.Header.Name;
DATA.currentexpt = 1;
DATA.xprobes = [];
DATA.plot.showwave = 0;
DATA.vstep = 1;
DATA.plot.syncoverlay = 0;
DATA.densityplot = 0;
DATA.plot.setptsize = 0;
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
