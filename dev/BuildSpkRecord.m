 function [res, details] = BuildSpikeRecord(path, expts, varargin)
Expts = {};
Clusters = {};
EM = [];
LFP = [];

showspks = [];
state.showcells = [];
state.marksaccades = 0;
state.plotmu = 0;
state.showsacsdf=0;
state.sdftype = 'switch';
state.getlfp = 0;
state.eyetrigger = 1;
j = 1;
while j <= length(varargin)
    if iscell(varargin{j}) && isfield(varargin{j}{1},'Trials')
        Expts = varargin{j};
        
    elseif iscell(varargin{j}) && isfield(varargin{j}{1},'mahal')
        Clusters = varargin{j};
    elseif isstruct(varargin{j}) && isfield(varargin{j},'Trials')
        if isfield(varargin{j}.Trials,'LFP')
            LFP = varargin{j};
        else
            EM = varargin{j};
        end
    elseif isstruct(varargin{j}) && isfield(varargin{j},'blklen')
        FullVData = varargin{j};
    elseif isnumeric(varargin{j}) && ndims(varargin{j}) == 3
        CellList = varargin{j};
    elseif strncmpi(varargin{j},'lfp',3) %make SDF for stim start with no nonimal stim change
        state.getlfp = 1;
    elseif strncmpi(varargin{j},'noswitch',4) %make SDF for stim start with no nonimal stim change
        state.sdftype = 'noswitch';
    elseif strncmpi(varargin{j},'nomu',4)
        state.plotmu = 0;
    elseif strncmpi(varargin{j},'plotmu',6)
        state.plotmu = 1;
    elseif strncmpi(varargin{j},'triglfp',5)
        state.lfptrigger = 1;
        state.eyetrigger = 0;
    elseif strncmpi(varargin{j},'sacsdf',5)
        state.showsacsdf = 1;
    elseif strncmpi(varargin{j},'showspks',5)
        j = j+1;
        showspks = varargin{j};
    elseif strncmpi(varargin{j},'showcells',5)
        j = j+1;
        state.showcells = varargin{j};
        plotmu = 1;
    end
    
    j = j+1;
end

if isstruct(path)  %plot results
    if isempty(state.showcells)
        state.showcells = 1:length(path.spktimes);
    end
    PlotBlock(path, [], state);
    return;
end
        

if length(expts) > 1
    for j = 1:length(expts)
        GetFigure(sprintf('Expt%d',expts(j)));
        [res{j}, details{j}] = BuildSpkRecord(path, expts(j), varargin{:});
    end
    return;
end
[a,b] = fileparts(path);
prefix = ['lem' b];
eid = expts(1);
res.exptno = expts;
ename = sprintf('%s/%s.%d.mat',path, prefix,eid);
eaname = sprintf('%s/%sA.%d.mat',path, prefix,eid);
emname = strrep(ename,'.mat','.em.mat');
lfpname = strrep(ename,'.mat','.lfp.mat');
cellname = sprintf('%s/CellList.mat',path );

cname = sprintf('%s/Expt%dClusterTimes.mat',path, eid);
if ~exist(ename)  && isempty(Expts)
    fprintf('Can''t find %s\n',ename);
    return;
end
if ~exist(cname) && isempty(Clusters)
    fprintf('Can''t find %s\n',cname);
    return;
end

if exist(cellname)
    load(cellname);
end

if isempty(Expts)
    [a, Expts] = APlaySpkFile(ename,'bysuffix');
end
if isempty(Clusters)
   [Clusters, FullVData] =  LoadClusters(path,eid);
end

if isempty(showspks) && isempty(state.showcells)
    showspks = 1:length(Clusters);
end
if exist(emname) && isempty(EM)
    load(emname);
    EM = CheckSaccades(Expt,'allsz');
end

Expt = Expts{1};
ibs = [Expt.Trials.IB];
new = find(abs(diff(ibs)) > 0)+1;
noswitchid = setdiff(1:length(ibs),new);
stimtime = [Expt.Trials(new).Start]./10000;
res.stimids = ibs(new);
res.starttime = [Expt.Trials(noswitchid).Start]./10000;
ns = length(stimtime);
details.Expt = Expt;
res.spktimes = {};
for j = 1:length(FullVData.blkstart)
    res.blocktimes(1,j) = FullVData.blkstart(j);
    res.blocktimes(2,j) = FullVData.blkstart(j)+FullVData.blklen(j).*FullVData.samper;
end


if exist('CellList','var')
    e = find(CellDetails.exptids == eid);
    doneprobes = zeros(1,size(CellList,2));
    for j = 1:length(state.showcells)
        [pid, cid] = find(squeeze(CellList(e,:,:)) == state.showcells(j));
        if length(pid)
        if cid == 1
            res.spktimes{j} = Clusters{pid}.times;
        else
            if ~isfield(Clusters{pid}.next{cid-1},'times')
                Clusters{pid}.next{cid-1}.times = Clusters{pid}.t(Clusters{pid}.clst == cid+1);
            end
            res.spktimes{j} = Clusters{pid}.next{cid-1}.times;
        end
        if plotmu
            row = pid;
        else
            row = j;
        end
        doneprobes(pid) = 1;
        res.suprobes(j) = pid;
        res.cellids(j) = state.showcells(j);
        end
    end
else
    res.spktimes = {};
    doneprobes = zeros(1, length(Clusters));
end
id = find(doneprobes ==0);
res.muprobes = id;
for j = 1:length(id)
    res.mutimes{j} = Clusters{id(j)}.times;
end


endtime = [Expt.Trials(new).End];
IBid = [Expt.Trials(new).Start];
%PlotTrialOnOff(stimtime,endtime);
res.stimtime = stimtime;

details.nsdf = PlotBlock(res, EM, state, LFP);


function nsdf = PlotBlock(B, EM, state, varargin)
nsdf = [];
LFP = [];


j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'Trials')
        LFP = varargin{j};
    end
    j = j+1;
end
hold off;
for j = 1:length(B.stimtime)
    plot([B.stimtime(j) B.stimtime(j)],[0 1]);
    hold on;
    text(B.stimtime(j),0.5,sprintf('%d',B.stimids(j)));
end

for j = 2:size(B.blocktimes,2)
    x(1) = B.blocktimes(2,j-1);
    x(2) = B.blocktimes(1,j);
    y = [0 1];
    plot(x,y,'r-','linewidth',3);
    y = [1 0];
    plot(x,y,'r-','linewidth',3);
end

for j = 1:length(state.showcells)
    k = find(B.cellids == state.showcells(j));
    if length(k) == 1
    plot(B.spktimes{k},B.suprobes(k),'r.');
    end
end
if state.plotmu
    for j = 1:length(B.mutimes)
        plot(B.mutimes{j},B.muprobes(j),'k.');
    end
end
if isfield(EM,'Trials')
    for j = 1:length(EM.Trials);
        t = EM.Trials(j).Start./10000 + [1:length(EM.Trials(j).Eyevals.rh)].*EM.Header.CRrates(1);
        plot(t,EM.Trials(j).Eyevals.rh);
    end
end
oldfig = gcf;


GetFigure('sdfs');
hold off;
colors = mycolors;
Trials.Trigger = B.stimtime .* 10000;
if state.showsacsdf && isfield(EM,'Trials')
    Trials.Trigger = [];
    ts = -5000:10:5000;
    for j = 1:length(EM.Trials)
        id = find([EM.Trials(j).Saccades.size] > 0.2);
        Trials.Trigger = cat(2, Trials.Trigger, EM.Trials(j).Start + [EM.Trials(j).Saccades(id).start]);
    end
else
    ts = -10000:10:5000;
    Trials.Trigger = B.stimtime.*10000;
    if strmatch(state.sdftype,'noswitch')
        Trials.Trigger = B.starttime.*10000;
    end
end


if state.plotmu
    spktimes = [B.spktimes(:) B.mutimes(:)];
else
    spktimes = B.spktimes;
end
a = minmax(Trials.Trigger);
for j = 1:length(spktimes)
    Trials.Spikes = (spktimes{j} .*10000)';
    sdf = trigsdfa(Trials,100,ts,'freetimes');
    nsdf(j,:) = sdf./mean(sdf);
    plot(ts,sdf,'color',colors{j});
    hold on;
end
details.nsdf = nsdf;
if state.eyetrigger
    GetFigure('EyeSpeed');
    hold off;
    tsmp = -1200:1200;
    for j = 1:length(EM.Trials)
        id = find(Trials.Trigger > EM.Trials(j).Start & Trials.Trigger < EM.Trials(j).End)
        if length(id)
        sid = round((Trials.Trigger(id) - EM.Trials(j).Start)./(EM.Header.CRrates(1).*10000));
        speed = sqrt(sum(diff(EM.Trials(j).EyeData)'.^2));
        tval = repmat(tsmp',1,length(sid))+ repmat(sid,length(tsmp),1);
        tval(tval < 1) = 1;
        tval(tval > length(speed)) = length(speed); 
        ms(j,:) = mean(speed(tval),2);
        end
    end
    plot(tsmp.*EM.Header.CRrates(1),mean(ms,1));
    details.eyespeed = mean(ms,1);
elseif state.lfptrigger
    tsmp = 0:200;
    for j = 1:length(LFP.Trials)
        id = find(Trials.Trigger > LFP.Trials(j).Start & Trials.Trigger < LFP.Trials(j).End);
        if length(id)
            sid = round((Trials.Trigger(id) - LFP.Trials(j).Start)./(LFP.Header.CRsamplerate.*10000));
            tval = repmat(tsmp',1,length(sid))+ repmat(sid,length(tsmp),1);
            tval(tval < 1) = 1;
            tval(tval > size(LFP.Trials(j).LFP,1)) = size(LFP.Trials(j).LFP,1);
            for p = 1:size(LFP.Trials(j).LFP,2)
                x = LFP.Trials(j).LFP(:,p);
                ms(j,p,:) = mean(x(tval),2);
            end
        end
    end
    details.triglfp = squeeze(mean(ms,1));
    imagesc(tsmp.*LFP.Header.CRsamplerate,1:size(ms,2),details.triglfp);

end
figure(oldfig);
if state.marksaccades
    for j = 1:length(Trials.Trigger)
        plot([Trials.Trigger(j) Trials.Trigger(j)]./10000,[-10 1],'g:');
    end
end


function PlotTrialOnOff(S, E)


ns = length(S);
x(1:4:ns*4) = S;
x(2:4:ns*4) = S;
y(1:4:ns*4) = 0;
y(2:4:ns*4) = 1;
x(3:4:ns*4) = E;
x(4:4:ns*4) = E;
y(3:4:ns*4) = 1;
y(4:4:ns*4) = 0;
plot(x,y);
