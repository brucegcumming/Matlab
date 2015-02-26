function [DATA, counts] = CountSpikes(DATA, expid, varargin)
%CountSpikes(DATA, expid, ...  expid is index id DATA.exptnos
%Counts spikes in data files and puts them into correct trials of Expts

if isfield(DATA,'AllClusters')
    %can be using AllData.Spikes, but have no spikes
elseif ~isempty(DATA.AllData.Spikes) && length(DATA.AllData.Spikes.codes) == 0
    counts = 0;
    return;
end
replot = 0;
reclassify = DATA.state.recount;
j = 1;
while j <= nargin -2
    if strncmp(varargin{j},'reclassify',3)
        reclassify = 1;
    elseif strncmp(varargin{j},'replot',3)
        replot = 1;
    end
    j = j+1;
end

nt = 1;
if DATA.state.recut
    ctype = 2;
else
    ctype = 1;
end
if  DATA.state.online == 2
    return;
end

%This is all done below anyway.

if 0 && DATA.Expts{expid}.gui.classified == 0  || reclassify
    DATA = SetExptSpikes(DATA, expid, 0,'quick');
end

xcl = [];

spikelist = cmb.WhichClusters(DATA.toplevel);

if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
elseif isfield(DATA,'AllClusters')
    if iscell(DATA.AllClusters)
        if DATA.listbycell && isfield(DATA.CellDetails,'exptids')
            eid = find(DATA.CellDetails.exptids == DATA.exptnos(expid));
            if DATA.listbycell == 2
                pid = DATA.probe;
                if isfield(DATA,'muclusters')
                    id = find(DATA.expid == eid); %this was wrong. jus use eid
                    cid = DATA.muclusters(eid);
                else
                    cid = 1;
                end
            else
                [pid, cid] = find(squeeze(DATA.CellList(eid(1),:,:)) == DATA.probe);
            end
        else
            pid = DATA.probe;
            cid = cmb.WhichClusters(DATA.toplevel);
        end
        [DATA, D] = cmb.CheckClusterLoaded(DATA, expid, pid);
        if length(pid) == 1
            Spks = DATA.AllClusters{expid}(pid(1));
            spikelist = cid;
            if isfield(Spks,'excludetrialids')
                if cid(end) <= length(Spks.excludetrialids) && cid(1) > 0
                    xcl = Spks.excludetrialids{cid};
                end
            end
        elseif length(pid) > 1
            Spks = DATA.AllClusters{expid}(pid(1));
            spikelist = cid(1);
            if cid(1) <= length(Spks.excludetrialids)
                xcl = Spks.excludetrialids{cid(1)};
            end
        else
            Spks = [];
        end
        if isfield(Spks,'missingtrials') && ~isempty(Spks.missingtrials)
            xcl = [xcl Spks.missingtrials];
        end
        if isempty(pid)
            counts = [];
            return;
        end
        newcode = 0;
        if newcode %causes trouble with recursion
            [DATA, counts] = cmb.CountSpikes(DATA, expid, varargin{:});
            [DATA.Expts{expid}.Trials.excluded] = deal(0);
            [id, aid] = intersect([DATA.Expts{expid}.Trials.id],xcl);
            [DATA.Expts{expid}.Trials(aid).excluded] = deal(1);
            return;
        end
    else
        Spks = DATA.AllClusters(DATA.probe);
        xcl = [];
    end
elseif DATA.state.online
    [havecodes, spkcache] = cmb.SpkCache(DATA,expid, DATA.probe,'check');

    if expid == DATA.AllData.Spikes.exptid
        Spks = DATA.AllData.Spikes;
%if the Current AllData.Spikes matches the stored list, but has not codes set, its
%probably been reloaded, so add back the old codes
        if sum(Spks.codes(:,2) == 0) && havecodes && size(Spks.codes,1) == length(spkcache.codes)           DATA.AllData.Spikes.codes(:,2)  = spkcache.codes;
           Spks = DATA.AllData.Spikes;
           if havecodes == 1
               fprintf(' applying previous cut');
           elseif havecodes == 2
                fprintf('No events in range');
           end
        else
            needspikes =1;
        end
    elseif havecodes 
            fprintf(' Using previous cut');            
            Spks.times = spkcache.times;
            Spks.codes = spkcache.codes;
    else
        counts = 0;
        return;
    end
else
    Spks = DATA.AllData.Spikes;
end

if isempty(Spks) | ~isfield(Spks,'codes')
    for nt = 1:length(DATA.Expts{expid}.Trials)
        DATA.Expts{expid}.Trials(nt).count = 0;
        DATA.Expts{expid}.Trials(nt).Spikes = [];
    end
    
    counts = 0;
    return;
end
alli = [];
%restrict finds to the range of this expt - speeds things up
espk = ExptSpikeListAll(DATA,expid,Spks.times);
nt = 1;
if size(Spks.codes,2) == 1
    ctype = 1;
end
%ispk = find(ismember(Spks.codes(espk,ctype),spikelist));
%espk = espk(ispk);
    if ismember(DATA.Expts{expid}.Trials(nt).id,xcl)
        DATA.Expts{expid}.Trials(nt).excluded = 1;
    else
        DATA.Expts{expid}.Trials(nt).excluded = 0;
    end
    
    if isempty(Spks.codes)
        counts = [];
        return;
    end
excluded = ismember([DATA.Expts{expid}.Trials.id],xcl);
iscl = ismember(Spks.codes(espk,ctype),spikelist);
cellspks = espk(iscl); %cell events that are in this expt
muspks = espk(~iscl); %mu events that are in this expt
cellt = Spks.times(cellspks);
mut = Spks.times(muspks);

for trial = [DATA.Expts{expid}.Trials]
    if isempty(Spks.times)
        DATA.Expts{expid}.Trials(nt).Spikes = [];
        DATA.Expts{expid}.Trials(nt).count = 0;
        DATA.Expts{expid}.Trials(nt).OSpikes = [];
    else
        ispk = find(cellt > trial.Start(1)-DATA.state.preperiod & cellt < trial.End(end)+DATA.state.postperiod);
        spks = double(round(cellt(ispk)) - trial.Start(1));
        DATA.Expts{expid}.Trials(nt).Spikes = spks;
        %    DATA.Expts{expid}.Trials(nt).Spikes = double(round(Spks.times(ispk) - trial.Start(1)));
        DATA.Expts{expid}.Trials(nt).count = sum(spks > 500);
        muispk = find(mut > trial.Start(1)-DATA.state.preperiod & mut < trial.End(end)+DATA.state.postperiod);
        DATA.Expts{expid}.Trials(nt).OSpikes = double(round(mut(muispk)) - trial.Start(1));
        alli = [alli ispk(:)' muispk(:)'];
        DATA.Expts{expid}.Trials(nt).Ocodes = Spks.codes(muspks(muispk),ctype);
    end
    DATA.Expts{expid}.Trials(nt).excluded = excluded(nt);
    nt = nt+1;
end
if replot
    cmb.SetFigure(DATA.tag.dataplot,DATA);
    args = cmb.PlotArgs(DATA, DATA.Expts{expid});
    PlotExpt(DATA.Expts{expid},args{:});
    t = get(get(gca,'title'),'String');
    title([t sprintf(' P%d',DATA.probe) sprintf(' Cl %d',spikelist)]);
end
DATA.Expts{expid}.gui.counted = sum(spikelist);
if ~isfield(DATA.Expts{expid}.gui,'setispk') || length(alli) ~= DATA.Expts{expid}.gui.setispk
    DATA.Expts{expid}.gui.spks = alli;
    DATA.Expts{expid}.gui.setispk = 2;
end
DATA.spikelist = spikelist;

if isempty(Spks.codes)
    nc = 0;
    counts(1) = 0;
else
    nc = max(Spks.codes(alli,ctype));
    counts(1) = sum(Spks.codes(alli,ctype) == 0);
end
if nc
    for j = 1:nc
        counts(j+1) = sum(Spks.codes(alli,ctype) == j);
    end
end
DATA.Expts{expid}.gui.spkcounts = counts;
DATA.Expts{expid}.gui.nclusters = nc;

