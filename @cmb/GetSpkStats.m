function spkstats = GetSpkStats(DATA)

% summary stats for a probe by expt, to try and track drift
spkstats = [];
ids = get(DATA.elst,'value');
exid = get(DATA.clst,'value');
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{DATA.probe};
else
    Spks = DATA.AllData.Spikes;
end
for j = 1:length(ids)
    eid = DATA.expid(ids(j));
    Expt =  DATA.Expts{eid};
    
    mins = min(Spks.values(Expt.gui.spks,:)');
    maxs = max(Spks.values(Expt.gui.spks,:)');
    if length(mins) > 10 && length(maxs) > 10
        spkstats.max(j,:) = prctile(maxs,[50 90 99]);
        spkstats.min(j,:) = prctile(maxs,[50 10 1]);
        spkstats.h(j,:) = prctile(maxs-mins,[50 90 99]);
        spkstats.k(j) = moment(maxs-mins,3);
    else
        spkstats.max(j,:) = [0 0 0];
        spkstats.min(j,:) = [0 0 0];
        spkstats.h(j,:) = [0 0 0];
        spkstats.k(j) = 0;
    end
    if isempty(Expt.gui.spks)
        spkstats.trange(j,:) = [0 0];
    else
        spkstats.trange(j,:) = [min(Spks.times(Expt.gui.spks)) ...
            max(Spks.times(Expt.gui.spks))];
        ispk = Expt.gui.spks;
        if ~isfield(DATA,'Spikes') || max(ispk) > length(DATA.Spikes.cx) || sum(DATA.Spikes.cx(ispk)) == 0
            DATA = CalcClusterVars(DATA,ispk,'expt',eid);
        end
        if length(ispk) > 10000
            xid = find(DATA.Spikes.cx(ispk) > prctile(DATA.Spikes.cx(ispk),99.9));
        elseif length(ispk) > 100
            v = sort(DATA.Spikes.cx(ispk),'descend');
            xid = find(DATA.Spikes.cx(ispk) > v(100));
        else
            xid = 1:length(ispk);
        end
        spkstats.V = mean(Spks.values(ispk(xid),:));
        xid = find(Spks.codes(ispk,2) == 1);
        if length(xid)
            spkstats.cV = mean(Spks.values(ispk(xid),:));
        end
    end
end

