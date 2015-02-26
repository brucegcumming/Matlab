function TestCombine(DATA, dlgs)
%try different ways to count spikes/load expts for spead


function LocalCombine(DATA)

spikelist = [];
ids = [];
j = 1;
while  j <= length(varargin)
    if strncmpi(varargin{j},'cluster',5)
        j = j+1;
        spikelist = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        j = j+1;
        ids = varargin{j};
    end
    j = j+1;
end
if isempty(ids)
    ids = get(DATA.elst,'value');
end
if isempty(spikelist)
     spikelist = DATA.spikelist;
end

suffixes = [];
plotres = [];
needfields = {};
oldf = gcf;
ts = now;
if DATA.state.online < 2
    [DATA, Stimvals, needfields, ok] = cmb.CheckCombine(DATA, dlgs);
    if ok < 0
        Expt = [];
        return;
    end
    figure(gcf);
end
if DATA.extype ==1
    msgbox('You must select an Expt type before combining','Selection Error');
    Expt = [];
    plotres = [];
    return;
end

if DATA.listbycell
    celleid = [];
    if DATA.listbycell == 2
        %cid is probe number, which might vary if adjusting with drift
        [eid, cid] = cmb.FindNoCell(DATA, DATA.probe, DATA.expid); 
        xid = find(DATA.muclusters > 1); %this only works with no drift adjustment
        if ~isempty(xid)
            alleid = union(eid, xid);
            ia = ismember(alleid,eid);
            ib = ismember(alleid,xid);
            allcid(ia) = 1; %cluster id
            allcid(ib) = DATA.muclusters(xid);
            eid = alleid; 
            clid = allcid;
        end
        if isfield(DATA,'CellList')
            [celleid, cellclid] = find(squeeze(DATA.CellList(:,DATA.probe,:)) > 0);
            if ~isempty(DATA.muCellList) %May have some cells marked, both others marked as mu. 
                [xeid, xcclid] = find(squeeze(DATA.muCellList(:,DATA.probe,:)) > 0);
                [celleid, ia] = setdiff(celleid,xeid); %expts we want mu for
                cellclid = cellclid(ia);
            end
            [a, ida, idb] = intersect(DATA.exptnos(DATA.expid(ids)), DATA.CellDetails.exptids(celleid));
            celleid = celleid(idb); %expts without cell
            [a, ida, idb] = intersect(DATA.exptnos(DATA.expid(ids)), DATA.CellDetails.exptids(eid));
            ids = ids(ida);
            cellprobes = cid(idb);
        else
            cellprobes = ones(size(ids)).* DATA.probe;
            ida = DATA.expid(ids);
        end
    else
        [eid, cid, clid] = cmb.FindCell(DATA, DATA.probe);
        [a, ida, idb] = intersect(DATA.exptnos(DATA.expid(ids)), DATA.CellDetails.exptids(eid));
        ids = ids(ida);
        cellprobes = cid(idb);
    end
    probepos = mean(cellprobes);
    if isempty(ida)
        Expt = [];
        return;
    end
    figure(gcf);
else
    cellprobes = ones(size(ids)).* DATA.probe;
end

if DATA.listbycell == 2
    for j = 1:length(ids)
        e = DATA.expid(ids(j));
        if isempty(DATA.AllClusters{e}(cid(j)).codes)
            meanrate(j) = NaN;
            maxrate(j) = NaN;
        else
            meanrate(j) = sum(DATA.AllClusters{e}(cid(j)).codes(:,1) == 1)./sum([DATA.Expts{e}.Trials.dur]);
            maxrate(j) = length(DATA.AllClusters{e}(cid(j)).codes)./sum([DATA.Expts{e}.Trials.dur]);
        end
    end
    meanrate = meanrate .*10000;
    maxrate = maxrate .*10000;
    tgtrate = nanmean(meanrate);
    if ~isempty(celleid)
        for j = 1:length(celleid)
            e = DATA.CellDetails.exptids(celleid(j));
            cmeanrate(j) = sum(DATA.AllClusters{e}(cid(j)).codes(:,1) == cellclid(j))./sum([DATA.Expts{e}.Trials.dur]);
            cellid(j) = DATA.CellList(celleid(j),DATA.probe,cellclid(j));
        end
        cmeanrate = cmeanrate .*10000;
        tgtrate = nanmean(cmeanrate);
    end

    for j = 1:length(ids)
        if isnan(meanrate(j))
        else
            e = DATA.expid(ids(j));
            p = 100 .*tgtrate./maxrate(j);
            if p > 100
                p = 100;
            end
            th = prctile(DATA.AllClusters{e}(cid(j)).cx,100-p);
            id = find(DATA.AllClusters{e}(cid(j)).cx > th);
            DATA.AllClusters{e}(cid(j)).codes(id,2) = 1;
            id = find(DATA.AllClusters{e}(cid(j)).cx <= th);
            DATA.AllClusters{e}(cid(j)).codes(id,2) = 0;
        end
    end
end
%
% when combining expts with certain differences, make sure these are
% recorded for each trial, so that they can be treared differently
% add to this list any parameters shown manually with DATA.show
CheckF = {'sx' 'Fr' 'rb' 'ap' 'jx' 'ns' 'sf' 'or' 'backxo' 'bo' 'backyo' 'co' 'Bc' 'bh' 'SpikeGain' 'dd' 'puA' 'USd' 'USp' 'puF' 'dx'};
f = fields(DATA.show);
for j = 1:length(f)
    if isempty(strmatch(f{j},{'ed'},'exact')) && DATA.show.(f{j})
        CheckF = {CheckF{:} f{j}};
    end
end
for f = 1:length(CheckF)
    vals = [];
    for j = ids(1:end)
        if isfield(DATA.Expts{DATA.expid(j)}.Stimvals,CheckF{f})
            vals(j) = DATA.Expts{DATA.expid(j)}.Stimvals.(CheckF{f});
        else
            vals(j) = Inf;
        end
    end
    if length(unique(vals(ids))) > 1
        for j = ids(1:end)
            if ~isfield(DATA.Expts{DATA.expid(j)}.Trials,CheckF{f})
                DATA.Expts{DATA.expid(j)}= FillTrials(DATA.Expts{DATA.expid(j)},CheckF{f});
            end
        end
    end
end
tfields = {};
Frs = [];
Comments = [];
for j = ids(1:end)
    tfields = union(tfields,fields(DATA.Expts{DATA.expid(j)}.Trials));
    if isfield(DATA.Expts{DATA.expid(j)}.Stimvals,'Fr')
        Frs(j) = DATA.Expts{DATA.expid(j)}.Stimvals.Fr;
    else
        Frs(j) = 0;
    end
    Comments = [Comments  DATA.Expts{DATA.expid(j)}.Comments];
end
tfields = union(tfields,needfields);
if length(unique(Frs(ids))) > 1
    for j = ids(1:end)
        if ~isfield(DATA.Expts{DATA.expid(j)}.Trials,'Fr')
            DATA.Expts{DATA.expid(j)}= FillTrials(DATA.Expts{DATA.expid(j)},'Fr');
        end
    end
end
j = DATA.expid(ids(end));
DATA.combineids = DATA.expid(ids);
if length(DATA.extype) == 1
    DATA.allcombineids{DATA.extype} = DATA.combineids;
end
Expt.Header = DATA.Expts{j}.Header;
Trials = [];
Pulses = [];
nb = 1;
Clusters = {};
if DATA.plot.autoclustermode < 3
    DATA.autocut = 1;
else
    DATA.autocut = 0;
end
counts = {}; nspk = 0;
excludeids = [];
MeanPulses = {};
for ij = 1:length(ids);
    j = ids(ij);
    eid = DATA.expid(j);
    E = DATA.Expts{eid};
    BlockStart(nb) = E.Trials(1).Trial;
    BlockStartid(nb) = E.Trials(1).id;
    if isfield(E.Header,'suffix')
        suffixes(nb) = E.Header.suffix;
    end
    exped(nb) = GetEval(DATA.Expts{DATA.expid(j)},'ed');
    Header = DATA.Expts{DATA.expid(j)}.Header;
    if isfield(DATA,'AllClusters')
        Clusters{nb} = cmb.SmallCluster(DATA.AllClusters,DATA.expid(j),cellprobes(ij));
        if isfield(Clusters{nb},'excludetrialids')
            excludeids = cat(1,excludeids(:),Clusters{nb}.excludetrialids{1}(:));
        end
        if isfield(Clusters{nb},'missingtrials')
            excludeids = cat(1,excludeids(:),Clusters{nb}.missingtrials(:));
        end
    elseif isfield(DATA.Expts{DATA.expid(j)},'Cluster')
        Clusters{nb} = DATA.Expts{DATA.expid(j)}.Cluster;
        recount = 0;
    end
    if isfield(DATA.Expts{DATA.expid(j)},'MeanPulse')
        MeanPulses{nb} = DATA.Expts{DATA.expid(j)}.MeanPulse;
        PulseSDs{nb} = DATA.Expts{DATA.expid(j)}.MeanPulse;
    end
    if DATA.state.nospikes
        recount = 1; %% for chagning probes
    elseif DATA.autocut & (size(Clusters) < nb | size(Clusters{nb},2) < DATA.probe ...
            | isempty(Clusters{nb}) | isempty(Clusters{nb}{1,DATA.probe})...
            | ~isfield(Clusters{nb}{1,DATA.probe},'x')...
            | (isfield(Clusters{nb}{1,DATA.probe},'autocut') & Clusters{nb}{1,DATA.probe}.autocut > 0  &DATA.state.redoautocut))
        DATA = cmb.AutoCut(DATA,DATA.expid(j),j);
        Clusters{nb} = DATA.Expts{DATA.expid(j)}.Cluster;
        recount = 1;
    else
        recount = 0;
    end
    if DATA.state.psychonly
        recount = 0;
    elseif DATA.state.online == 0 | recount
        %            if cmb.iscluster(DATA.Expts{eid}.Cluster,DATA.currentcluster,DATA.probe)
        %                DATA.cluster = DATA.Expts{eid}.Cluster;
        %            end
        [DATA, counts{j}] = cmb.CountSpikes(DATA,DATA.expid(j)); %% recount in case cluster # changed
    end
    %tfields is the combined fields of all expts (above). Sao any differences
    %mean it is missing from the current Expt;
    newf = setdiff(tfields, fields(DATA.Expts{DATA.expid(j)}.Trials));
    for k = 1:length(newf)
        if strcmp(newf{k},'FalseStart')
            [DATA.Expts{DATA.expid(j)}.Trials.(newf{k})] = deal(0);
        else
            DATA.Expts{DATA.expid(j)} = FillTrials(DATA.Expts{DATA.expid(j)},newf{k});
        end
    end
    newtrials = DATA.Expts{DATA.expid(j)}.Trials;
    Trials = [Trials newtrials];
    if isfield(DATA.Expts{DATA.expid(j)},'ExcludeId')
    end
    if isfield(DATA.Expts{DATA.expid(j)},'Pulses')
        Pulses = [Pulses DATA.Expts{DATA.expid(j)}.Pulses];
    end
    if Header.trange(1) < Expt.Header.trange(1)
        Expt.Header.trange(1) = Header.trange(1);
    end
    if Header.trange(2) > Expt.Header.trange(2)
        Expt.Header.trange(2) = Header.trange(2);
    end
    nb = nb+1;
end
for j = 1:length(counts)
    for k = 1:length(counts{j})
        nspk(j,k) = counts{j}(k);
    end
end
eid = DATA.expid(ids(end));
Expt.Stimvals = Stimvals;
Expt.Trials = Trials;
Expt.Header.Name = cmb.BuildName(DATA.Expts{eid}.Header.Name);
Expt.Header.BlockStart = BlockStart;
Expt.Header.BlockStartid = BlockStartid;
Expt.Header.depths = exped;
Expt.Header.Clusters = Clusters;
Expt.Header.Combined = ids;
Expt.Header.Combineids = DATA.expid(ids);
Expt.Header.Spikelist = DATA.spikelist;
Expt.Header.CombineDate = now;
Expt.Header.bysuffix = DATA.bysuffix;
Expt.Header.combineversion = DATA.progversion;
if DATA.bysuffix && ~isempty(suffixes)
    Expt.Header.suffixes = suffixes;
else
    Expt.Header.suffixes = 0;
end
Expt.Header.nspk = sum(nspk);
if ~isempty(MeanPulses)
    Expt.MeanPulses = MeanPulses;
    Expt.PulseSDs = PulseSDs;
end
if isfield(DATA,'AllClusters')
    if DATA.listbycell  == 1
        Expt.Header.probe = probepos;
        Expt.Header.cellnumber = DATA.probe;
        [elst, cid, clid] = cmb.FindCell(DATA, DATA.probe);
        cellprobe(elst) = cid;
        cellclust(elst) = clid;
        Expt.Header.muforcell = 0;
    elseif DATA.listbycell  == 2
        Expt.Header.probe = DATA.probe;
        Expt.Header.probes = cellprobes;
        Expt.Header.cellnumber = 0;
        if ~isempty(celleid)
           Expt.Header.muforcell = cellid;
        end
    else
        Expt.Header.probe = DATA.probe;
    end
    for j = 1:length(ids)
        eid = DATA.expid(ids(j));
        if DATA.listbycell ~= 1
            ndip = length(DATA.AllClusters{eid}(DATA.probe).dips);
        end
        if DATA.listbycell == 1
            probe = cellprobe(DATA.cellexid(eid));
            cid = cellclust(DATA.cellexid(eid));
            
            if cid > 1
                try
                    Expt.Header.Clusters{j}.cluster = cid;
                    dips = DATA.AllClusters{eid}(probe).next{cid-1}.dips;
                    Expt.Header.dropi(j) = DATA.AllClusters{eid}(probe).next{cid-1}.dropi;
                catch err
                    fprintf('Missing next or field for E%dP%dC%d(%s)\n',eid,probe,cid,err.message);
                    dips = DATA.AllClusters{eid}(probe).dips;
                    Expt.Header.dropi(j) = NaN;
                end
            else
                dips = DATA.AllClusters{eid}(probe).dips;
                Expt.Header.dropi(j) = DATA.AllClusters{eid}(probe).dropi;
                Expt.Header.Clusters{j}.cluster = 1;
            end
            if length(dips) < 4
                cprintf('errors','Missing Quantification E%dP%dC%d(%d values)\n',eid,probe,cid,length(dips));
                dips(4) = NaN;
            end
            Expt.Header.dips(j,1:length(dips)) = dips;
            
        elseif ndip > 0
            Expt.Header.dips(j,1:ndip) = DATA.AllClusters{eid}(DATA.probe).dips;
            Expt.Header.dropi(j) = DATA.AllClusters{eid}(DATA.probe).dropi;
        end
    end
end

if ~isempty(Comments)
    Expt.Comments = Comments;
end
if isfield(DATA,'Comments') & isfield(DATA.Comments,'Peninfo')
    Expt.Header.Peninfo = DATA.Comments.Peninfo;
    id = strfind(DATA.Comments.Peninfo.trode,'Contact');
    if length(id)
        x = id(1);
        id = strfind(DATA.Comments.Peninfo.trode(id:end),' ');
        x = sscanf(DATA.Comments.Peninfo.trode(id+x:end),'%d');
        Expt.Header.probesep = x;
    end
end

if ~DATA.state.online & ~isfield(DATA,'AllClusters') & ~DATA.state.psychonly
    Expt.Header.SpkStats = cmb.GetSpkStats(DATA);
end
if isfield(Expt.Trials,'FalseStart')
    id = find([Expt.Trials.FalseStart] > 2);
    if ~isempty(id)
        for j = 1:length(id)
            Expt.Trials(id(j)).Start = Expt.Trials(id(j)).Start - Expt.Trials(id(j)).FalseStart;
        end
        msgbox(sprintf('%s,%s: %d Trials with long delays',Expt.Header.Name,Expt.Header.expname,length(id)));
    end
    id = find([Expt.Trials.FalseStart] ==1);
    if ~isempty(id) && isfield(Expt.Header,'AplayVer')
        msgbox(sprintf('%s,%s: %d Trials missing frametime',Expt.Header.Name,Expt.Header.expname,length(id)));
    end
end
if ~isempty(Pulses)
    Expt.Pulses = Pulses;
end
if ~isfield(Expt.Header,'Spike2Version')
    Expt.Header.Spike2Version = 1.0;
end
[a, Expt.Header.fileprefix] = cmb.CombinedName(DATA, DATA.extype, 1);
Expt.Header.nprobes = length(DATA.probelist);
%excludeids is list of excluded trials in celllist (from Plotclsuters)
%Expt.trials.excluded has excluded trias fomr AllVPcs in AllClusters
if isfield(Expt.Trials,'excluded')
    xid = find([Expt.Trials.excluded] ==1);
    Expt.Header.excludeids = union(excludeids,[Expt.Trials(xid).id]);
end
id = find(~ismember([Expt.Trials.id],excludeids));
Expt.Trials = Expt.Trials(id);

if strmatch(DATA.exptypelist{DATA.extype(1)},'TwoCylDispXhxXbhPTWO','exact')
    id = find([Expt.Trials.bh] > 0.01);
    Expt.Trials = Expt.Trials(id);
    Expt.Stimvals.e3 = 'e0';
elseif strmatch(DATA.exptypelist{DATA.extype(1)},'TwoCylDispXhxXbhPDID','exact')
    id = find([Expt.Trials.bh] > 0.01);
    Expt.Trials = Expt.Trials(id);
    Expt.Stimvals.e3 = 'e0';
    Expt.Stimvals.e2 = 'Id';
    Expt.Stimvals.et = 'dx';
end

if DATA.state.uselfp & length(ids) > 0  %reload LFP to match lengths etc
    Expt = LoadSpike2LFP(Expt,'reload');
end
id = strmatch(DATA.Expts{eid}.Header.expname,DATA.expstrs,'exact');
if Stimvals.BlockedStim
    DATA.outname = cmb.Expt2FileName(DATA, Expt, DATA.spikelist(1));
    set(DATA.saveitem,'string',DATA.outname);
end
fprintf('Combine took %.1f sec\n',mytoc(ts));
[Expt, plotres] = cmb.PlotCombined(DATA, Expt);
if DATA.state.autofit
    fit = FitExpt(plotres,'plotfit','vonmises');
    plotres.fit = fit;
    Expt.fit = fit;
    DATA = cmb.AddFitToData(DATA, plotres,fit);
    if strcmp(plotres.type{1},'stimxy')
        popRF('expts',Expt);
    end
end

