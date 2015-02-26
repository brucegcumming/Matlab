function [id, DATA] = ReTrigger(DATA, varargin)
%Apply new Trigger and Rebuild AllV
maxrate = NaN;  %will need to get these from varargin
maxspikeset = NaN;

if DATA.spkrate > 0
    DATA.setnspk = 0;
    DATA.Trigger = 0;
end
smv = AllV.mygetappdata(DATA,'TriggerV');
Vall = AllV.mygetappdata(DATA,'Vall');



[id,  DATA.Trigger, D] = AllV.TriggerV(DATA, smv);


if length(id) > 1e6
    if DATA.interactive  < 0
        if ~isfield(DATA,'chspk')
            DATA.chspk = 0;
        end
        fprintf('E%dP%d> 1 million events. (Trigger %f) Exiting AllVPcs\n',DATA.exptno,DATA.chspk(1),DATA.Trigger);
        return;
    end
    yn =  questdlg('>1Million events. Proceed?','trigger','Yes','No','No');
    if strcmpi(yn,'no')
        return;
    end
end
DATA.cluster.eventrate = D.nevents./length(smv);
DATA.cluster.trigparams = CopyFields(DATA.cluster.trigparams,D,{'skew' 'roc' 'qq'});

if DATA.spkrate <= 0 %not set. Match previous trigger rate
    DATA.spkrate = DATA.cluster.eventrate./Vall.samper;
end

%remove any spikes at very beginning or end where there isn't enough
%data to include the whole spike
if size(id,1) > 1
    id = id';
end

if DATA.verbose > 0 fprintf('E%dP%d %dEvents (%.1fHz)\n',DATA.exptno,AllV.ProbeNumber(DATA),length(id),length(id)./DATA.duration); end
ignoreid = [];
missedtrials = [];
if isfield(DATA, 'Expt') && isfield(DATA.Expt,'Trials') && DATA.usetrials
    if isfield(DATA.Expt.Header,'expname')
        res.expname = DATA.Expt.Header.expname;
    end
    iid = [];
    piid = [];
    if isfield(DATA.cluster,'excludetrialids')
        xct = DATA.cluster.excludetrialids;
    else
        xct = [];
    end
    missedtrials = [];
    blkend = (Vall.blkstart + Vall.blklen.*Vall.samper).*10000;
    blkstart = Vall.blkstart .* 10000;
    
    for j = 1:length(DATA.Expt.Trials)
        tid = find(blkstart < DATA.Expt.Trials(j).Start(1) & blkend > DATA.Expt.Trials(j).End(end));
        if isempty(tid)
            missedtrials = [missedtrials DATA.Expt.Trials(j).id];
        end
        oid = find(Vall.t(id) > DATA.Expt.Trials(j).Start(1)./10000 - DATA.preperiod & ...
            Vall.t(id) < DATA.Expt.Trials(j).End(end)/10000 + DATA.postperiod);
        pid = find(Vall.t(id(oid)) > DATA.Expt.Trials(j).End(end)/10000); %in postperiod
        iid = [iid oid];
        piid = [piid oid(pid)];
    end
    ntrials = j;
    uid = unique(iid);
    if length(DATA.clst) == length(id)
        DATA.clst = DATA.clst(uid);
        if recluster == 4 && length(DATA.xy{1}) >= max(uid)
            DATA.xy{1} = DATA.xy{1}(uid,:);
        end
    end
    ignoreid = setdiff(id,id(uid));
    id = id(uid);
    [a, pid] = ismember(piid,uid);
    DATA.postevents = pid; %index to event index, not FullV times
    if isfield(DATA.Expt.Header,'trialdur')
        DATA.duration = (DATA.Expt.Header.trialdur+ntrials*(DATA.preperiod+DATA.postperiod))./10000;
    else
        disp(sprintf('Header Missing trialdur E%dP%d\n',DATA.exptno,DATA.probe(1)));
        DATA.duration = ntrials * (2+DATA.preperiod+DATA.postperiod);
    end
    setnspk = DATA.duration * DATA.spkrate;
else
    DATA.postevents = [];
end
DATA.missedtrials = missedtrials;
if length(missedtrials)
    DATA = AllV.AddErr(DATA,'Missing %d/%dTrials%s\n',length(missedtrials),ntrials,sprintf(' %d',missedtrials));
end

if maxrate > 0
    maxspikeset = DATA.duration .* maxrate;
end


if isempty(id)
    DATA = AllV.AddErr(DATA,'No Spikes in Expt\n');
    res.t = [];
    res.Trigger = DATA.Trigger;
    return;
else
    fprintf('%d events in Trials',length(id));
end

if isfield(D,'trigid') && length(D.trigid) >= max(uid)
    DATA.rV = smv(:,D.trigid(uid));
elseif isfield(D,'triggerV') %loaded from file
    DATA.rV = D.triggerV;
else
    DATA.rV = smv(:,id);
end