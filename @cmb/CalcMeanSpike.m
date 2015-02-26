function DATA = CalcMeanSpike(DATA,expid)

p = DATA.probe;
Spks = cmb.GetSpikeData(DATA, p);

if ~isfield(DATA,'TemplateScores')
DATA.TemplateScores = [];
end
if ~isfield(Spks, 'values') %Can happen if bysuffix
return;
end
if size(Spks.values,2) == size(DATA.Templates,2)
for j = 1:size(DATA.Templates,1)
scores(j,:) = Spks.values * DATA.Templates(j,:)';
end
else
scores = [];
DATA.TemplateScores = [];
end
if size(scores,1) > size(DATA.TemplateScores,2) %% new templates
DATA.TemplateScores(:,size(scores,1),:) = 0;
end
empties = [];  %trial with no spikes, not missing trials
for eid = expid;
[DATA, ispk] = SetExptSpikes(DATA, eid,0);
if length(ispk) > 100
if length(ispk) > 10000
xid = find(DATA.Spikes.cx(ispk) > prctile(DATA.Spikes.cx(ispk),99.9));
else
v = sort(DATA.Spikes.cx(ispk),'descend');
xid = find(DATA.Spikes.cx(ispk) > v(100));
end
DATA.MeanSpike.v(eid,p,:) = mean(Spks.values(ispk(xid),:));
DATA.MeanSpike.sd(eid,p,:) = std(Spks.values(ispk(xid),:));
end
trials = [DATA.Expts{eid}.Trials.Trial];
for j = 1:length(DATA.Expts{eid}.Trials)
times = [DATA.Expts{eid}.Trials(j).Start(1) DATA.Expts{eid}.Trials(j).End(end)];
tspk = find(Spks.times(ispk) > times(1) & ...
Spks.times(ispk) < times(2));
cspk = find(Spks.codes(ispk(tspk),2) > 0);
if isempty(cspk)
sv = sort(DATA.Spikes.cx(tspk));
if length(sv) > 5
DATA.TrialVar.cx(p,trials(j)) = -mean(sv(end-5:end));
else
DATA.TrialVar.cx(p,trials(j)) = -mean(DATA.Spikes.cx(tspk));
end
else
cspk = ispk(tspk(cspk));
DATA.TrialVar.cx(p,trials(j)) = mean(DATA.Spikes.cx(cspk));
end
if length(tspk) && ~isempty(scores)
DATA.TemplateScores(p,:,trials(j)) = max(scores(:,ispk(tspk)),[],2);
else
empties = [empties, trials(j)];
end
end
trigt = 9;
if isfield(DATA.Expts{eid},'Cluster') & cmb.iscluster(DATA.Expts{eid}.Cluster,1,p) ...
& isfield(DATA.Expts{eid}.Cluster{1,p},'dprime')
DATA.MeanSpike.dprimes(p,eid) = DATA.Expts{eid}.Cluster{1,p}.dprime;
spks = find(Spks.codes(ispk,2) == 1);
%only average the group that triggered in the same direction
sgns = sign(Spks.values(ispk(spks),trigt));
if mean(sgns) > 0
cspk = ispk(spks(find(sgns > 0)));
else
cspk = ispk(spks(find(sgns < 0)));
end
DATA.MeanSpike.Cluster(eid,p,:) = mean(Spks.values(cspk,:));
DATA.MeanSpike.ClusterSD(eid,p,:) = std(Spks.values(cspk,:));
if isfield(DATA.Expts{eid}.Cluster{1,p},'autocut')
DATA.MeanSpike.autocut(p,eid) = DATA.Expts{eid}.Cluster{1,p}.autocut;
else
DATA.MeanSpike.autocut(p,eid) = 0;
end
else
DATA = cmb.AutoCut(DATA, eid, 1,'noplot');
%codes have changed in DATA, so reset Spks
Spks = cmb.GetSpikeData(DATA, p);
spks = find(Spks.codes(ispk,2) == 1);
sgns = sign(Spks.values(ispk(spks),trigt));
if mean(sgns) > 0
cspk = ispk(spks(find(sgns > 0)));
else
cspk = ispk(spks(find(sgns < 0)));
end
DATA.MeanSpike.Cluster(eid,p,:) = mean(Spks.values(cspk,:));
DATA.MeanSpike.ClusterSD(eid,p,:) = std(Spks.values(cspk,:));
if length(ispk) < 10
DATA.MeanSpike.dprimes(p,eid) = NaN;
else
DATA.MeanSpike.dprimes(p,eid) = DATA.Expts{eid}.Cluster{1,p}.dprime;
end
DATA.MeanSpike.autocut(p,eid) = 1;
end
DATA.MeanSpike.nspk(p,eid) = length(cspk);
if DATA.MeanSpike.Cluster(eid,p,trigt) > 0
nspk = find(Spks.values(ispk,trigt) > 0 & Spks.codes(ispk,2) ~= 1);
else
nspk = find(Spks.values(ispk,trigt) < 0 & Spks.codes(ispk,2) ~= 1);
end
DATA.MeanSpike.NotCluster(eid,p,:) = mean(Spks.values(ispk(nspk),:));
end
id = find(DATA.TemplateScores(p,end,:) == 0);

