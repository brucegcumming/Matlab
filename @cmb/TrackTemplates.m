function  TrackTemplates(DATA)

DATA.TemplateScores = [];
probelist = DATA.probelist;
for p = 1:length(probelist)
fprintf('Probe %d..',p);
Spks = cmb.GetProbeFiles(DATA,p,DATA.subprobe);
scores = [];
for j = 1:size(DATA.Templates,1)
scores(j,:) = Spks.values * DATA.Templates(j,:)';
end
for eid = DATA.exabsid;
trials = [DATA.Expts{eid}.Trials.Trial];
espk = find(Spks.times> DATA.Expts{eid}.Trials(1).Start(1) & ...
Spks.times < DATA.Expts{eid}.Trials(end).End(end));
for j = 1:length(DATA.Expts{eid}.Trials)
times = [DATA.Expts{eid}.Trials(j).Start(1) DATA.Expts{eid}.Trials(j).End(end)];
tspk = find(Spks.times(espk) > times(1) & ...
Spks.times(espk) < times(2));
if length(tspk) && ~isempty(scores)
DATA.TemplateScores(p,:,trials(j)) = max(scores(:,espk(tspk)),[],2);
end
end
end
end
GetFigure('TemplateScores');
set(DATA.toplevel,'UserData',DATA);


