function PlotTemplateScores(DATA,ti,varargin)

trials = [1:size(DATA.TemplateScores,3)]';
trialrange = [1 length(trials)];
j = 1;

while j <= length(varargin)
if strncmpi(varargin{j},'Trials',5)
j = j+1;
trialrange = varargin{j};
end
j = j+1;
end
if isfield(DATA.Comments,'Peninfo')
probesep = DATA.Comments.Peninfo.probesep;
else
probesep = 75;
end
nprobes = size(DATA.TemplateScores,1);
id = find(squeeze(sum(DATA.TemplateScores(:,ti,:),1)) > 0 & trials > trialrange(1));
imagesc(trials(id),[1 nprobes],squeeze(DATA.TemplateScores(:,ti,id)));
Trials = [];
eds = []
Expts = []; nx = 1;
for ex = 1:length(DATA.Expts)
[a,b] = ismember(DATA.Expts{ex}.Trials(1).Trial,id);
if a
Expts(nx).Start = id(b); %id # of starting trial
nx = nx+1;
end
Trials = [Trials [DATA.Expts{ex}.Trials.Trial]];
if isfield(DATA.Expts{ex}.Trials,'ed')
eds([DATA.Expts{ex}.Trials.Trial]) = [DATA.Expts{ex}.Trials.ed];
else
eds([DATA.Expts{ex}.Trials.Trial]) = DATA.Expts{ex}.Stimvals.ed;
end
end
hold on;
if size(DATA.TemplateInfo,1) >= ti
estart = GetEval(DATA.Expts{DATA.TemplateInfo(ti).exid},'ed');
ed = DATA.TemplateInfo(ti).probe + (eds(id)-estart) .* 1000./probesep;
plot(trials(id),ed,'w');
for j = 1:length(Expts)
plot([Expts(j).Start Expts(j).Start],[1 nprobes],':');
end
end


