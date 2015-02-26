function Expt = AddComments(Expt, varargin)

d = GetName(Expt,'folder');
cmfile = [d '/Comments.mat'];
if ~exist(cmfile,'file')
    return;
end
load(cmfile);
if iscell(Expt)
    for j = 1:length(Expt)
        Expt{j} = AddExptComment(Expt{j},Comments);
    end
else
    Expt = AddExptComment(Expt,Comments);
end

function Expt = AddExptComment(Expt, Comments)

e = GetExptNumber(Expt);
id = [];
for j = 1:length(Comments)
    if ismember(e,Comments(j).exptno)
        id(end+1) = j;
    end
end
for j = id(:)'
    Expt.Comments.text{end+1} = ['cm=' Comments(j).comment];
    Expt.Comments.times(end+1) = Expt.Trials(end).End(end)+10000;
end

[~,uid] = unique(Expt.Comments.text);
[a,b] = sort(Expt.Comments.times(uid));
uid = uid(b);
Expt.Comments.text = Expt.Comments.text(uid);
Expt.Comments.times = Expt.Comments.times(uid);