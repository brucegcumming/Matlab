function [a,b] = TrialRange(DATA)
eid = get(DATA.elst,'value');
if isfield(DATA,'expid')
exps = DATA.expid(eid);
else
exps = 1;
end
if DATA.firsttrial > 1
a = DATA.Expts{exps(1)}.Trials(DATA.firsttrial).Trial;
else
a = DATA.Expts{exps(1)}.Trials(1).Trial;
end
if DATA.lasttrial >= DATA.firsttrial && DATA.lasttrial > 0
b = DATA.Expts{exps(end)}.Trials(DATA.lasttrial).Trial;
else
b = DATA.Expts{exps(end)}.Trials(end).Trial;
end


