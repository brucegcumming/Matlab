function SpoolAllExpts(a,b, toplevel, type)
DATA = get(toplevel,'UserData');
if type == 1
skip = 10;
elseif type == 2
skip = 20;
elseif type == 3
skip = 50;
end    
h = 0;
stop = findobj(DATA.svfig,'Tag','StopSpool');

for j = 1:length(DATA.Expts)
DATA.currentexpt = j;
trials = [1:skip:length(DATA.Expts{j}.Trials) length(DATA.Expts{j}.Trials)];
for k = trials
cmb.PlayOneTrial(DATA,k,0,'nosetgui');
if h == 0
h = text(10,0,'trial','color','r');
else
set(h,'string',sprintf('%d',DATA.Expts{j}.Trials(k).Trial));
end
end
if get(stop,'value') > 0 
break;
end
end
delete(h);

