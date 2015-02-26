function MakeCellTemplate(a,b)
DATA = GetDataFromFig(a); 
cn = 1;
eid = get(DATA.elst,'value');
exps = DATA.expid(eid);
t(1) = DATA.Expts{exps(1)}.Trials(1).Start(1);
t(2) = DATA.Expts{exps(end)}.Trials(end).End(end);
spkt = DATA.AllData.Spikes.times;
ispk = find(DATA.AllData.Spikes.times > t(1) & DATA.AllData.Spikes.times < t(2)...
& DATA.AllData.Spikes.codes(:,2) == cn);
mnspk = mean(DATA.AllData.Spikes.values(ispk,:));
subplot(2,1,1);
plot(mnspk);
subplot(2,1,2);

it = findobj(get(a,'parent'),'Tag','CellNumber');
if ~isempty(it)
cell = get(it(1),'value');
else 
cell = 1;
end
DATA.Templates(cell,:) = mnspk;
tmpl = DATA.AllData.Spikes.values * mnspk';
for j = length(DATA.Expts):-1:1
for k = length(DATA.Expts{j}.Trials):-1:1
ispk = find(spkt > DATA.Expts{j}.Trials(k).Start(1) & spkt < DATA.Expts{j}.Trials(k).End(end));
if length(ispk)
score(DATA.Expts{j}.Trials(k).Trial) = max(tmpl(ispk));
end
end
end
set(DATA.toplevel,'UserData',DATA);

%cmb.SaveCellList(DATA);
plot(score(find(score > 0)));

