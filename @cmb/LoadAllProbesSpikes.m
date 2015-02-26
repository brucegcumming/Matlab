function DATA = LoadAllProbesSpikes(DATA)


set(DATA.toplevel,'Name','Loading Spikes');
eid = DATA.currentexpt(1);
for k = 1:length(DATA.currentexpt)
eid = DATA.currentexpt(k);
for j = 1:length(DATA.probelist)
DATA.probe = DATA.probelist(j);
DATA = cmb.LoadSpikes(DATA, eid);
DATA.AllData.Spikes;
AllSpikes{eid,j}.energy  = sum(DATA.AllData.Spikes.dVdt'.^2);
AllSpikes{eid,j}.svar = var(DATA.AllData.Spikes.values');
AllSpikes{eid,j}.times = DATA.AllData.Spikes.times;
end
T = DATA.Expts{eid}.Trials;
for j = 1:length(T)
for p = 1:length(DATA.probelist)
t = AllSpikes{eid,p}.times - T(j).Start(1);
tmax = T(j).End(end)- T(j).Start(1)+2000;
tid = find(t > -1000 & t < tmax);
T(j).AllSpikes(p).t = t(tid);
T(j).AllSpikes(p).energy = AllSpikes{eid,p}.energy(tid);
T(j).AllSpikes(p).var = AllSpikes{eid,p}.svar(tid);
end
end
DATA.Expts{eid}.Trials = T;     
end
setappdata(DATA.toplevel','AllSpikes',AllSpikes);
cmb.NotBusy(DATA);

