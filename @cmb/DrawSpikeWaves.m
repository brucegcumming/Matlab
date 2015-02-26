function DrawSpikeWaves(DATA, ispk, nclusters, ctype)
if  (DATA.svfig == 0 || DATA.state.nospikes) 
return;
end

for j = 1:nclusters+1
vs{j} = [];
xs{j} = [];
end
for j = 1:length(DATA.svh)
vs{j} = [];
xs{j} = [];
end
needdv = 0;

if isfield(DATA,'AllSpikes')
splen = size(DATA.AllSpikes{DATA.probe}.values,2).*size(DATA.AllSpikes{DATA.probe}.values,3);
adc = DATA.AllSpikes{DATA.probe}.values(ispk,:);
codes = DATA.AllSpikes{DATA.probe}.codes(ispk,:);
if needdv
dvdt = DATA.AllSpikes{DATA.probe}.dVdt(ispk,:);
end
else
splen = size(DATA.AllData.Spikes.values,2).*size(DATA.AllData.Spikes.values,3);
adc = DATA.AllData.Spikes.values(ispk,:);
codes = DATA.AllData.Spikes.codes(ispk,:);
dvdt = DATA.AllData.Spikes.dVdt(ispk,:);
end
for spk = 1:length(ispk);
j = codes(spk, ctype)+1;
if DATA.plot.dvdt
vs{j} = [vs{j} dvdt(spk,:) NaN];
xs{j} = [xs{j} [1:splen-1] NaN];
elseif DATA.plot.voltxy
vs{j} = [vs{j} adc(spk,33:64) NaN];
xs{j} = [xs{j} adc(spk,65:96) NaN];
else
vs{j} = [vs{j} adc(spk,:) NaN];
xs{j} = [xs{j} [1:splen] NaN];
end
end
set(0,'CurrentFigure',DATA.svfig);
nc = min([nclusters+1 length(DATA.svh)]);
for j = 1:max([nc length(DATA.svh)])
if ~isempty(xs{j})
if ishandle(DATA.svh(j))
set(DATA.svh(j),'Xdata' , xs{j}, 'Ydata', vs{j});
else
DATA.svh(j) = line('Xdata' , xs{j}, 'Ydata', vs{j});
end
elseif ishandle(DATA.svh(j)) %no spikes, wipe clean
set(DATA.svh(j),'Xdata' , 0, 'Ydata', 0);
end
end
drawnow;

