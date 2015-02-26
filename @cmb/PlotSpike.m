function PlotSpike(DATA, ispk, probe)

set(0,'CurrentFigure',DATA.svfig);
if nargin == 2
j = DATA.AllData.Spikes.codes(ispk,2)+1;
if DATA.plot.dvdt == 2
set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',DATA.AllData.Spikes.values(ispk,2:end));
elseif DATA.plot.dvdt
set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.dVdt(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.dVdt,2)]);
else
set(DATA.svh(j), 'Ydata', DATA.AllData.Spikes.values(ispk,:),'Xdata',[1:size(DATA.AllData.Spikes.values(ispk,:),2)]);
end
title(sprintf('%d: Cl %d at %.3f',ispk,j-1,DATA.AllData.Spikes.times(ispk)/10000)); 
else
j = DATA.AllSpikes{DATA.probe}.codes(ispk,2)+1;
set(DATA.svh(j), 'Ydata', DATA.AllSpikes{probe(1)}.values(ispk,:)-DATA.plot.SpikeMaxV/2,'Xdata',[1:size(DATA.AllSpikes{probe(1)}.values,2)]);
id = find(ispk == DATA.sids{1});
if length(id)
jspk = DATA.sids{2}(id);
set(DATA.svh(j+DATA.nclusters+1), 'Ydata', DATA.AllSpikes{probe(2)}.values(jspk,:)+DATA.plot.SpikeMaxV/2,'Xdata',[1:size(DATA.AllSpikes{probe(2)}.values,2)]);
end
title(sprintf('%d: Cl %d at %.3f',ispk,j-1,DATA.AllSpikes{DATA.probe}.times(ispk)/10000));
end
drawnow;

