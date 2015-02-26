function RescaleClusterPlot(a,b)
DATA = GetDataFromFig(a);
it = findobj(DATA.xyfig,'Tag','ClusterXmax');
if it
DATA.plot.clusterXrange(2) = str2num(get(it,'string'));
end

it = findobj(DATA.xyfig,'Tag','ClusterYmax');
if it
DATA.plot.clusterYrange(2) = str2num(get(it,'string'));
end
it = findobj(DATA.xyfig,'Tag','ClusterZmax');
if it
DATA.plot.clusterZrange(2) = str2num(get(it,'string'));
end
ax = findobj(DATA.xyfig,'type','ax');

%
%if the user has typed in a value, surely they want to manually scales
DATA.plot.autoscale = 0;
set(findobj(DATA.xyfig,'Tag','AutoScale'),'value',DATA.plot.autoscale);
if DATA.plot.autoscale == 0
set(ax,'Xlim',DATA.plot.clusterXrange);
set(ax,'Ylim',DATA.plot.clusterYrange);
end



if DATA.densityplot
caxis([0 DATA.plot.clusterZrange(2)]);
end
set(DATA.toplevel,'UserData',DATA);

