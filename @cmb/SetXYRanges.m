function DATA = SetXYRanges(DATA, cx, cy)


if nargin == 1
cx = DATA.Spikes.cx(DATA.spklist);
cy = DATA.Spikes.cy(DATA.spklist);
end

if DATA.plot.autoscale == 0
return;
elseif DATA.plot.autoscale == 1 || length(cx) == 0
DATA.plot.clusterXrange  = get(gca,'Xlim');
DATA.plot.clusterYrange  = get(gca,'Ylim');
elseif ismember(DATA.plot.autoscale,[2 3 4])
if min(cx) > 0
minx = 0;
else
minx = prctile(cx,1 * 1.1);
end
if min(cy) > 0
miny = 0;
else
miny = prctile(cy,1 * 1.1);
end
if DATA.plot.autoscale == 2
DATA.plot.clusterYrange = [miny prctile(cy,99.9).*1.2];
DATA.plot.clusterXrange = [minx prctile(cx,99.9).*1.2];
elseif DATA.plot.autoscale == 3
DATA.plot.clusterYrange = [miny prctile(cy,99.5).*1.2];
DATA.plot.clusterXrange = [minx prctile(cx,99.5).*1.2];
else
DATA.plot.clusterYrange = [miny prctile(cy,95).*2];
DATA.plot.clusterXrange = [minx prctile(cx,95).*2];
end
if min(cx) < 0
DATA.plot.clusterXrange(1) = prctile(cx,1) * 1.1;
end
if min(cy) < 0
DATA.plot.clusterYrange(1) = prctile(cy,1) * 1.1;
end
elseif DATA.plot.autoscale == 4
DATA.plot.clusterYrange = [min(cy) max(cy)];
DATA.plot.clusterXrange = [min(cx) max(cx)];
end


