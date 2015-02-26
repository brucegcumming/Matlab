function PlotDDF(DATA)
cspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
plottype = 0;
nclusters = size(DATA.cluster,1);
p = DATA.probe;
for cl = nclusters:-1:1
id = find(DATA.AllData.Spikes.codes(cspks,2) == cl);
nid = find(DATA.AllData.Spikes.codes(cspks,2) ~= cl);
sx = std(DATA.Spikes.cx(cspks(id)));
sy = std(DATA.Spikes.cy(cspks(id)));
mx = mean(DATA.Spikes.cx(cspks(id)));
my = mean(DATA.Spikes.cy(cspks(id)));
if p <= size(DATA.cluster,2)
C = DATA.cluster{cl,p};
x = (DATA.Spikes.cx(cspks) - mx)./sx;
y = (DATA.Spikes.cy(cspks) - my)./sy;
xr = x .* cos(C.angle) + y .* sin(C.angle);
yr = y .* cos(C.angle) - x .* sin(C.angle);
d = sqrt(((yr.^2 + xr.^2)));
dprime = (mean(d(id)) - mean(d(nid)))./sqrt(mean([var(d(id)) var(d(nid))]));
x = (DATA.Spikes.cx(cspks) - C.x(1))./C.x(2);
y = (DATA.Spikes.cy(cspks) - C.y(1))./C.y(2);
ddf = sqrt((y).^2 + (x).^2);
GetFigure('DDF');
if plottype == 1
hold off;
plot(xr(id),yr(id),'r.');
hold on;
plot(xr(nid),yr(nid),'.');
axis('image');
elseif plottype == 2
hist(ddf,500);
else
hold off;
[y,x] = smhist(ddf,'sd',0.1,'xval',[0:0.1:10]);
plot(x,y);
hold on;
[y,x] = smhist(ddf(id),'sd',0.1,'xval',[0:0.1:10]);
plot(x,y,'r');
ddfprime = (mean(ddf(id)) - mean(ddf(nid)))./sqrt(mean([var(ddf(id))  var(ddf(nid))]));

x = (DATA.Spikes.cx(cspks) - C.x(1))./C.x(3);
y = (DATA.Spikes.cy(cspks) - C.y(1))./C.y(3);
xr = x .* cos(C.angle) + y .* sin(C.angle);
yr = y .* cos(C.angle) - x .* sin(C.angle);
ddf = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
drprime = (mean(ddf(id)) - mean(ddf(nid)))./sqrt(mean([var(ddf(id))  var(ddf(nid))]));

[y,x] = smhist(d,'sd',0.1,'xval',[0:0.1:10]);
plot(x,y,':');
hold on;
[y,x] = smhist(d(id),'sd',0.1,'xval',[0:0.1:10]);
plot(x,y,'r:');
end
title(sprintf('D = %.4f (%.4f,%.4f)',-dprime,-drprime,-ddfprime))
id = cspks(find(d < 1));
ddf = (yr./C.y(2)*C.y(3)).^2 + (xr./C.x(2)*C.x(3)).^2;
end
end

