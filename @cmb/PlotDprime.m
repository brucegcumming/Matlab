function PlotDprime(DATA)

hold off;
Cx = DATA.Spikes.cx;
Cy = DATA.Spikes.cy;
Spks = DATA.AllData.Spikes;
cspks = DATA.Expts{DATA.currentexpt(1)}.gui.spks;
C = DATA.cluster{DATA.currentcluster,DATA.probe};
x = (Cx(cspks) - C.x(1))./C.x(3);
y = (Cy(cspks) - C.y(1))./C.y(3);
xr = x .* cos(C.angle) + y .* sin(C.angle);
yr = y .* cos(C.angle) - x .* sin(C.angle);
id = find(Spks.codes(cspks,2) == DATA.currentcluster);
nid = find(Spks.codes(cspks,2) ~= DATA.currentcluster);
xc = mean(xr(id));
yc = mean(yr(id));
sy = ((yr-mean(yr(id)))./std(yr));
sx = ((xr-mean(xr(id)))./std(xr));
o =C.dprimepar(1);
d = sx .* cos(o) + sy.* sin(o);
sd = std(d);
[y,x] = smhist(d,'sd',sd/10);
h = plot(sx,sy,'.','markersize',DATA.ptsize);
axis('equal');
refline(tan(o));
hold on;
yscale = diff(get(gca,'ylim'))./max(y);
rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
plot(rx,ry);
[y,x] = smhist(d(id),'sd',sd/10);
rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
plot(rx,ry,'r');
[y,x] = smhist(d(nid),'sd',sd/10);
rx = (x-xc).*cos(-o) + (y-yc).*sin(-o)*yscale;
ry = (y-yc).*cos(-o)*yscale - (x-xc).*sin(-o);
plot(rx,ry,'g');


