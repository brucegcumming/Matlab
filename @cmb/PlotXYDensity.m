function  PlotXYDensity(energy,vw)

if length(vw) < 10
return;
end
if length(vw) > 100000
lprc = 0.01;
hprc = 99.99;
sx=3;
sy=3;
elseif length(vw) > 10000
lprc = 0.1;
hprc = 99.9;
sx=3;
sy=3;
elseif length(vw) > 1000
lprc = 1;
hprc = 99;
sx=5;
sy=5;
else
lprc = 5;
hprc = 95;
sx=8;
sy=8;
end    

erange = [prctile(energy,lprc) prctile(energy,hprc)];
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
mode = 2;
tic;
if mode ==1 % add real gaussian to grid for each
z = zeros(size(x));
sx = (diff(erange)/100)^2;
sy = (diff(vrange)/100)^2;
for j=1:length(energy)
z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
end
elseif mode ==2 %build fine 2-D histogram, then smooth
[gx,gy] = meshgrid(-10:10,-10:10);

G = exp(-(gx).^2/sx - (gy).^2/sy);
G = G./sum(G(:));
z = zeros(size(x));
vi = 1+floor(nbins * (vw-vrange(1))/diff(vrange));
ei = 1+floor(nbins * (energy-erange(1))/diff(erange));
idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
for j =idx
z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
end
z = conv2(z,G,'same');
end
toc
pcolor(x,y,z);
shading('interp')


