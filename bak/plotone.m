function fig = plotgrid(file, interp, doshade )

if nargin < 3
  doshade = 0;
  interp = 0;
elseif nargin < 2
  doshade = 0;
end

opd = dlmread(file, ' ');
[rows, cols] =  size(opd);
orows = opd(1,3);
ntab = opd(1,4);
ori = opd(1,5);
X = opd(2:orows+1,1:cols);
Y = opd(2+orows:2*orows+1 ,1:cols);
Z = opd(2+orows*2:3 * orows+1,1:cols);
DX = opd(2+orows*3:4*orows+1,1:cols);
DY = opd(2+orows*4:5*orows+1 ,1:cols);

subplot(1,2,1);
if interp
 minx = min(X,[],2);
 maxx = max(X,[],2);
 xi = linspace(minx(1), maxx(1), 40);
 miny = min(Y,[],1);
 maxy = max(Y,[],1);
 yi = linspace(miny(1), maxy(1),40);
 [xxi, yyi] = meshgrid(xi,yi);
%Lances interpolation instead
% zi = Interpf(X,Y,Z,xxi,yyi,1,1);
%matlab interolation
 zi = interp2(X,Y,Z,xxi,yyi,'cubic');
 fig = pcolor(xxi,yyi,zi);
%use this for Greyscale pictures
%colordata = [0:0.01:1 ; 0:0.01:1; 0:0.01:1]';
%set (gcf,'colormap',colordata);
 if(doshade)
  shading('interp');
 end
 axis('image');
else
  fig = pcolor(X,Y,Z);
end

axis('image');
xlabel('Orthogonal Disparity');
ylabel('Paralell Disparity');

subplot(1,2,2);
if interp
minx = min(min(DX));
maxx = max(max(DX));
xi = linspace(minx(1), maxx(1), 40);
miny = min(min(DY));
maxy = max(max(DY));
yi = linspace(miny(1), maxy(1),40);
[xxi, yyi] = meshgrid(xi,yi);
zi = griddata(DX,DY,Z,xxi,yyi,'cubic');
fig = pcolor(xxi,yyi,zi);
axis('image');
if(doshade)
shading('interp');
end
else
fig = pcolor(xxi,yyi,zi);
end
axis('image');
ylabel('Vertical Disparity');
xlabel('Horizontal Disparity');

