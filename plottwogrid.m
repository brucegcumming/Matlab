function fig = plottwogrid(file, interp, gaussfit , greyscale, fitfft, theplot, addfit,vgc)

ocolor = [0.5 0.5 0.5];
final = 1;
hold off

if nargin < 4
  theplot = 3;
  greyscale = 0;
  fitfft = 1
end

if nargin < 3
  gaussfit = 0;
elseif nargin < 2
  gaussfit = 0;
  doshade = 0;
  interp = 0;
end

if interp
   doshade = 1;
else
  doshade = 0;
end

opd = dlmread(file, ' ');
[rows, cols] =  size(opd);
orows = opd(1,3);
ntab = opd(1,4);
ori = opd(1,5);
orirad = ori * pi/180;
X = opd(2:orows+1,1:cols);
Y = opd(2+orows:2*orows+1 ,1:cols);
Z = opd(2+orows*2:3 * orows+1,1:cols);
DX = opd(2+orows*3:4*orows+1,1:cols);
DY = opd(2+orows*4:5*orows+1 ,1:cols);
if(rows > orows*8)
 ERR = opd(2+orows*7:8*orows+1 ,1:cols)
end

greys = [1 1 1; 0.84 0.84 0.84; 0.67 0.67 0.67; 0.5 0.5 0.5; 0.33 0.33 0.33; 0.17 0.17 0.17 ; 0 0 0; 0 0 0];

%Make copies of data with duplicates of the final row and column
%Then pcolor will plot _ALL_ of the data
xvals = max(X);
yvals = Y(:,1)
maxx = max(X,[],2);
minx = min(X,[],2);
xinc = (maxx(1)-minx(1))/(orows-1);
miny = min(Y,[],1);
maxy = max(Y,[],1);
yinc = (maxy(1)-miny(1))/(cols-1);
Xall = [X;X(orows,:)];
Xall = [Xall (Xall(:,cols) + xinc)];
Yall = [Y;(Y(orows,:)+yinc)];
Yall = [Yall Yall(:,cols)];
Zall = [Z;Z(orows,:)];
Zall = [Zall Zall(:,cols)];
zrange = [ max(max(Z)), min(min(Z))]
xrange = [ min(min(X)), max(max(X))];
yrange = [ min(min(Y)), max(max(Y))];
pfile = strrep(file,'.oxp','.gpm');
fp = dlmread(pfile, ' ')
% X,Y, SX,SY, F Phase Ori Amp base
% 1 2   3  4  5   6    7   8   9

if greyscale > 0
  colordata = [0:0.01:1 ; 0:0.01:1; 0:0.01:1]';
  set (gcf,'colormap',colordata);
else
  colormap default;
end

clf;
if fitfft
  ffile = strrep(file,'.oxp','.fitft');
  fopd = dlmread(ffile, ' ');
  [frows, fcols] =  size(fopd);
     forows = fopd(1,1);
  radius_pval = fopd(2,4);
  spori = fopd(2,2);
% If the fit was required a gaussian peak away from DC
% Use the angle to the peak for Orienation. Otherwise
% Use the axis of orientation of the Gaussian Blob.
  if radius_pval >= 0.05 | gaussfit
    fitori = fopd(1,3);
    FITZ = fopd(3+forows*4:5*forows+2,1:fcols);
  else
    fitori = fopd(2,1);
    FITZ = fopd(3+forows*3:4*forows+2,1:fcols);
  end
  fftdori = fitori + ori*pi/180;
  if ntab > 5
    nrpt = opd(6*orows+3,6);
  end
elseif ntab > 5
%Gabor fit in space.
 FITZ = opd(2+orows*5:6*orows+1 ,1:cols);
if(rows > orows*7)
 VMR = opd(2+orows*6:7*orows+1 ,1:cols);
 showvmr = 1;
 ntab = 7;
else
 ntab = 6;
end
 sx = opd(ntab*orows+2,3);
 sy = opd(ntab*orows+2,4);
 freq = opd(ntab*orows+2,5);
 fitori = opd(ntab*orows+2,7);
 residual = opd(ntab*orows+3,5);
 nrpt = opd(ntab*orows+3,6);
 aspect = opd(ntab*orows+3,2);
 fitori = fitori * 180/pi;
  aori = fitori;
 if abs(sx) < abs(sy) | abs(.2/freq) < abs(sy)
  aori = aori + 90;
 end
 dori = aori - 90 + ori;
 dorirad = pi * dori/180;
end %fitfft or ntab > 5


%Build interpolated matrices even if not interp - these are used
%for the contour plots
range = [ max(max(Z)), min(min(Z))];
 maxx = max(X,[],2);
 minx = min(X,[],2);
 xinc = (maxx(1)-minx(1))/(orows-1);
 xi = linspace(minx(1), maxx(1), 40);
 yc = (maxy(1)+miny(1))/2; 
 xc = (maxx(1)+minx(1))/2;
 yi = linspace(miny(1), maxy(1),40);
 [xxi, yyi] = meshgrid(xi,yi);
%Gauassian smoothing interpolation instead
if(xinc < yinc)
 smoothing = yinc/2;
else
 smoothing = xinc/2;
end
 zi = Interpf(X,Y,Z,xxi,yyi,1,smoothing);
range = [ max(max(zi)), min(min(zi))];

% X,Y, SX,SY, F Phase Ori Amp base
% 1 2   3  4  5   6    7   8   9

xfit = (xxi - fp(1,1)) * cos(fp(1,7)) + (yyi - fp(1,2)) * sin(fp(1,7));
yfit = (yyi - fp(1,2)) * cos(fp(1,7)) - (xxi - fp(1,1)) * sin(fp(1,7));
xg = exp(-(xfit .* xfit)/(2 * fp(1,3) * fp(1,3)));
yg = exp(-(yfit .* yfit)/(2 * fp(1,4) * fp(1,4)));
sfit = cos((xfit * 2 * pi * fp(1,5)) + fp(1,6));
zfit = fp(1,9) + fp(1,8) * (xg .* yg .* sfit);


if interp
 %matlab interolation
% zi = interp2(X,Y,Z,xxi,yyi,'cubic');
if (theplot == 1)
% subplot(2,2,2);
 fig = pcolor(xxi,yyi,zi);
%use this for Greyscale pictures
 if(doshade)
  shading('interp');
 end
%  hold on
%  contour(xxi,yyi,zi,5);
%  hold off
 axis('image');
caxis([min(range),max(range)]);
  set(gca,'fontsize',18);
title(sprintf('%s Or %.1f',file,ori));
xlabel('Orthogonal Disparity');
ylabel('Paralell Disparity');
axis('image');
end %if theplot ==1 

if ntab > 5
 fzi = Interpf(X,Y,FITZ,xxi,yyi,1,smoothing);
% fzi = interp2(X,Y,FITZ,xxi,yyi,'cubic');
if(theplot == 2)
% subplot(2,2,1)
 title('Fit');
 fig = pcolor(xxi,yyi,zfit);
 if doshade
  shading('interp');
 end % if doshade
 axis('image');
end %if(theplot == 2
end %if interp
else %not if ntab 
if (theplot == 1)
%  subplot(2,2,2) 
 fig = pcolor(Xall,Yall,Zall);
%   imagesc(xrange,yrange,Z);
 caxis([min(zrange),max(zrange)]);
  hold on


zirange = [ max(max(zi)), min(min(zi))];
nlines = 7;
for line = 1:nlines
  zval = max(zirange) - (line-1) * (max(zirange) - min(zirange))/nlines
  if(addfit)
    [C, h] = contour(xxi+xinc/2,yyi+yinc/2,zfit,[zval zval]);
  else
    [C, h] = contour(xxi+xinc/2,yyi+yinc/2,zi, [zval zval]);
  end
   for i = 1:length(h)
%    set(h(i),'edgeColor', greys(line,:),'LineWidth',2);
    set(h(i),'edgeColor', 'w','LineWidth',2);
   end
end
   caxis([min(zrange),max(zrange)]);
   hold off
   if(doshade)
      shading('interp');
   end
axis('image');
%set(fig,'EdgeColor', 'none');
set(gca,'fontsize',18);
title(sprintf('%s Or %.1f',file,ori));
xlabel('Orthogonal Disparity');
ylabel('Paralell Disparity');
end %if(theplot ==1)


myca = [1 0 1; 0 1 0; 0.6 0.5 0; 1 0 0;  0 0 1; 0 0 0.5 ; 0.5 0 0.5; 0 0.5 0];

if(theplot == 5)
  h = plot3(X',Y',Z','o-');
  cca = get(gca, 'ColorOrder');
  set(gca,'xgrid','on','ygrid','on','zgrid','on', 'Xlim', xrange, ...
	  'Ylim',yrange,'Xtick',xvals,'View',[35 15],'linewidth', ...
	  2,'fontsize',18,'view',[35 45]);
%34,15 is more like the default  
  hold on;
for j = 1:orows
  set(h(j),'MarkerFaceColor', myca(j,:), 'Color', myca(j,:),'LineWidth',3);

%show error bars
  for k= 1:cols
    E = [Z(j,k) + ERR(j,k); Z(j,k)-ERR(j,k)];
    YI = [Y(j,k); Y(j,k)];
    XI = [X(j,k); X(j,k)];
    plot3(XI,YI,E, 'Color', myca(j,:),'LineWidth',3);
  end

  % Add a line underneath in the same color
    XI = [X(j,1); X(j,1); X(j,cols); X(j,cols)];
    YI = [Y(j,1); Y(j,1); Y(j,cols); Y(j,cols)];
    E = [Z(j,1); 0; 0; Z(j,cols)];
    plot3(XI,YI,E, 'Color', myca(j,:),'LineWidth',2);
end
  hold off;


  if (final == 0)
    title(sprintf('%s Or %.1f',file,ori));
  end
xlabel({'Orthogonal Disparity','(degrees)'},'Rotation',340);
ylabel({'Paralell Disparity','  (degrees)'},'Rotation',35);
zlabel('Firing rate');
end %if theplot == 5

if (theplot == 2)
  if ntab > 5
    axis('image');
%    subplot(2,2,1);
    fig = pcolor(X,Y,FITZ); 
    set(fig,'EdgeColor', 'none');
 end
end %(if theplot == 2)
end

if(theplot ==2)
axis('image');
if(fitfft)
  title(sprintf('Fit Ori %.3f = %.3f',fitori*180/pi,spori*180/pi));
else
   title(sprintf('Fit: O %.2f:%.2f(%.2f) S %.3f %.3f F %.2f R %.1f',dori,fitori,aori,sx,sy,freq,residual));
end
end %(if theplot == 2)

if(theplot == 3)
if ntab > 5
%subplot(2,2,4);
else
%subplot(2,2,2);
end
ylabel('Vertical Disparity');
xlabel('Horizontal Disparity');
minx = min(min(DX));
maxx = max(max(DX));
miny = min(min(DY));
maxy = max(max(DY));
yc = (maxy(1)+miny(1))/2; 
xc = (maxx(1)+minx(1))/2; 

%Calculate intepolated grid and fit on this grid regardless of
%whether intepolated data is shown. - Need it for the fit as a contour

xi = linspace(minx(1), maxx(1), 40);
yi = linspace(miny(1), maxy(1),40);
[xxi, yyi] = meshgrid(xi,yi);

%Now for each xxi,yyi, convert back into dP,dO before calculating the
%fit

ooi = xxi * sin(orirad) - yyi * cos(orirad);
ppi = yyi * sin(orirad) + xxi * cos(orirad);
xfit = (ooi - fp(1,1)) * cos(fp(1,7)) + (ppi - fp(1,2)) * sin(fp(1,7));
yfit = (ppi - fp(1,2)) * cos(fp(1,7)) - (ooi - fp(1,1)) * sin(fp(1,7));
xg = exp(-(xfit .* xfit)/(2 * fp(1,3) * fp(1,3)));
yg = exp(-(yfit .* yfit)/(2 * fp(1,4) * fp(1,4)));
sfit = cos((xfit * 2 * pi * fp(1,5)) + fp(1,6));
zfit = fp(1,9) + fp(1,8) * (xg .* yg .* sfit);


subplot(2,1,1);
if interp

  zi = Interpf(DX,DY,Z,xxi,yyi,1,smoothing);
  %zi = griddata(DX,DY,Z,xxi,yyi,'cubic');
  fig = pcolor(xxi,yyi,zi);
  axis('image');
  if(fitfft)
   title(sprintf('Disp Ori %.3f (%.3f)',fftdori,fftdori*180/pi));
  else
   title(sprintf('Disparity Orientation %.0f',dori));
   tca = text(xc,yc,'monocular orientation tuning');
  end
  if(doshade)
   shading('interp');
  end
else
  fig = pcolor(DX,DY,Z);
  hold on
  if(addfit)
    [C , h] = contour(xxi+xinc/2,yyi+yinc/2,zfit,7);
  else
    [C, h] = contour(xxi+xinc/2,yyi+yinc/2,zi,7);
  end

  for i = 1:length(h)
     set(h(i),'edgeColor', 'w','LineWidth',2);
  end

  hold off
  title(sprintf('Disparity Orientation %.0f',dori));
  set(gca,'fontsize',18);
  tca = text(xc,yc,{'orientation tuning','(monocular)'});
  caxis([min(zrange),max(zrange)]);
  set(tca,'color',ocolor);
  if(doshade)
    shading('interp');
    caxis([min(range),max(range)]);
  end
end %if/else interp

set(fig,'EdgeColor', 'none');
axis('image');
xlabel('Horizontal Disparity');
ylabel('Vertical Disparity');
%rotate(fig,[0 90],45);

hold on

if(fitfft)
OMARKX = [minx, maxx, miny/tan(fftdori),maxy/tan(fftdori)] +xc;
OMARKY = [minx * tan(fftdori), maxx * tan(fftdori), miny, maxy]+yc;
else
OMARKX = [minx, maxx, miny/tan(dorirad),maxy/tan(dorirad)]+xc;
OMARKY = [minx * tan(dorirad), maxx * tan(dorirad), miny, maxy]+yc;
end
%plot(OMARKX, OMARKY, 'w');

otfile = strrep(file,'rds.OxPD.oxp','grating.prefOT.ma');
if exist(otfile)
OR =  dlmread(otfile,' ');
A = OR(:,2);
R = OR(:,1) * pi/180.0;
[XX, YY] = pol2cart(R,A);
if(abs(cos(orirad)) > 0.3 & abs(sin(orirad)) > 0.3) %Oblique
scale = (maxx - minx)/( 3 * max(A));
else
scale = (maxx - minx)/( 2 * max(A));
end
XX = (XX * scale) + xc;
YY = (YY * scale) + yc;
plot(XX, YY,'color', ocolor, 'linewidth', 2);
%plot(XX, YY);
%axis('image');
end
hold off
subplot(2,1,2);
hist(vgc);
axis([minx maxx 0 100]);
axis('auto y');
axis('square');
end %if(theplot == 3)


FZ = fft2(Z);
FZ(1,1) = 0;
SFZ = fftshift(FZ);
AZ=abs(SFZ);
if(theplot == 4)
%subplot(2,2,3);
%pcolor(X,Y,AZ);

outfile = strrep(file,'.oxp','.fft');
fout = fopen(outfile,'w');
for n = 1:orows
for m = 1:cols
fprintf(fout, '%.3f %.3f %.3f\n',X(n,m),Y(n,m),AZ(n,m));
if AZ(n,m) == 0
centre(1) = m;
centre(2) = n;
end
end
end
fclose(fout);


%
% one day must do this with matrices. 
% r = 1:rows;
% x = centre(1) - r
% ?Generate two matrices of fractions

[rows, cols] =  size(AZ);
showpolar = 1;
if(showpolar)

A = [];
R = [];
for angle= -pi/4:0.02:pi/4
  z = 0;
n = 0;
for i = 1:cols;
  x = i - centre(1);
  y = x * tan(angle);
  j = centre(2) + floor(y);
  k = centre(2) + ceil(y);
if(k <= rows & j > 0)
  frac = y + centre(2) -j;
  z = z + (AZ(j,i) * (1-frac) + AZ(k,i) * frac);
  n = n+1;
end
end

%fprintf('%f %f\n',angle,z/n);
A = [A  angle];
R = [R  z/n];
end

for angle= pi/4:0.02:(3*pi)/4
  z = 0;
n = 0;
for i = 1:rows;
  y = i - centre(2);
  x = y * cot(angle);
  j = centre(1) + floor(x);
  k = centre(1) + ceil(x);
if(k <= cols & j > 0)
  frac = x + centre(1) -j;
  z = z + (AZ(i,j) * (1-frac) + AZ(i,k) * frac);
  n = n+1;
end
end

A = [A  angle];
R = [R  z/n];
end

if(fitfft)
[peak, peaki] = max(R);
fprintf('peak at %f\n',A(peaki) * 180.0/pi);
C = [A(peaki) A(peaki)+pi];
D = [peak peak];
B = A + pi;
A = [A B];
R = [R R];
hold off;
pcolor(X,Y,AZ);
title(sprintf('2D-FFT: Peak at %.1f',A(peaki) * 180/pi));
h1 = gca;
hold on;
h2 = axes('Position',get(h1,'Position'));
[XX, YY] = pol2cart(A,R);
[CX, CY] = pol2cart(C,D);
plot(XX, YY, 'w');
axis('image');
axis('off');
set(h2,'Color','none');
h2 = axes('Position',get(h1,'Position'));
plot(CX, CY, 'w');
axis('image');
axis('off');
set(h2,'Color','none');
hold off;
else
if(showvmr)
  range = [ max(max(VMR)), min(min(VMR))];
  pcolor(X,Y,VMR);
  title(sprintf('%s %.1f - %.1f',file,min(range),max(range)));
end
end
end
if(interp)
fprintf( '%.3f %.3f %.3f\n',xinc,yinc,smoothing);
end
end %if(theplot == 4)

hold off
