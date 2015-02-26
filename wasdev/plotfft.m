function fig = plotgrid(file, interp, doshade )

if nargin < 3
  doshade = 0;
  interp = 0;
elseif nargin < 2
  doshade = 0;
end

pori = -1;
opd = dlmread(file, ' ');
[rows, cols] =  size(opd);
orows = opd(1,1)
fitori = opd(1,3);
spori = opd(2,2);
X = opd(3:orows+2,1:cols);
Y = opd(3+orows:2*orows+2 ,1:cols)
Z = opd(3+orows*2:3 * orows+2,1:cols)
FZ = opd(3+orows*3:4*orows+2,1:cols);
range = [ max(max(Z)), min(min(Z)), max(max(FZ)), min(min(FZ))];

subplot(1,2,2);
fig = pcolor(X,Y,FZ);
axis('square');
caxis([min(range),max(range)]);
title(sprintf('Fit Ori %.3f at %.3f,%.3f',fitori*180/pi,opd(1,4),opd(1,5)));
tx = get(gca,'Title');
tp = get(tx,'Position')
%text(min(min(X)),max(max(Y))*1.25,sprintf('SD %.3f,%.3f',opd(1,6),opd(1,7)));
text(tp(1),tp(2)+0.2,sprintf('SD %.3f,%.3f',opd(1,6),opd(1,7)));

subplot(1,2,1);
fig = pcolor(X,Y,Z);
axis('image');
caxis([min(range),max(range)]);

[zmin, minrow] = min(Z)
[tmp, mincol] = min(zmin);
centre(2) = minrow(mincol);
centre(1) = mincol;
centre
[rows, cols] =  size(Z)

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
  z = z + (Z(j,i) * (1-frac) + Z(k,i) * frac);
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
  z = z + (Z(i,j) * (1-frac) + Z(i,k) * frac);
  n = n+1;
end
end
A = [A  angle];
R = [R  z/n];
end
[peak, peaki] = max(R);
%title(sprintf('Data'));
title(sprintf('Data: Peak at %.1f',A(peaki) * 180/pi));
fprintf('peak at %f\n',A(peaki) * 180.0/pi);
C = [A(peaki) A(peaki)+pi];
D = [peak peak];
B = A + pi;
A = [A B];
R = [R R];
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
axis('image');
plot(CX, CY, 'w');
axis('image');
axis('off');
set(h2,'Color','none');
hold off;

