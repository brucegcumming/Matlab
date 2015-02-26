function ellipse(ym, cc)
%Draws an ellipse interactively
global mousept;

ff = GetFigure('figtest');
if nargin == 0
    ym = 1;
end
plot([0 1],[0 ym]);
axis('manual');
set(gca,'ylim',[0 ym]);
hold on;
if nargin < 2
mousept.start = [];
mousept.lasth = [];
mousept.angle = 0;
mousept.r = [1 1];
else
    mousept.c = [0.5 ym/2];
    mousept.start = mousept.c;
    mousept.r = cc(1,1:2);
    mousept.angle = cc(1,3);
    mousept.mode = 5;
    myellipse(mousept);
end
mousept.down = 0;
mousept.ratio = ym;
for j = 1:size(cc,1)
    e.r(1) = cc(j,1);
    e.r(2) = cc(j,2);
    e.angle = -cc(j,3);
    e.c = mousept.c;
    DrawEllipse(e);
   [x,y] = exy(e.angle,e.r(1),e.r(2),ym);
    fprintf('%.3f %.3f %.3f l = %.3f\n',e.angle,x,y,sqrt(x^2+y^2));
 end
if size(cc,1) > 1
    return;
end
%for angle = [0 pi/16 pi/8 3 * pi/8 11 * pi/24 pi/2]
for angle = [pi/20 pi/4 pi/2]
%    angle = atan(tan(angle * 10));
    [x,y] = exy(angle,mousept.r(1),mousept.r(2));
    plot(x+mousept.c(1),y+mousept.c(2),'o');
    fprintf('%.3f %.3f %.3f l = %.3f\n',angle,x,y,sqrt(x^2+y^2));
    hold on;
    [a,b] = exy(angle,mousept.r(1),mousept.r(2));
    e.r(1) = a;
    e.r(2) = b;
    e.angle = -angle;
    e.c = mousept.c + [0.1 0.1];
    DrawEllipse(e);
end


function [x,y] = exy(t,a,b,kax)

tt = tan(t) * kax;
tt = tan(t);
%tt = tan(t) * sqrt(kax);
x = sqrt((kax^2 * a^2)/((tt^2 + kax^2)));
y = kax * sqrt(a^2-x^2);
x = sqrt(x^2+y^2);
tt = tan(t+pi/2);
%tt = tan(t+pi/2)/(sqrt(kax));
ax = sqrt(b^2/(tt^2 + kax^2));
ay = sqrt(b^2 - kax^2*ax^2);
y = sqrt(ax^2+ay^2);;
plot(ax,ay,'ro');

function [x,y] = oldexy(t,a,b)

kax = 1;
x = sqrt(b^2/((tan(t)^2) + b^2/a^2));
y = sqrt(b^2 - x^2*b^2/a^2);
x = sqrt(x^2+(kax * y)^2);
ax = sqrt(b^2/((tan(t+pi/2)^2) + b^2/a^2));
ay = sqrt(b^2 - ax^2*b^2/a^2);
y = sqrt(ax^2+(ay/kax)^2);

function ButtonPressed(src, data)
global mousept;
mousept.lasth = 0;
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'extend' 'alternate' 'open'})
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');
end
mousept.down = 1;

function ButtonReleased(src, data)
global mousept;
mousept.down = 0
    pt = get(gca,'CurrentPoint');
   mousept= myellipse(mousept,pt);


function ButtonDragged(src, data)

global mousept;

if mousept.down;
    pt = get(gca,'CurrentPoint');
    if mousept.lasth
        delete(mousept.lasth);
        delete(mousept.xl);
        delete(mousept.yl);
    end
    mousept= myellipse(mousept,pt);
end


function h = ellipse(mousept, finish)

start = mousept.start;
circle = rsmak('circle');
if mousept.mode == 1
    a = finish(1,1)-start(1,1);
    b = finish(2,2)-start(2,2);
    mousept.r = [a b]
elseif mousept.mode == 2
    mousept.angle = atan(diff(finish)/diff(start));
    [a,b] = mousept.r;
elseif mousept.mode == 3
    mousept.angle = atan(diff(finish)/diff(start));
    a = mousept.r(1);
    b = mousept.r(2);
else
    a = finish(1,1)-start(1,1);
    b = finish(2,2)-start(2,2);
    mousept.r = [a b]
end
ellipse = fncmb(circle,[a 0; 0 b]);
sn = sin(mousept.angle);
cn = cos(mousept.angle);
ellipse = fncmb(ellipse,[cn -1/cn; cn 1/cn]);
ellipse = fncmb(ellipse,[start(1,1); start(2,2)]);
%ellipse = fncmb(ellipse,[start(1); start(2)]);
pts = fnplt(ellipse);
h = plot(pts(1,:),pts(2,:));
%drawnow;
%h = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)]);

function mousept = myellipse(mousept, finish)

start = mousept.start;
if mousept.mode == 1
    a = (finish(1,1)-start(1,1))/2;
    b = (finish(2,2)-start(2,2))/2;
    mousept.r = [a b];
    mousept.c = [finish(1,1)+start(1,1) finish(2,2)+start(2,2)]/2;
elseif mousept.mode == 2
    a = (finish(1,1)-start(1,1))/2;
    b = (finish(2,2)-start(2,2))/2;
    mousept.r = [a b];
elseif mousept.mode == 3
    
    dy = (finish(2,2)-mousept.c(2));
    dx = (finish(1,1)-mousept.c(1));
    t = 0.1 * -dy/dx;
    mousept.angle = atan(0.1 * -dy/dx);
    [a,b] = exy(mousept.angle,mousept.r(1),mousept.r(2));
    fprintf('%.3f %.3f %.3f\n',mousept.angle,a,b);
%    a = (mousept.r(1) * cos(mousept.angle)) + mousept.ratio * mousept.r(1) * sin(mousept.angle)
%    b = (mousept.r(2) * cos(mousept.angle)) - (mousept.r(2) * sin(mousept.angle)/mousept.ratio)
    %b = mousept.r(2);
elseif mousept.mode == 4  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
else
    a = mousept.r(1);
    b = mousept.r(2);
end
moustpt = DrawEllipse(mousept);



function e = DrawEllipse(e, varargin)

j = 1;
lineonly = 0;
while j <= nargin-1
    if strncmpi(varargin{j},'lineonly',5)
        lineonly = 1;
    elseif strncmpi(varargin{j},'linex',5)
        lineonly = 2;
    end
    j= j+1;
end
sn = sin(e.angle);
cn = cos(e.angle);
%sn = sin(pi/4);
%cn = cos(pi/4);
a = e.r(1);
b = e.r(2);
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn) + e.c(1);
yr = (y .* cn - x .*sn) + e.c(2);
if lineonly ==  0
    e.lasth = plot(xr,yr);
end
x = [-a a];
y =  sqrt(b.^2 - (x.*b/a).^2);
xr = (x .* cn + y .*sn) + e.c(1);
yr = (y .* cn - x .*sn) + e.c(2);
e.xl = plot(xr,yr);
y = [b -b];
x =  sqrt(a.^2 - (y.*a/b).^2);
xr = (x .* cn + y .*sn) + e.c(1);
yr = (y .* cn - x .*sn) + e.c(2);
if ismember(lineonly,[0 1 3])
    e.yl = plot(xr,yr,'r');
end
%drawnow;
%h = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)]);

