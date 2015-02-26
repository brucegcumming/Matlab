function figtest(ym, cc)
global mousept;

ff = GetFigure('figtest');
if nargin == 0
    ym = 1;
end
plot([0 1],[0 ym]);
set(ff, 'WindowButtonDownFcn',@ButtonPressed);
set(ff, 'WindowButtonMotionFcn',@ButtonDragged);
set(ff, 'WindowButtonUpFcn',@ButtonReleased);
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
    mousept.r = cc(1:2);
    mousept.angle = cc(3);
    mousept.mode = 5;
    myellipse(mousept);
end
mousept.down = 0;
mousept.ratio = ym;
mousept.xrange = 1;
mousept.yrange = ym;

for angle = 0:pi/100:pi/2
    [x,y] = exy(angle,mousept.r(1),mousept.r(2));
    plot(x+mousept.c(1),y+mousept.c(2),'o');
    fprintf('%.3f %.3f %.3f\n',angle,x,y);
    hold on;
end


function [x,y] = exy(t,a,b)

x = sqrt(b^2/(tan(t)^2 + b^2/a^2));
y = sqrt(b^2 - x^2*b^2/a^2);

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
    
    dy = (finish(2,2)-mousept.c(2))./mousept.yrange;
    dx = (finish(1,1)-mousept.c(1))./mousept.xrange;
    t = -dy/dx;
    mousept.angle = atan(-dy/dx);
    a = mousept.r(1);
    b = mousept.r(2);
%   [a,b] = exy(mousept.angle,mousept.r(1),mousept.r(2));
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

a = a./mousept.xrange;
b = b./mousept.yrange;
sn = sin(mousept.angle);
cn = cos(mousept.angle);
%sn = sin(pi/4);
%cn = cos(pi/4);
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = mousept.xrange .* (x .* cn + y .*sn) + mousept.c(1);
yr = mousept.yrange .* (y .* cn - x .*sn) + mousept.c(2);
mousept.lasth = plot(xr,yr);
x = [-a a];
y =  sqrt(b.^2 - (x.*b/a).^2);
xr = mousept.xrange.*(x .* cn + y .*sn) + mousept.c(1);
yr = mousept.yrange.*(y .* cn - x .*sn) + mousept.c(2);
mousept.xl = plot(xr,yr);
y = [b -b];
x =  sqrt(a.^2 - (y.*a/b).^2);
xr = mousept.xrange .* (x .* cn + y .*sn) + mousept.c(1);
yr = mousept.xrange .* (y .* cn - x .*sn) + mousept.c(2);
mousept.yl = plot(xr,yr);
%drawnow;
%h = plot([start(1,1) finish(1,1)],[start(2,2) finish(2,2)]);

