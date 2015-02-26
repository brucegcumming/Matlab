function mousept = myellipse(mousept, finish)


%
%
% mousept modes
%  1 draw ellipse
%  2 re-size both dimensions
%  3 rotate
%  4 move ellipse
%  5 change left edge (move center and h radius)
if nargin > 1
    start = mousept.start;
else
    mousept.mode = 4;
end
if mousept.mode == 1
    a = (finish(1,1)-start(1,1))/2;
    b = (finish(2,2)-start(2,2))/2;
    mousept.r = abs([a b]);
    mousept.c = [finish(1,1)+start(1,1) finish(2,2)+start(2,2)]/2;
elseif mousept.mode == 2
    a = (finish(1,1)-mousept.c(1));
    b = (finish(2,2)-mousept.c(2));
    mousept.r = abs([a b]);
elseif mousept.mode == 3
    
    dy = (finish(2,2)-mousept.c(2))./mousept.yrange;
    dx = (finish(1,1)-mousept.c(1))./mousept.xrange;
    t = -dy/dx;
    mousept.angle = atan(-dy/dx);
    a = mousept.r(1);
    b = mousept.r(2);
    %   [a,b] = exy(mousept.angle,mousept.r(1),mousept.r(2));
    %    fprintf('%.3f %.3f %.3f\n',mousept.angle,a,b);
    %    a = (mousept.r(1) * cos(mousept.angle)) + mousept.ratio * mousept.r(1) * sin(mousept.angle)
    %    b = (mousept.r(2) * cos(mousept.angle)) - (mousept.r(2) * sin(mousept.angle)/mousept.ratio)
    %b = mousept.r(2);
elseif mousept.mode == 4  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
elseif mousept.mode == 11  %% move ellipse
    a = mousept.r(1);
    b = mousept.r(2);
    mousept.c(1) = finish(1,1)-mousept.offset(1);
    mousept.c(2) = finish(2,2)-mousept.offset(2);
elseif mousept.mode == 6 %% move R boundary
    a = (finish(1,1)-mousept.c(1));
    mousept.c(1) = mousept.c(1) + (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 7 %% move L boundary
    a = abs(mousept.c(1) - finish(1,1));
    mousept.c(1) = mousept.c(1) - (a-mousept.r(1))/2;
    mousept.r(1) = abs(finish(1,1) - mousept.c(1));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 8 %% move Top boundary
    b = (finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) + (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);
elseif mousept.mode == 9 %% move Bottom boundary
    b = abs(finish(2,2)-mousept.c(2));
    mousept.c(2) = mousept.c(2) - (b-mousept.r(2))/2;
    mousept.r(2) = abs(finish(2,2) - mousept.c(2));
    b = mousept.r(2);
    a = mousept.r(1);
    
else %just draw what we have
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
mousept.lasth = plot(real(xr),real(yr),'color',mousept.color);


