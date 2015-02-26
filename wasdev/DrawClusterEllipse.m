function h= DrawEllipse(E,varargin)

if E.shape == 1
    h = DrawLine(E,varargin{:});
    return;
end
if isfield(E,'pos')
a = (E.pos(3)-E.pos(1))/2; %radius 
b = (E.pos(4)-E.pos(2))/2;
cntr(1) = mean(E.pos([1 3]));
cntr(2) = mean(E.pos([2 4]));
elseif isfield(E,'xyr')
    a = E.xyr(3);
    b = E.xyr(4);
    cntr(1) = E.xyr(1);
    cntr(2) = E.xyr(2);
end
    
sn = 0;
cn = 1;

x = linspace(0,a);
sn = sin(E.angle);
cn = cos(E.angle);
if isfield(E,'aspectratio') && E.aspectratio > 0
    b  = b ./E. aspectratio;
    y =  sqrt(b.^2 - (x.*b/a).^2);
else
    y =  sqrt(b.^2 - (x.*b/a).^2);
end

x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+cntr(1);
if isfield(E,'aspectratio') && E.aspectratio > 0
    y = yr.*E.aspectratio+cntr(2);
else
    y = yr+cntr(2);
end

if ~isfield(E,'color')
    E.color = 'r';
end

if isfield(E,'h') & ishandle(E.h)  & get(E.h,'parent') == gca
%    fprintf('Using handle %f\n',get(E.h,'parent'));
    set(E.h,'Xdata',x,'Ydata',y,'color',E.color);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),'color',E.color,varargin{:});
    hold off;
end

function h= DrawLine(E,varargin)

x = [E.pos(1) E.pos(3)];
y = [E.pos(2) E.pos(4)];
%fprintf('%.3f,%.3f -> %.3f,%.3f\n',x(1),y(1),x(2),y(2));
if isfield(E,'h') & ishandle(E.h) & get(E.h,'parent') == E.axis; 
    set(E.h,'Xdata',x,'Ydata',y);
    h = E.h;
else
    hold on;
    h = plot(real(x),real(y),varargin{:});
    hold off;
end
