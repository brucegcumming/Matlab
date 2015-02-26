function [h, y] = DrawEllipse(xyr, varargin)
%DrawEllipse(xyr, varargin) draws ellipse
%at xyr(1),xyr(2). Radii xyr(3),xyr(4) rotated through angle xyr(5);
%if xyr has length 3, draws circle radius xyr(3)
%[x, y] = DrawEllipse(xyr, 'noplot') returns a set of points, but does no
%plotting

onoff = {'off' 'on'};
noplot = 0;
holdaxis = 1; 
if length(xyr) < 4
    xyr(4) = xyr(3); %circle
end
if length(xyr) < 5
    xyr(5) = 0; 
end
holdstate = ishold;
skip = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'add',3)
        hold on;
        skip = [skip j];
    elseif strncmpi(varargin{j},'noplot',5)
        noplot = 1;
    end
    j = j+1;
end
varargin = varargin(setdiff(1:length(varargin),skip));
aspectratio = 1;
a = xyr(3); %x radius in nomalized units
b = xyr(4);
angle = xyr(5);
sn = 0;
cn = 1;
x = linspace(0,a);
b  = b ./aspectratio;
y =  sqrt(b.^2 - (x.*b/a).^2);

sn = sin(angle);
cn = cos(angle);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = (x .* cn + y .*sn);
yr = (y .* cn - x .*sn);
x = xr+mean(xyr(1));
y = yr.*aspectratio+mean(xyr(2));
h = x;
if noplot
    return;
end
xl = get(gca,'xlim');
yl = get(gca,'ylim');
h = plot(real(x),real(y),varargin{:});    
if holdaxis
    set(gca,'xlim',xl,'ylim',yl);
end
hold(onoff{holdstate});