function DrawEllipses(a,b,varargin)


xl = get(gca,'xlim');
xrange = range(xl);
yl = get(gca,'ylim');
yrange = range(yl);
c = [mean(xl) mean(yl)];
a = a./xrange;
b = b./yrange;
hold off; 
for angle = 0:pi/8:pi;
sn = sin(angle);
cn = cos(angle);
%sn = sin(pi/4);
%cn = cos(pi/4);
x = linspace(0,a);
y =  sqrt(b.^2 - (x.*b/a).^2);
x = [x fliplr(x) -x fliplr(-x)];
y = [y fliplr(-y) -y fliplr(y)];
xr = xrange .* (x .* cn + y .*sn) + c(1);
yr = yrange .* (y .* cn - x .*sn) + c(2);
plot(real(xr),real(yr));
hold on;
end
