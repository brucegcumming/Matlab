function screen(vscale, hscale)


arrowplot =1;
[x,y] = meshgrid(-1000:200:1000,-1000:200:1000);
xa = x + x .* y .* (1-vscale(1));
ya = y + y .* y .* (1-vscale(2));
xb = x + x .* y .* (1-1/vscale(1));
yb = y + y .* y .* (1-1/vscale(2));
xc = xa + x .* x .* (1-1/hscale(1));
yc = ya + y .* x .* (1-1/hscale(2));
hold off;
plot(x,y,'sk');
hold on;
plot(xc,yc,'sr');

xd = xc + xc .* yc .* (1-1/vscale(1));
yd = yc + yc .* yc .* (1-1/vscale(2));
%plot(xd,yd,'sg');

if arrowplot
    hold off;
    for j = 1:prod(size(x))
        hl = max([50 sqrt((x(j)-xc(j)).^2 + (x(j)-xc(j)).^2)/5]);
        arrow([x(j) xc(j)], [y(j) yc(j)], 30, hl);
        hold on;
    end
end

title(sprintf('Tz = %.7f %.7f %.7f %.7f',vscale(1),vscale(2),hscale(1),hscale(2)));