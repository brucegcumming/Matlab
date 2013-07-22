function [x,y] = TiltElectrode(Z, L, a, b)
%TiltElectrode(Z, L, a, b)
ca = cos(a.*pi/180);
ta = tan(a.*pi/180);
cb = cos(b.*pi/180);
tb = tan(b.*pi/180);
sa = sin(a.*pi/180);
sb = sin(b.*pi/180);
plottype = 1;
x = Z.*ca.*cb + L.*sa.*cb - L .* ca .* sb + Z .* sa .*sb;
y = (L.*ca  - Z .* sa + x .* sb)./cb;

if plottype == 1
plot([0 Z.* ca Z.*ca+L.*sa],[0 -Z .* sa L .* ca - Z.*sa]);
hold on;
plot([0 x.* cb x.*cb+y.*sb],[0 -x .* sb y .* cb - x.*sb],'r');
set(gca,'ydir','reverse');
axis image;
else
   x = [Z x];
   r = 11.5;
   h(1) = (r .* ca -x(1))./ta;
   h(2) = (r .* cb -x(2))./tb;
   plot([r r],[0 30]);
   hold on;
   plot([0 r],[0 r .* -tan(b.*pi/180)],'r');
   plot([0 r],[0 r .* -tan(a.*pi/180)]);
   plot([x(1).* ca r],[-x(1) .* sa h(1)]);
   plot([x(2).* cb r],[-x(2) .* sb h(2)],'r');
set(gca,'ydir','reverse');
axis image;
end