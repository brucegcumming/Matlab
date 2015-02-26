function [xr, yr] = ClusterRange(C, p)

for j = 1:size(C,1)
if isfield(C{j,p},'x') & C{j,p}.x(2) > 0
xh(j) = C{j,p}.x(1)+C{j,p}.x(2);
xl(j) = C{j,p}.x(1)-C{j,p}.x(2);
yh(j) = C{j,p}.y(1)+C{j,p}.y(2);
yl(j) = C{j,p}.y(1)-C{j,p}.y(2);
else
xh(j) = 0;
xl(j) = 0;
yl(j) = 0;
yh(j) = 0;
end
end
xr = [min(xl) max(xh)];
yr = [min(yl) max(yh)];

