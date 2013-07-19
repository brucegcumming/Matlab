function [maxs, mins] = Extrema2d(Z, varargin)
%[maxs, mins] = Extrema2d(Z, ...)
%finds local maxima/minima for 2-d function
% Extrema2d(Z, 'min', [a b]) sets threhsoding. only minima < a, and maxima > b, count 
zmin = min(Z(:));
zmax = max(Z(:));

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'min',3)
        j = j+1;
        zmax = varargin{j}(1);
        if length(varargin{j}) > 1
            zmin = varargin{j}(2);
        end
    end
    j = j+1;
end

lx = size(Z,1);
ly = size(Z,2);
dy = diff(Z,[],2);
dx = diff(Z);
yz = diff(sign(dy),[],2);
xz = diff(sign(dx));
maxid = find(yz(2:lx-1,:) < 0 & xz(1:198,2:ly-1) < 0 & Z(2:lx-1,2:ly-1) > zmin);
[py,  px] = ind2sub([lx-2 ly-2], maxid);
px = px+1;
py = py+1;
maxs = [px py];

minid = find(yz(2:lx-1,:) > 0 & xz(1:198,2:ly-1) > 0 & Z(2:lx-1,2:ly-1) < zmax);
[py,  px] = ind2sub([lx-2 ly-2], minid);
px = px+1;
py = py+1;
mins = [px py];

