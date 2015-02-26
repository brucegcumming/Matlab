function y = expand(z, scale, midpt)

x = [1:length(z)] -midpt;
xi = min(x):diff(minmax(x))/1000:max(x);
newx = x ./ scale;
xi = MatchInd(newx, xi, 'nearest');
y = interp1(x,z,xi);