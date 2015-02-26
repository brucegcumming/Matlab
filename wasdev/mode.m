function m  = mode(x)
%deprecated. Use prctile(x,50);
if isempty(x)
    m = [];
    return;
end
vals = unique(x(find(~isnan(x))));
j = 1;
n = [];
%for val = vals
for j = 1:length(vals)
  n(j) = length(find(x == vals(j)));
end

if isempty(n) %x was all NaNs
    m = NaN;
else
[nmax, imax] = max(n);
m = vals(imax);
end


