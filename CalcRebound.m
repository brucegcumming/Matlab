function [sz, lat, peak] = CalcRebound(X, varargin)


if min(size(X) == 1)
    [sz,lat, peak] = Rebound(X);
    return;
end
    
for j = 1:size(X,1)
    [sz(j),lat(j), peak(j)] = Rebound(X(j,:));
end


function [sz, lat, peak] = Rebound(X)

[a,b] = min(X);
peak = max(X(b:end))./mean(X);
base = mean(X(1:20));
id = find(X(b:end) > base);
if isempty(id)
    sz = 0;
    lat = 0;
else
    lat = b+id(1)-1;
    sz = sum(X(lat:end)-base)./mean(X);
end