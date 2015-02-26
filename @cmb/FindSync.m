function [aid, bid, n] = FindSync(at,bt,dt,varargin)
%from two lists of times, find events synchronous witthin dt
n = 0;
k = 1;
aid = [];
bid = [];
if isempty(at) || isempty(bt)
return;
end
if length(bt) < 101
bts = [1:length(bt)];
else
bts = [1:101];
end
nb = length(bt);
for j = 1:length(at)
t = bt(bts)-at(j);
while(t(end) < 0 & bts(end) < nb)
if bts(end)+100 > nb
bts = bts(end):nb;
else
bts = bts+100;
end
t = bt(bts)-at(j);
end
[a,b] = min(abs(t));
if a <= dt
n = n+1;
aid(n) = j;
bid(n) = bts(b);
end
end

