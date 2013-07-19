x = -pi:pi/100:pi;

both = [];
for j = 1:100;
a = sin(x);
both(j,:) = a+ sin(2 * x +x(j));
end

%find point with max peak-peak amplitude.
[maxa, besti] = max(max(b,[],2) - min(bin,[],2));
