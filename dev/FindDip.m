function [dip, details] = FindDip(v, varargin)
%Find a dip in a histogram of x;

sigma = std(v);
sm = sigma/10;
r = diff(minmax(v));
x = min(v) - r/10:r/100:max(v)+r/10;
[y,x] = smhist(v,'sd',sm,'xval',x);
peaks = find(diff(sign(diff(y))) < 0);
while length(peaks) < 2
    sm = sm ./ 2;
    [y,x] = smhist(v,'sd',sm,'xval',x);
    peaks = find(diff(sign(diff(y))) < 0);
end
while length(peaks) > 2
    sm = sm .* 1.5;
    [y,x] = smhist(v,'sd',sm);
    peaks = find(diff(sign(diff(y))) < 0);
end
while length(peaks) < 2
    sm = sm * 0.9;
    [y,x] = smhist(v,'sd',sm);
    peaks = find(diff(sign(diff(y))) < 0);
end
peaks = peaks([1 end]);
dip = x(peaks);
rawx = x;
rawy = y;
    
a = find(y > max(y)/40);
id = a(1):a(end);
dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
while length(dpeaks) > 4
    sm = sm  .* 2;
    y = smhist(v,'sd',sm,'xval',x);
    dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
end
signpeaks = diff(sign(diff(y(id),2)));
signpeaks = signpeaks(dpeaks);

    
lastbig = 0;
mx = 0.8;
j = 0;
while length(dpeaks) ~= 3 && length(dpeaks) ~= 4
    j = j+1;
    if length(dpeaks) > 3
        sm = sm ./ mx;
        big = 1;
    else
        sm = sm.*mx; %less smoothing
        big = 0;
    end
    sms(j) = sm;
    if lastbig ~=big;
        mx = sqrt(mx);
    end
    lastbig = big;
    y = smhist(v,'sd',sm,'xval',x);
    dpeaks = find(abs(diff(sign(diff(y(id),2)))) > 0);
    signpeaks = diff(sign(diff(y(id),2)));
    signpeaks = signpeaks(dpeaks);
    np(j) = length(dpeaks);
    %        fprintf('sm %.3f, %d peaks\n',sm,length(dpeaks));
end
fy = smhist(v,'sd',sm/2,'xval',x);

dy = diff(fy(id));
[a,b] = max(y);  %1 or 3
ns = sum(x(id(dpeaks)) < x(b));
if length(dpeaks) == 3
    if signpeaks(1) > 0
        [a,b] = min(abs(dy(dpeaks(1):dpeaks(2))));
        dip(4) = x(id(dpeaks(1)+b-1));
    else
        [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
        dip(4) = x(id(dpeaks(2)+b-1));
    end
elseif ns > 2
    [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
    dip(4) = x(id(dpeaks(2)+b-1));
else
    [a,b] = min(abs(dy(dpeaks(2):dpeaks(3))));
    dip(4) = x(id(dpeaks(2)+b-1));
end
details.fy = fy;
details.x = x;