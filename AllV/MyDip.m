function result = MyDip(values, varargin)
%X = MyDip(X, varargin) Heuristic for finding best dip in a distribution
%if X is a cell array of Gaussian Fits, find the dip

if iscell(values) && isfield(values{1},'amp')
    result = GaussDip(values);
    return;
end

ids = [];
cl = 2;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'idlist',4)
        j = j+1;
        ids = varargin{j};
    end
    j = j+1;
end
v = sort(values);
if ~isempty(ids)  
    l(1) = mean(values(ids ==cl));
    l(2) = mean(values(ids ==1));
    if diff(l) < 0
        l = fliplr(l);
    end
    x = linspace(l(1),l(2));
    id = find(v > l(1) & v < l(2));
    sm = prctile(diff(v(id)),10).*200;
    y = smhist(v,'sd',sm,'xval',x);
    minima = find(diff(sign(diff(y))) >0);
    while length(minima) > 1
        sm = sm*2;
        y = smhist(v,'sd',sm,'xval',x);
        minima = find(diff(sign(diff(y))) >0);
        while isempty(minima)
            sm = sm.*0.9;
            y = smhist(v,'sd',sm,'xval',x);
            minima = find(diff(sign(diff(y))) >0);
            if length(minima) > 1
                minima = minima(1);
            end
        end
    end
    result.x = x;
    result.y = y;
    dip(1) = minima(1);
    dip(2) = y(minima(1));
    dip(3) = max(y(1:dip(1)));
    dip(4) = max(y(dip(1):end));
    result.dip = dip;
    result.d = (x(dip(1))-x(1))./diff(minmax(x));

    return;
end

sigma = std(v);
sm = sigma/10;
minpeakratio = 20;
r = range(values);
x = min(values) - r/10:r/100:max(values)+r/10;
[y,x] = smhist(v,'sd',sm,'xval',x);
peak = max(y);
peaks = find(diff(sign(diff(y))) < 0 & y(2:end-1) > peak/minpeakratio);
xx = 2;
up = 0;
while length(peaks) ~= 2 && xx > 1.01
    if length(peaks) < 2
        sm = sm ./ xx;
        if up == 2
            xx = sqrt(xx);
        end
        up = 1;
    else
        sm = sm .* xx;
        if up == 1
            xx = sqrt(xx);
        end
        up = 2;
    end
    [y,x] = smhist(v,'sd',sm,'xval',x);
    peak = max(y);
    peaks = find(diff(sign(diff(y))) < 0 & y(2:end-1) > peak/minpeakratio);
end
result.peaks = peaks([1 end]);
[a,b] = min(y(peaks(1):peaks(end)));
result.dip(1) = peaks(1)+b-1;
result.dipsize(1) = y(result.dip(1))./min(y(result.peaks));
finey = smhist(v,'sd',sm/2,'xval',x);
[a,b] = min(finey(peaks(1):peaks(end)));
result.dip(2) = peaks(1)+b-1;
result.dipsize(2) = finey(result.dip(2))./min(finey(result.peaks));
finey = smhist(v,'sd',sm/4,'xval',x);
result.finey = finey;

[a,b] = min(finey(peaks(1):peaks(end)));
result.dip(3) = peaks(1)+b-1;
result.dipsize(3) = finey(result.dip(3))./min(finey(result.peaks));
dy = diff(y);
[a,b] = min(dy(peaks(1):peaks(end)-1));
id = b+peaks(1)-1; %min in dy after first peak
xid = find(diff(sign(diff(dy)))<0)+1;
if xid(end) < peaks(1)+5
    result.dip(4) = xid(end);
else
    xid = xid(xid >= peaks(1)+5);
    accel = diff(diff(dy));
    [a,b] = min(accel(xid-1));
    result.dip(4) = xid(b);
end
result.x = x;
result.y = y;

function result = GaussDip(fits, varargin)


xr = cat(2,fits{1}.state.meanlimit,fits{2}.state.meanlimit);
fx = linspace(min(xr(:)),max(xr(:)),100);
         fya = FitGauss(fx, fits{1}.params, 'eval');
         fyb = FitGauss(fx, fits{2}.params, 'eval');
         fy = fya+fyb;
         id = find(fx > fits{1}.params(1) & fx < fits{2}.params(1));
         if ~isempty(id)
             [a,b] = min(fy(id));
             minxpt = fx(id(b));
             result.dip(1) = id(b);
         else
             result.dip(1) = NaN;
             minxpt = NaN;
         end
result.x = fx;
result.y = fy;