function result = MyDip(values, varargin)
%X = MyDip(X, varargin) Heuristic for finding best dip in a distribution
%if X is a cell array of Gaussian Fits, find the dip
%result.dipsize = min at dip / smallest peak; 
%result.dipsize(1,2,3) are for less and less smoothed histograms.

if iscell(values) && isfield(values{1},'amp')
    result = GaussDip(values);
    return;
end

quickmode = 0;
ids = [];
cl = 2;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'idlist',4)
        j = j+1;
        ids = varargin{j};
    elseif strncmpi(varargin{j},'quick',4)
        quickmode = 1;
    end
    j = j+1;
end

if isempty(values)
    result.dip = 1;
    result.x = 0;
    return;
end
if quickmode
    result = QuickDip(values, varargin{:});
    return;
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

function result = QuickDip(values,varargin)
%do a quick version for rapid autocutting, esp online spikes.

v = sort(values);

result.sigma = std(v);
result.nvals = length(v);
minpeakratio = 20;
sm = 2 * result.sigma/sqrt(length(v));
result.sm(1) = sm;
r = range(values);
x = min(values) - r/10:r/100:max(values)+r/10;
if length(values) < 3
    peak = 1;
    result.x = v;
    result.x = ones(size(v));
    result.dip(1) = 1;
    result.peaks = v;
    return;
end
    
[y,x] = smhist(v,'sd',sm,'xval',x);
peak = max(y);
peaks = find(diff(sign(diff(y))) < 0 & y(2:end-1) > peak/minpeakratio);
xx = 1.5;
up = 0;
while length(peaks) > 10
    sm = sm .* 2;
    result.sm(end+1) = sm;
    [y,x] = smhist(v,'sd',sm,'xval',x);
    peak = max(y);
    peaks = find(diff(sign(diff(y))) < 0 & y(2:end-1) > peak/minpeakratio);
end
while length(peaks) < 2
    sm = sm ./ xx;
    result.sm(end+1) = sm;
    [y,x] = smhist(v,'sd',sm,'xval',x);
    peak = max(y);
    peaks = find(diff(sign(diff(y))) < 0 & y(2:end-1) > peak/minpeakratio);
end
valleys = find(diff(sign(diff(y))) > 0);
valleys = valleys+1;
peaks = peaks+1;
id = find(valleys > peaks(1) & valleys < peaks(end));
valleys = valleys(id);
[a,b] = min(y(valleys));
if (b >1 & b < length(valleys)) || length(valleys) < 2
    result.dip(1) = valleys(b);
    p = min(y(peaks));
    dipsize(1) = y(valleys(b))./p;
    valleypts = 1:length(valleys);
else
    k =1;
    for j = 1:length(peaks)-1
        if valleys(k) < peaks(j)
            k = k+1;
        end
        p(j) = min([max(y(peaks(1:j))) max(y(peaks(j+1:end)))]);
        dipsize(j) = y(valleys(k))./p(j);
        valleypts(j) = k;
    end
end
    

%at least two qualifying peaks here.
result.peaks = peaks;
[a,b] = min(dipsize);
result.dip(1) = valleys(valleypts(b));
result.dipsize(1) = a;


result.x = x;
result.y = y;


