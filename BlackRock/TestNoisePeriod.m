function X = TestNoisePeriod(FullV, varargin)

block =1;
pts = [];
nperiods = 10;
j = 1;
while j<= length(varargin)
    if strncmpi(varargin{j},'block',5)
        j = j+1;
        block = varargin{j};
    elseif strncmpi(varargin{j},'pts',3)
        j = j+1;
        pts = varargin{j};
    end
    j = j+1;
end

if ~isempty(pts)
elseif block > 1
    a = cumsum(FullV.blklen(1:block));
    pts = a:a+FullV.blklen(block)-1;
else
pts = 1:length(FullV.V);
pts = 1:FullV.blklen(1);
end

%pts = pts+length(pts)*offset;
periods = (4.942*1):0.0002:(4.947*1);
periods = periods * nperiods;
for p = 1:length(periods);
    period = periods(p);
    x = round(1 * mod(pts./period,1).*period);
    bins = unique(x);
    for j = 1:length(bins)
        id = find(x == bins(j));
        avg(j) = mean(FullV.V(id));
    end
    amp(p) = std(avg);
    avgs(p,1:length(avg)) = avg;
end
imagesc(avgs);
X.periods = periods;
X.amps = amp;
X.avgs =avgs;
