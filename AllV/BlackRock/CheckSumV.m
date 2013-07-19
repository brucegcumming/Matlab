function CheckSumV(sumv, varargin)

FullV = [];
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'intscale')
        FullV = varargin{j};
    end
    j = j+1;
end

if ~isempty(FullV)
    sumsq = sumv*sumv';
    FullV.V = double(FullV.V);
    plot(FullV.V);
    FullV.sumscale = (sumv*FullV.V')./sumsq;
    FullV.V = FullV.V - sumv .* FullV.sumscale;
    hold on;
    plot(FullV.V,'r');
    return;
end

sd = std(sumv);
sgn = diff(sign(diff(sumv)));
peaks = sumv(1+find(sgn ~= 0));
th = sd*3:sd/2:sd*10;
for j= 1:length(th)
    y(j) = sum(abs(peaks) > th(j));
end
plot(th,y,'o-');  
%hist(peaks,1000);