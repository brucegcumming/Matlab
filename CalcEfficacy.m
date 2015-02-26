function [e, nspk, details] = CalcEfficacy(x,y,varargin)
% [e, nspk] = CalcEfficacy(x,y,varargin)
%find peak in xcorrelogram over +- 1ms range.
%then count sync spikes +- 0.2ms from this
%e(1) is %of spikes in y within 0.2ms of x
%e(2) is %of spikes in x within 0.2ms of y
%e(1:2) con-incidence #s are adjusted by baseline rate.
%e(3) and e(4) are the same calculation for spikes <1ms apart but more than
%0.2ms. So if e(3) is as big as e(1) its not true synchrony
method = 1;
details = {};
rate(1) = median(diff(x));
rate(2) = median(diff(y));
p = 1./(rate * 1000); %p(spike) in 1ms bin
nspk = [length(x) length(y)];
adjust(1)  = length(x) .* p(2);
adjust(2)  = length(y) .* p(1);

if isempty(x) || isempty(y)
    e = [0 0 0 0];
    return;
end
if method == 1 %slightly quicker

    a = bsxfun(@minus,x,y');
    xi = find(abs(a) < 0.001);
    adjust = adjust * 2; %+-1ms is a 2ms bin
    delay = median(a(xi));
    if nargout > 2
        [details.id{2}, details.id{1}] = find(abs(a - delay) < 0.0002);
        nsync = length(details.id);
    else
        nsync = sum(abs(a(:) - delay) < 0.0002);
    end
    nsyncw = length(xi)-nsync;
    if length(xi) == nsync %all spikes in central ms are in same bin
        if nsync == 1
            nsyncw = 1;
        else
            xc = histc(a(:)-delay,[-0.01:0.0002:0.01]);
            dx = sort(xc);
            nsyncw = dx(end-1); %second largest peak. 
        end
    end
else
    X = repmat(x(:),1,length(y));
    Y = repmat(y(:)',length(x),1);
    d = X-Y;
    nsync = sum(abs(d(:)) < 0.001);
end
nsync = nsync-(adjust/5);
e = nsync./[length(x)-adjust(1)/5 length(y)-adjust(2)/5];
e([3 4]) = nsyncw./[length(x)-adjust(1).*0.8  length(y)-adjust(2).*0.8];

function t = absdiff(a,b)
t = abs(a-b) < 0.001;

