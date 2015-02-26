function [m, nt] = WeightedSum(x, n, varargin)
%
% [m, n] = WeightedSum(x, n, varargin)
% calculates the Weighted mean of x, wieghted by
% n. x can be a 2-D matrix.
dosqrt = [0 0];
j = 1;
while j <= nargin-2
    if strncmpi(varargin{j},'sqrt',4)
        dosqrt = [1 0];
    end
        
    j = j+1;
end
if dosqrt(1)
    x = x.^0.5;
end
if dosqrt(2)
    n = n.^0.5;
end

nt = sum(n);
if iscell(x)
    nt = sum(nt(:));
    m = zeros(size(x{1}));
    for k = 1:size(x,1)
    for j = 1:size(x,2)
        if n(j) > 0 & length(x{k,j})
            m = m + x{k,j}.*n(k,j);
        end
    end
    end
    m = m./nt;
elseif size(x,2) == size(n,1)
x(isnan(x)) = 0;  %NaNs have n = 0
    m = (x * n)./ sum(n);
elseif size(x,2) == size(n,2) & size(x,2) > 1
x(isnan(x)) = 0;  %NaNs have n = 0
    m = (x * n')./ sum(n);
elseif size(x,1) == size(n,2)
    n = repmat(n',1,size(x,2));
    id = find(isnan(x));
    x(id) = 0;
    n(id) = 0;
    m = sum(x.*n)./sum(n);
elseif length(x) == length(n)
    id = find(~isnan(x));
    m = sum(x(id) .* n(id))./sum(n(id));
elseif size(x,1) == length(n)
    id = find(~isnan(sum(x,2)));
    m = (n(id)' * x(id,:))./sum(n(id));   
else
    m = ones(size(x)) .* NaN;
end

if dosqrt(1)
    m = m.^2;
end