function [count, vals] = Counts(x, varargin)
%[counts, vals] = Counts(x) 
%return counts for each unique value of x
%[counts, vals] = Counts(x, vals)
%return counts for each value in vals. N.B. not a histogram, 
% only counts exact matches

minval= NaN;
if isempty(x)
    count = [];
    vals  = NaN;
    return;
end
sorttype = '';

j = 1;
while j <= length(varargin)
    if sum(strncmpi(varargin{j},{'ascend' 'descend'},5))
        sorttype = varargin{j};
    elseif strncmpi(varargin{j},'minval',6)
        j = j+1;
        minval = varargin{j};
    end
    j = j+1;
end

if length(varargin) & isnumeric(varargin{1})
    vals = varargin{1};
else
vals = unique(x);
end

if iscellstr(x)
    for j = 1:length(vals)
        count(j) = sum(strcmp(vals{j},x(:)));
    end
else
    for j = 1:length(vals)
        count(j) = sum(x(:) == vals(j));
    end
end

if ~isnan(minval);
    count = count(vals >= minval);
    vals = vals(vals >= minval);
end

if ~isempty(sorttype)
    [count,b] = sort(count, sorttype);
    vals = vals(b);
end