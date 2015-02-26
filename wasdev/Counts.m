function [count, vals] = Counts(x, varargin)
%[counts, vals] = Counts(x) 
%return counts for each unique value of x
%[counts, vals] = Counts(x, vals)
%return counts for each value in vals. N.B. not a histogram, 
% only counts exact matches
% if x is a cell array of strings, vals is returns a cell array
% if x is a cell array of numeric values, count returns a vector, with counts
%  counts(1) is the number of unique elements. counts(2:n+1) is the count of unique values of
%  the nth element of each element of x. 
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
elseif iscellstr(x) 
    vals = unique(x);
elseif iscell(x)
    for j = 1:length(x)
        lens(j) = length(x{j});
        for k = 1:length(x{j})
            allvals(j,k) = x{j}(k);
        end
    end
    count(1) = length(unique(lens));
    if max(lens) > 1
    for j = 1:size(allvals,2);
    count(j+1) = length(unique(allvals(:,j)));
    vals{j} = unique(unique(allvals(:,j)));
    end
    else
        vals = unique(cat(1,x{:}));
    end
else
    vals = unique(x);
end

if iscellstr(x)
    for j = 1:length(vals)
        count(j) = sum(strcmp(vals{j},x(:)));
    end
elseif iscell(x)
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