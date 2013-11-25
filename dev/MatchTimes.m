function ids = MatchTimes(ta, tb, varargin)

ids = [];
method = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'method',5)
        j = j+1;
        method = varargin{j};
        fprintf('Using method %d ', method');
    elseif strncmpi(varargin{j},'selftest',5)
        for k = 1:6
            ids = MatchTimes(ta, tb, 'method', k);
        end
        return;
    end
    j = j+1;
end
tic;
ids = zeros(size(ta));
if method == 1
    for j = 1:length(ta)
        x = find(tb > ta(j),1);
        if isempty(x)
            ids(j) = NaN;
        else
        ids(j) = x;
        end
    end
elseif method == 2
    for j = 1:length(ta)
        a = find(tb > ta(j));
        if ~isempty(a)
            ids(j)  = a(1);
        end
    end
elseif method == 3
    [x,y] = meshgrid(ta,tb);
    [b, ids] = min(abs(x-y));
elseif method == 4
    x = repmat(ta, 1, length(tb));
    y = repmat(tb', length(ta),1);
    [b, ids] = min(abs(x-y));
elseif method == 5
    last = 1;
    for j = 1:length(ta)
        x = find(tb(last:end) > ta(j),1);
        if isempty(x)
            ids(j) = NaN;
        else
            ids(j) = x+last-1;
            last = ids(j);
        end
    end
elseif method == 6
    for j = 1:length(ta)
        x = find(tb > ta(j));
        if isempty(x)
            ids(j) = NaN;
        else
        ids(j) = x(1);
        end
    end
end
toc