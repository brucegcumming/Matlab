function ids = MatchTimes(ta, tb, varargin)

ids = [];
method = 1;

tic;
if method == 1
    for j = 1:length(ta)
        ids(j) = find(tb > ta(j),1);
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
end
toc