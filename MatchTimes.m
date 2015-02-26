function [ids, xids, taid, nearest] = MatchTimes(ta, tb, varargin)
%ids = MatchTimes(ta, tb, tw) find times in tb that match ta within +- dw
% [ids, missed] = MatchTimes(ta, tb, tw) also returns elements in tb with
% no match
% [ids, missed, ia] = MatchTimes(ta, tb, tw) returns indices of ta that
% match
% [ids, missed, ia, nearest] = MatchTimes(ta, tb, tw) returns nearest
% matches in ta to th enon-matches in tb



ids = [];
xids = [];
taid = [];
nearest = [];

method = 1;
j = 1;
if ~isempty(varargin) && isnumeric(varargin{1})
    tw = varargin{1};
    j = 2;
else
    tw = 1;
end
while j <= length(varargin)
    if strncmpi(varargin{j},'method',5)
        j = j+1;
        method = varargin{j};
        fprintf('Using method %d ', method');
    elseif strncmpi(varargin{j},'selftest',5)
        for k = 1:7
            ids = MatchTimes(ta, tb, tw, 'method', k);
        end
        return;
    end
    j = j+1;
end
%tic;
ids = zeros(size(ta));
if isempty(ta) || isempty(tb)
    return;
end
if method == 7
    d = bsxfun(@minus,ta(:),tb(:)');
    dt = min(abs(d));
    ids = find(dt < tw); %match within 1 sec
elseif method == 1
    for j = length(ta):-1:1;
        x = find(tb > ta(j)-tw,2);
        if isempty(x)
          ds(j) = NaN;
        else
          [ds(j),a]  = min(abs(tb(x)-ta(j)));
          bid(j) = x(a);
        end
    end
    taid = find(ds < tw);
    ids = bid(taid);
    [ids, ia] = unique(ids);
    if nargout > 2
        taid = taid(ia);
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
    x = repmat(ta(:), 1, length(tb));
    y = repmat(tb(:)', length(ta),1);
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

if nargout > 1
    xids = setdiff(1:length(tb),ids);
end
if nargout > 3
    d = bsxfun(@minus,ta(:),tb(xids));
    [dt, nearest] = min(abs(d));
end
%toc;