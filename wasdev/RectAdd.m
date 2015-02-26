function M = RectAdd(a,b, varargin)
%M = RectAdd(a,b,...) creates a matrix M from adding
%the colomns (rows) of a to the rows(columns) of b;
%RectAdd(a,b,'range',[min max]) makes sure the elements of M
% are >= min and <= max;

range = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'range',5)
        j = j+1;
        range = varargin{j};
    end
    j = j+1;
end

if size(a,1) == 1 && size(b,2) == 1
    M = repmat(a,size(b,1),1) + repmat(b,1,size(a,2));
elseif size(a,2) == 1 && size(b,1) == 1
    M = repmat(a,size(b,1),1) + repmat(b,1,size(a,2));
elseif size(a,1) == 1 && size(b,1) == 1
    M = repmat(a,size(b,2),1) + repmat(b',1,size(a,2));
elseif size(a,2) == 1 && size(b,2) == 1
    M = repmat(a,size(b,1),1) + repmat(b',1,size(a,1));
end

if length(range) == 2
    M(M < range(1)) = range(1);
    M(M > range(2)) = range(2);
end
