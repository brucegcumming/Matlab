function mm = minmax(x, varargin)

if isempty(x)
    mm = [];
    return;
end
mm(2) = max(x(:));
mm(1) = min(x(:));