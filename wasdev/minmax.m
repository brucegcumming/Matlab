function [mm, id] = minmax(x, varargin)

if isempty(x)
    mm = [];
    id = [];
    return;
end
[mm(2) id(2)] = max(x(:));
[mm(1) id(1)] = min(x(:));