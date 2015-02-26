function t = TimeRange(Ch, varargin)
t = [];


if isfield(Ch,'start')
    t(1) = Ch.start;
    t(2) = Ch.start + Ch.interval .* Ch.length;
end