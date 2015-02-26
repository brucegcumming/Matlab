function d = timediff(t, varargin)
%timediff(t, varargin)
%calculates differences between time values, in seconds (default)
%if t is a vector of times, returns differeces relative to earliest

reftime = [];
units = 'sec';
if isempty(reftime)
    reftime = min(t);
end
if strcmp(units,'sec')
d = (t-reftime) .* 24 * 60 * 60;
end