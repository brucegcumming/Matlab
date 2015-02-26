function tr = TimeRange(x)

tr = [];
if isfield(x,'interval') && isfield(x,'length') %A spike2 Channel
    tr(1) = x.start;
    tr(2) = x.start + x.length.*x.interval;
end