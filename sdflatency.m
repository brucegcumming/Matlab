function [latency, details] = sdflatency(sdf, presamples, varargin)

%[latency, details] = sdflatency(sdf, presamples, varargin) estmate latency
% of sdf. 
% returns a latency estimate from a spike density function, sdf
% presamples is the number of samples in the sdf guaranteed not to
% contain a response.
% offset, if defined, determines a starting location in the buffer


baseline = NaN;
offset = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'baseline',4)
        j = j+1;
        baseline = varargin{j};
    elseif strncmpi(varargin{j},'offset',4)
        j = j+1;
        offset = varargin{j};
    end
    j = j+1;
end


if presamples+offset > length(sdf)
    latency = NaN;
    return;
end

if size(sdf,2) == 1
    sdf = sdf';
end

t(1) = presamples+offset;
prerate = mean(sdf(1+offset:presamples+offset));

[maxv, maxt] = max(sdf);
halfv = prerate + (maxv-prerate)/2;
id = find(sdf(1:maxt) < halfv);
halft = id(end);
v = prerate + (maxv-prerate)/5;
id = find(sdf(1:halft) < v);
vt = id(end);

%guess latency
x(1) = vt - (halft-vt)/1.5;
if x(1) < t(1)
    x(1) = t(1);
end
x(2) = (halfv-v)./(halft-vt); %guess slope

if isnan(baseline)
    x(3) = prerate;
else
    state.baseline = baseline;
end
if halft < t(1)+presamples/10;
    halft = t(1)+presamples/2;
end

state.sdf = sdf(1:halft);
options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
[fittedparams, fval,exitflag, output] = fminsearch(@FitStart, x,options, state);
[a, b] = FitStart(fittedparams,state);
if ~isempty(b)
details.fitted = b.fit;
end
latency = x(1);



function [ssq, details] = FitStart(x, state)

details = [];
x(1) = round(x(1));
if x(1) >= length(state.sdf) || x(1) < 1
    ssq = NaN;
    return;
end
if length(x) == 3
    base = x(3);
else
    base = state.baseline;
end
pred(1:x(1)) = base;
pred(x(1):length(state.sdf)) = base + (x(2) .* 0:length(state.sdf)-x(1));
diffs = state.sdf-pred;
ssq = sum(diffs.^2);
details.fit = pred;

