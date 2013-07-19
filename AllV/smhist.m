function [smoothed, x] = smhist(dat, varargin)
%
% [smoothed, x] = smhist(dat, ...)
%
% returns a smoothed density histogram for a set of values.
% x returns the values of dat associated with each density value.
%
% the raw data are smoothed with a Gaussian whose values are estimated
% based on the data.
% smhist(x,'smoother', z)
% Changes the sigma to is originalSD/z
%
%  smhist(x, 'sd', sigma) uses a SD of sigma for smoothing.
%  (...'nsmp', n) returns a function with n samples
%  (...'xval', X) returns a function evaluated at X (vector)
%  
nsmp = 100;
dat = dat(find(~isnan(dat)));

%sd = 4 * range(dat)/sqrt(length(dat));
sd = std(dat)/sqrt(length(dat));
step = (2 * sd + diff(minmax((dat))))/nsmp;
drange = [];
period = 0;
setx = [];
smoothed = [];
x = [];
j= 1;
while j < nargin
    if strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    elseif strncmpi(varargin{j},'circular',3)
        period = 2 * pi;
        dat = [dat; dat + period];
    elseif strncmpi(varargin{j},'period',3)
        j = j+1;
        period = varargin{j};
        dat = [dat; dat + period];
    elseif strncmpi(varargin{j},'nsmp',2)
        j = j+1;
        nsmp = varargin{j};
    elseif strncmpi(varargin{j},'range',2)
        j = j+1;
        drange = varargin{j};
    elseif strncmpi(varargin{j},'step',2)
        j = j+1;
        step = varargin{j};
    elseif strncmpi(varargin{j},'smooth',2)
        j = j+1;    
        if varargin{j} > 0
            sd = varargin{j} * 4 * range(dat)/sqrt(length(dat));
        end
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        if length(varargin) > j
            plotargs = varargin{j:end};
        else
            plotargs = {};
        end
    elseif strncmpi(varargin{j},'xval',2)
        j = j+1;    
        setx =  varargin{j};
    end
    j = j+1;
end
if step > 2 * sd 
    step = sd;
end
aw = 2 * sd.^2;

if sd == 0
    smoothed = dat;
    x = min(dat):range(dat)/(length(dat)+1):max(dat)
    return
end

if isempty(dat)
    return;
end
if isempty(drange)
    if period
        drange(1) = (min(dat)+max(dat) - period)/2;
        drange(2) = drange(1) + period;
    else
        drange(1) = min(dat)-sd;
        drange(2) = max(dat) + sd;
    end
end
if length(drange) > 2
    x = drange;
else
    x = drange(1):step:drange(2)+0.8 * step;
end

if ~isempty(setx)
    x = setx;
end

boxcar = 1;
if boxcar == 1
    j = 1;
    for j = 1:length(x)
        smoothed(j) = sum(exp(-(dat-x(j)).^2/aw));
    end
elseif boxcar == 2
    [b,a] = hist(dat);

k = trapz(a,b)/trapz(x,smoothed);
smoothed = smoothed * k;

else
    j = 1;
    for j = 1:length(x)
        smoothed(j) = sum(abs(dat-x(j))< sd);
    end
end

if exist('plotargs','var')
    plot(x,smoothed,plotargs{:});
end
