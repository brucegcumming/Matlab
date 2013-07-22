function [dps, details] = GaussNMF(sd, varargin)
%[dps, details] = GaussNMF(sd, varargin)
% calculate dprime for Gaussian tuning functions assuming poisson
% statisitcs (var = k.mean)
%
%GaussNMF([30 30],'gains',[55 55 * 0.64],'offsets',[45 12.3])
%approximates data of straight line fit in Figure S6 of Lee et al 2012

gains = [];
offsets = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'gains',4)
        j = j+1;
        gains = varargin{j};
    elseif strncmpi(varargin{j},'offsets',4)
        j = j+1;
        offsets = varargin{j};
    end
    j = j+1;
end
        
if isempty(offsets)
    offsets = zeros(size(sd));
end

if isempty(gains)
    gains = ones(size(sd));
end

for j = 1:length(sd)
x = -90:90;
y = Gauss([0 sd(j) gains(j) offsets(j)],x);
y(y <0) = 0;
details.fits(j,:) = y;

cv = sqrt((y(1:end-1)+y(2:end))/2);
dp = diff(y)./cv;
yx = y(2:end) - mean(diff(y))/2;
dps(j) = max(dp(yx>gains(j)/4));
end
if length(sd) == 1
    plot(dp,'o-');
else
    plot(sd,dps,'o-');
end