function rogs = RoG(k,sd, varargin)
if nargin < 2
    sd = 10;
end
if nargin < 1
    k = 10;
end
sds = [10:100 200:100:10000];


x = -999:1000;
%center Gaussian
Gc = exp(-x.^2./(2.*sd.^2));
Gc = Gc./sum(Gc);
%surround
Gs = exp(-x.^2./(3.*2.*sd.^2));
Gs = Gs./sum(Gs);

periods = [2:2:200];

for j = 1:length(periods)

    %make square wave
    period = periods(j);
    y = round(mod(x,period)./period)-0.5;
    %convolve with Gaussians
    zc = conv(y,Gc);
    zs = conv(y,Gs);
    %just take central period
    zc = zc(2500-period:2500+period);
    zs = zs(2500-period:2500+period);
    %do RoG
    ratio = abs(zc)./(1+(k.*abs(zs)));
    %look at ratio
    ratios{j} = ratio;
    mr(j) = max(zc)./max(zs);
    %record ROG resp, and each Gaussian
    rogs(j) = mean(abs(ratio));
    resp(j,1) = mean(zc(zc>0));
    resp(j,2) = mean(zs(zs>0));
end
plot(periods,resp(:,1));
hold on;
plot(periods,resp(:,2),'r');
plot(periods,rogs .* max(resp(:))./max(rogs),'k');
%what if do RoG after summing Gaussian convolutions
rog = resp(:,1)./(1+(k.*resp(:,2)));
plot(periods,rog .* max(resp(:))./max(rog),'k--');
