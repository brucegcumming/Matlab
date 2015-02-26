function fitp = DistanceVoltage(varargin)
%Simulate effect of quantized sampling in laminar arrays on
%distance-voltage relationship for 50uM spaced samples



x = -1000:1000; %in microns
lambda = 20; %true space constant in microns
d = exp(-abs(x)/lambda);
da = exp(-abs(x)/(lambda * 1.2));
for j = 1:50
    alldx{j} = [];
end
for j = 1:50
    di = (1+j):50:length(x);
    amps = d(di);
    [a,b] = max(amps);
    dxs = abs([1:length(amps)] -b);
    maxamp(j) = a;
    maxpt(j) = b;
    amps = amps./a;
    for k = 1:length(dxs)
        alldx{dxs(k)+1}(end+1) = amps(k);
    end
end


hold off
for j = 1:11
    plot(j-1, alldx{j},'o');
    means(j) = nanmean(alldx{j});
    distances(j) = (j-1)*50;
    hold on;
end
plot(0:10,means(1:11),'k-');
plot([0:500]/50,d(1001:1501),'r-');
%plot([0:500]/50,da(1001:1501),'g-');

p = 1.2;
fitp = fminsearch(@matchexp, p, [],  means, distances, lambda);
da = exp(-(abs(x/lambda))/(fitp(1)));
plot([0:500]/50,da(1001:1501),'b-');


p(1) = 1.2;  %space constant
p(2) = 1; %exponent

fitp = fminsearch(@matchexp, p, [],  means, distances, lambda);
da = exp(-(abs(x/lambda).^fitp(2))/(fitp(1)));
plot([0:500]/50,da(1001:1501),'g-');


function rss = matchexp(x,y, d, lambda)

if length(x) ==2
fity = exp(-((d/lambda).^x(2))/(x(1)));
else
fity = exp(-(d)/(x(1) * lambda));
end
rss = sum((y-fity).^2);