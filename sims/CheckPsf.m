function [fits, details] = CheckPSF(nloops, varargin)
showplot = 0;
details = [];
nresample = 1;
ntrials = 120;
nbins = 6;
x = [-4:1.6:4];
j = 1;


while j <= length(varargin)
    if strncmpi(varargin{j},'nbins',4)
        j = j+1;
        nbins = varargin{j};
    elseif strncmpi(varargin{j},'ntrials',4)
        j = j+1;
        ntrials = varargin{j};
    elseif strncmpi(varargin{j},'nresample',4)
        j = j+1;
        nresample = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'xvals',4)
        j = j+1;
        x = varargin{j};
    end
    j = j+1;
end

if isstruct(nloops)
    PlotPSFSim(nloops);
    return;
end
y = cdf('norm',x,0,1);
hold off;
if nbins == 0
    xr = range(x);
    xo = max(x);
else
plot(x,y);
end
for k = 1:nloops
if nbins == 0
rnd = rand(ntrials,1);
x = (rand(ntrials,1) .* xr) - xo;
y = cdf('norm',x,0,1);
for j = 1:length(x)
    pp(j).x = x(j);
    pp(j).n = size(rnd,2);
    pp(j).resp = rnd(j) < y(j);
end
else
rnd = rand(length(x),ntrials./length(x));
for j = 1:length(x)
    pp(j).x = x(j);
    pp(j).n = size(rnd,2);
    pp(j).resp = sum(rnd(j,:) < y(j));
end
end
resps(k,:) = [pp.resp]./[pp.n];
if showplot
fit = fitpsf(pp,'showfit');
drawnow;
else
fit = fitpsf(pp);
fity = cdf('norm',x,fit.fit(1),fit.fit(2));

if nbins == 0 
    rnd = rand(length(x),nresample);
    for n = 1:nresample
    for j = 1:length(x)
    pp(j).resp = rnd(j,n) > fity(j);
    end
    nfit = fitpsf(pp);
    rfits(n,:) = nfit.fit;
    end
else
    for n = 1:nresample
    rnd = rand(length(x),mean([pp.n]));
    for j = 1:length(x)
    pp(j).resp = sum(rnd(j,:) > fity(j));
    end
    nfit = fitpsf(pp);
    rfits(n,:) = nfit.fit;
    end
end

end
fits(k,:) = fit.fit;
if nresample
    details.rfits{k} = rfits;
end
end
hold on;
plot(x,mean(resps),'o');
details.psd = fits(:,2);
function PlotPSFSim(res)

for j = 1:length(res.rfits)
    psds(1:2,j) = prctile(abs(res.rfits{j}(:,2)),[2.5 97.5]);
end
plot(res.psd,psds(1,:),'o');
hold on;
plot(res.psd,psds(2,:),'ro');
plot(get(gca,'xlim'),[1 1],':');
np = sum(psds(1,:)> 1)./10;
pp = sum(psds(2,:)< 1)./10;
text(1.4,1.2,sprintf('%.1f%%',np));
text(0.4,0.9,sprintf('%.1f%%',pp));
xlabel('SD fitted to sample');
ylabel('Confidence interval');