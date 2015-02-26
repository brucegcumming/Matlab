function dp = SpikeShape(SPK, varargin)
%  dp = SpikeShape(SPK, varargin)
% Calculates Spike width, and isolation index
% width is difference between first two zero crossings, exluding any
% early ones that occur at very low dV/dt
%
% returns structure dp:
% dp.width peak-trough distance (#samples)
% dp.iwidth interopolated version i.e. finer resolution.
% dp.dprime  ROC measure of overlap between SU/MU parameter distribtions
% dp.extrema list of sample #s associated with extrema (may be > 2)
% SpikeShape(Expt.Spike, 'plot')  plots the spike shape, its temporal
% derivative, and marks extrema

type = 'roc';
thscale = 5;  %dV/dT > max(dV/dt)/thscale required before peak
showplot = 0;
j = 1;
while(j < nargin)
    if(strncmpi(varargin{j},'mean',4))
        type = 'mean';
    elseif(strncmpi(varargin{j},'plot',4))
        showplot = 1;
    elseif(strncmpi(varargin{j},'th',4))
        j = j+1;
        thscale = varargin{j};
    end
    j = j+1;
end

%SPK.DDF is a frequency histogram for the radial distances for
%spikes from the centre of the cluster.
%It has two colums first is the distance for all events. The second
%is the distances for just those events that were classified as
%belonging to the cluster


sumc = sum(SPK.DDF(:,2));
suma = sum(SPK.DDF(:,1) - SPK.DDF(:,2));

%detection rate
if sumc
    drate = cumsum(SPK.DDF(:,2)) ./ sumc;
    %false positive rate
    fpos = cumsum(SPK.DDF(:,1) - SPK.DDF(:,2)) ./ suma;
    %area gives ROC
    dp.dprime = trapz(fpos,drate);

    if size(SPK.DDF,2) > 2
        sumd = sum(SPK.DDF(:,4));
        suma = sum(SPK.DDF(:,3) - SPK.DDF(:,4));
        drate = cumsum(SPK.DDF(:,4)) ./ sumd;
    %false positive rate
    fpos = cumsum(SPK.DDF(:,3) - SPK.DDF(:,4)) ./ suma;
    dp.hdprime = 1 - trapz(fpos,drate);
    end
else
    dp.dprime = NaN;
end
dvdt = diff(SPK.Mean);
th = max(abs(dvdt));
minrate = th(1)/thscale;
dp.extrema = [];
dp.width = [];
if th(1) > 0
    while length(dp.extrema) < 2 && isempty(dp.width)
        start = find(abs(dvdt(:,1)) > minrate);
        zc = diff(sign(dvdt(:,1)));
        dp.extrema = find(abs(zc(start(1):end)) > 0) + start(1);
        if length(dp.extrema) > 1
            dp.width = diff(dp.extrema(1:2));
        end
        if minrate < th(1)/100 %% really is just one extrema, eg ruf673.0.SFM
            dp.width = NaN;
        end
        minrate = minrate/2;
    end
else
    dp.width = NaN;
end

%Also calculate interopolated extrema
spkmean = SPK.Mean(:,1);
for j = 1:length(dp.extrema)
    pp = dp.extrema(j);
    vdiffs = [abs(spkmean(pp+1)-spkmean(pp)) abs(spkmean(pp)-spkmean(pp-1))];
    dp.ix(j) = pp + (vdiffs(2) - vdiffs(1))/sum(vdiffs);
end

if length(dp.extrema) > 1
    dp.iw = diff(dp.ix(1:2));
else
    dp.iw = NaN;
end

[a, b] = max(spkmean);
if b == 1 | b == 32
    dp.peak = 1;
else
    vdiffs = [abs(spkmean(b+1)-spkmean(b)) abs(spkmean(b)-spkmean(b-1))];
    dp.peak = b + (vdiffs(2) - vdiffs(1))/sum(vdiffs);
end

[a, b] = min(spkmean);
if b == 1 | b == 32
    dp.valley = b;
else
    vdiffs = [abs(spkmean(b+1)-spkmean(b)) abs(spkmean(b)-spkmean(b-1))];
    dp.valley = b + (vdiffs(2) - vdiffs(1))/sum(vdiffs);
end
dp.vw = abs(dp.peak-dp.valley);
dp.h = range(SPK.Mean);
if showplot
    hold off;
    plot(SPK.Mean(:,1));
    spkrange = [min(SPK.Mean(:,1)) max(SPK.Mean(:,1))];
    hold on;
    for j = 1:length(dp.extrema)
        plot([dp.extrema(j) dp.extrema(j)],spkrange,':');
    end
    plot([dp.peak dp.peak],spkrange,'r:');
    plot([dp.valley dp.valley],spkrange,'r:');
    plot(dvdt(:,1) .* diff(spkrange)/range(dvdt(:,1)),'r');
%    plot(zc .* mean(abs(spkrange)),'g');
end


