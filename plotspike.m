function dp = plotspike(SPK, varargin)

type = 'roc';
showplot = 1;
j = 1;
while(j < nargin)
    if(strncmpi(varargin{j},'mean',4))
        type = 'mean';
    elseif(strncmpi(varargin{j},'noplot',4))
        showplot = 0;
    end
    j = j+1;
end

%SPK.DDF is a frequency histogram for the radial distances for
%spikes from the centre of the cluster.
%It has two colums first is the distance for all events. The second
%is the distances for just those events that were classified as
%belonging to the cluster

if showplot
    subplot(1,2,1);
    hold off;
    plot(SPK.DDF(:,1),'g');
    hold on;
    plot(SPK.DDF(:,2));
    plot(SPK.DDF(:,1)- SPK.DDF(:,2),'r');
end

sumc = sum(SPK.DDF(:,2));
suma = sum(SPK.DDF(:,1) - SPK.DDF(:,2));

%detection rate
if sumc
    drate = cumsum(SPK.DDF(:,2)) ./ sumc;
%false positive rate
fpos = cumsum(SPK.DDF(:,1) - SPK.DDF(:,2)) ./ suma;
%area gives ROC
dp = trapz(fpos,drate);
else  %% no events in cluster 0 = must be good
    drate = 1;
    fpos = 0;
    dp = 1;
end
if ~showplot
    return;
end

title(sprintf('Area under ROC %.5f',dp));


subplot(1,2,2);
hold off;
if strmatch(type,'roc')
  plot(fpos,drate);
else
  plot(SPK.Mean(:,1));
  hold on;
  plot(SPK.Mean(:,2),'r');
end

