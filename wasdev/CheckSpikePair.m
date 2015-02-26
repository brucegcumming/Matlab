function [indices, rawac, smoothed] = CheckSpikePair(E1, E2, varargin)


adcmode = 1;
dispw = 300;
j = 1;
indices = [];
rawac = [];
smoothed = [];
showplot = 1;




latency = 500;

showplot = 1;
j = 1;
while(j < nargin-1)
    if(strncmpi(varargin{j},'mean',4))
        type = 'mean';
    elseif(strncmpi(varargin{j},'noplot',4))
        showplot = 0;
    elseif(strncmpi(varargin{j},'adcc',4))
        adcmode = 2;
    elseif(strncmpi(varargin{j},'dispw',5))
        j = j+1;
        dispw = varargin{j};
    end
    j = j+1;
end

if(length(E1.Trials) == 0 | length(E2.Trials) == 0)
    if showplot
        delete(allchild(gca));
    end
    return;
end

[c, ai, bi] = intersect([E1.Trials.Start],[E2.Trials.Start]); % get matching trials
isis = [];
aisis = [];
bisis = [];
for j = 1:length(ai)
     Expt.Trials(j) = E1.Trials(ai(j));
     Expt.Trials(j).Spikes = sort([E1.Trials(ai(j)).Spikes; E2.Trials(ai(j)).Spikes]);
     isis = [isis diff(Expt.Trials(j).Spikes')];
     aisis = [aisis diff(E1.Trials(ai(j)).Spikes')];
     bisis = [bisis diff(E2.Trials(bi(j)).Spikes')];
end

duration = min([Expt.Trials.End] - [Expt.Trials.Start]);

rawac = autocorrelate(Expt.Trials,latency, duration);
rawac(1) = 0;
sd = 4;
w = 12;
kernel = gauss(sd,-w:w);
kernel(1:w) = 0; %% Half Gaussian so that smearing is only forward in time
kernel = kernel ./sum(kernel);
smoothed = conv(rawac(1:dispw),kernel);
smoothac{1} = smoothed(w:end-w);


subplot('position',[0.05 0.5 0.45 0.45]);
hold off;
plot(smoothac{1});
hold on;
rawac = autocorrelate(E1.Trials,latency, duration);
rawac(1) = 0;
smoothed = conv(rawac(1:dispw),kernel);
smoothac{2} = smoothed(w:end-w);
plot(smoothac{2},'r');
rawac = autocorrelate(E2.Trials,latency, duration);
rawac(1) = 0;
smoothed = conv(rawac(1:dispw),kernel);
smoothac{3} = smoothed(w:end-w);
plot(smoothac{3},'g');

xc = SpikeXcorr(E1.Trials,E2.Trials);
smoothed = conv(xc(1:dispw),kernel);
smoothac{4} = smoothed(w:end-w);
scale = sum(smoothac{1})/sum(smoothac{4});
plot(smoothac{4}.* scale,'k');


subplot('position',[0.05 0.05 0.45 0.45]);
hold off;
x = 0:10000;
isisum = conv(hist(isis,x),kernel);
isi{1} = isisum(w:end-w);
plot(isi{1}(1:dispw));
isisum = conv(hist(aisis,x),kernel);
isi{2} = isisum(w:end-w);
hold on;
plot(isi{2}(1:dispw),'r');
isisum = conv(hist(bisis,x),kernel);
isi{3} = isisum(w:end-w);
plot(isi{3}(1:dispw),'g');
plot(isi{3}(1:dispw)+isi{2}(1:dispw),':');

subplot('position',[0.5 0.05 0.45 0.9]);
hold off;
if adcmode == 1
    plot(E1.Spike.Mean(:,1),'r');
    hold on;
    plot(E2.Spike.Mean(:,1),'g');
else
    plot(E1.Spike.Mean(:,1),E2.Spike.Mean(:,1));
end
title(splitpath(E1.Header.Name));
