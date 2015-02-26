function [indices, rawac, smoothed] = CheckSpike(Expt, varargin)

showplot = 1;
if ischar(Expt)
    load(Expt);
end

duration = min([Expt.Trials.End] - [Expt.Trials.Start]);
latency = 500;

showplot = 1;
showisi = 0;

j = 1;
while(j < nargin)
    if(strncmpi(varargin{j},'mean',4))
        type = 'mean';
    elseif(strncmpi(varargin{j},'noplot',4))
        showplot = 0;
    elseif(strncmpi(varargin{j},'isi',3))
        showisi = 1;
    end
    j = j+1;
end

rawac = autocorrelate(Expt.Trials,latency, duration);
rawac(1) = 0;
sd = 4;
w = 12;
kernel = Gauss(sd,-w:w);
kernel(1:w) = 0; %% Half Gaussian so that smearing is only forward in time
kernel = kernel ./sum(kernel);
smoothed = conv(rawac(1:200),kernel);
smoothed = smoothed(w:end-w);

indices(1) = mean(rawac(1:15))/mean(rawac(100:200));
indices(2) = mean(rawac(1:15))/max(smoothed);
indices(3) = PlotSpike(Expt.Spike,'noplot');
if showisi
    isis = [];
    times = [];
    for j = 1:length(Expt.Trials)
        isis = [isis diff(Expt.Trials(j).Spikes)'];
        times = [times Expt.Trials(j).Start(1)+Expt.Trials(j).Spikes(2:end)'];
    end           
end
if(showplot)
    indices(3) = PlotSpike(Expt.Spike);
    plot(smoothed);
      uicontrol(gcf,'Style', 'pushbutton', 'Callback', ['ReplaySpikes('' ' Expt.Header.Name ' '');'],...
'String', 'Replay', 'Position', [10 10 40 20]);
else
    indices(3) = PlotSpike(Expt.Spike,'noplot');
end
