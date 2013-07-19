function PlotEm(Trials, varargin)

plottype = 0;
plotavg = 0;
figlabel{1} = 'EmPlotYT';
j = 1;
while j < nargin
    if strncmpi(varargin{j},'cv',2)
        plottype = 1;
    elseif strncmpi(varargin{j},'avg',2)
        plotavg = 1;
    elseif strncmpi(varargin{j},'xy',2)
        plottype = 2;
    end
    j = j + 1;
end


if plotavg
    buflen = NaN;
    for j = 1:length(Trials);
        buflen = min([buflen length(Trials(j).Evyevals.lh]);
    end
    lvs = [];
    rvs = [];
    lhs = [];
    rhs = [];
    for trial = 1:length(Expt.Trials)
        lvs(trial,:) = Expt.Trials(trial).Eyevals.lv(1:buflen)';
        rvs(trial,:) = Expt.Trials(trial).Eyevals.rv(1:buflen)';
        rhs(trial,:) = Expt.Trials(trial).Eyevals.rh(1:buflen)';
        lhs(trial,:) = Expt.Trials(trial).Eyevals.lh(1:buflen)';
    end
    ev.lh = mean(lhs);
    ev.rh = mean(rhs);
    ev.lv = mean(lvs);
    ev.rv = mean(rvs);
    Trials = [];
    Trials(1).Eyevals = ev;
    Trials(1).Trial = 0;
end

GetFigure(figlabel{1});

hold off;
for j = 1:length(Trials);
    ev = Trials(j).Eyevals;
    if plottype == 1
        plot(Trials(j).Eyevals.lh-Trials(j).Eyevals.rh,'r');
        hold on;
        plot(mean([Trials(j).Eyevals.lh Trials(j).Eyevals.rh],2),'g');
        plot(mean([Trials(j).Eyevals.lv Trials(j).Eyevals.rv],2),'b');
    elseif plottype == 2
        for k = 1:length(ev.lh);
            plot(ev.lh(1:k),ev.lv(1:k),'r');
            hold on;
            plot(ev.rh(1:k),ev.rv(1:k),'r');
            if mod(k,100) == 99
                drawnow;
            end
        end
    else
        plot(Trials(j).Eyevals.lh,'r');
        hold on;
        plot(Trials(j).Eyevals.rh,'g');
        plot(Trials(j).Eyevals.lv,'m');
        plot(Trials(j).Eyevals.rv,'c');
    end
end