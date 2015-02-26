function PlotExptSeq(Expt, varargin)
%Plot sequnece of events in an Expt. Default is to plot Stimon/off times

spktimes = [];
FullV = [];
j = 1;
plottype = 'stimon';

while j <= length(varargin)
    if isfield(varargin{j},'blkstart')
        FullV = varargin{j};
    elseif strncmpi(varargin{j},'Spikes',5)
        j = j+1;
        spktimes = varargin{j};
    elseif strncmpi(varargin{j},'ids',3)
        plottype = 'ids';
    elseif strncmpi(varargin{j},'reward',3)
        plottype = 'reward';
    end
    j = j+1;
end
T = Expt.Trials;
if strcmp(plottype,'stimon')
    for j = 1:length(Expt.Trials)
        start(j) = Expt.Trials(j).Start(1);
        ends(j) = Expt.Trials(j).End(end);
        plot([start(j) start(j) ends(j) ends(j)],[0 1 1 0]);
        hold on;
    end
    if isfield(FullV,'blkstart')
        for j = 1:length(FullV.blkstart)
            plot([FullV.blkstart(j) FullV.blkend(j)].*10000,[1.1 1.1],'r');
        end
    end
        
elseif strcmp(plottype,'ids')
    if isfield(T,'date')
        plot([T.date] ,[T.id]);
        datetick;        
    else
        plot(T.Start,rw);        
    end
elseif strcmp(plottype,'reward')
    if isfield(Expt.Trials,'Result')
        id = find([Expt.Trials.Result] > 0);
        T = Expt.Trials(id);
    end
    rw = cumsum([T.rw]);
    if isfield(T,'date')
        plot([T.date] ,rw);
        datetick;        
    else
        plot(T.Start,rw);        
    end
end

if isfield(Expt,'sessions')
    y = get(gca,'ylim');
    for j = 1:length(Expt.sessions)
        x = [Expt.sessions(j) Expt.sessions(j)];
        line(x,y);
    end
end
if ~isempty(spktimes)
    plot(spktimes, 1.1,'k.');
end