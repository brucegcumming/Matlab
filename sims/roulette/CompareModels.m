function [e,r,p] = CompareModels(varargin)

showplot = 1;
fits{1}.params = [0.3101 4.3749 1.3307 -0.0166 0.9084 1.4620];
fits{1}.type = 1;
smoothw = 20;
periods = [];

reload = 0;

j = 1;
while j <= nargin
    if strncmpi(varargin{j},'reload',3)
        reload = 1;
    elseif strncmpi(varargin{j},'smoothw',3)
        j = j+1;
        smoothw = varargin{j};
    elseif strncmpi(varargin{j},'period',3)
        j = j+1;
        periods = varargin{j};
    end
    j = j+1;
end

if showplot
    colors = mycolors;
    GetFigure('ShowDist');
    subplot(1,1,1);
    hold off;
    GetFigure('ShowEdge');
    subplot(1,1,1);
    hold off;
    GetFigure('ShowScale');
    subplot(1,1,1);
    hold off;
end


nmodels = 5;
labels = {'Gauss','2Gauss','+period','2+Skew','skew'};
if reload
    a = textread('new.txt');  %unwrapped, smoothed data
    fdist = a(:,2);
    sd = std(fdist);
    % robs units are 370 to the circle
    xvals = [1:length(fdist)] * 2 * pi/370;
    for nm = 1:nmodels
        fits{nm} = FitDist(xvals,fdist',nm-1);
%    angles = f2dist(fdist,1000,'xvals',[1:length(fdist)]);
    end
    save('RobFits.mat','fits','xvals','fdist');
else
    load('RobFits.mat');
end

if ~isempty(periods)
    stds = [];
    for j = 1:length(periods)
        px = mod(xvals,periods(j));
        [a,b] = sort(px);
        sm = smooth(fdist(b),50,'gauss');
        plot(sm,'color',colors{j});
        hold on;
        stds(j) = range(sm);
    end
    GetFigure('ShowEdge');
    plot(periods,stds);
    return;
end
xv = -8*pi:0.01:8*pi;
scales = [0.5:0.1:4];
for nm = 1:nmodels+2
    if nm == nmodels+1
        y = fdist';
        [m,b] = max(fdist);
        xv = xvals-xvals(b);
        labels = {labels{:} 'raw data'}
    elseif nm == nmodels+2
        y = smooth(fdist',smoothw,'gauss');
        [m,b] = max(fdist);
        xv = xvals-xvals(b);
        labels = {labels{:} 'smoothed data'}
    else
        y = FitDist(xv,fits{nm}.params,fits{nm}.type,'eval');
    end
    if showplot
        GetFigure('ShowDist');
        plot(xv,y,'color',colors{nm});
        hold on;
    end
    [e, r, p] = RunDists(scales,'dist',[y; xv]);
    
    if showplot
        GetFigure('ShowEdge');
        plot(r,e,'color',colors{nm});
        hold on;
        GetFigure('ShowScale');
        plot(scales,e,'color',colors{nm});
        hold on;
    end
end


j = nm+1;
y = fdist';
[e,r,p] = RunDists([0.5:0.1:4],'dist',[y; xvals-xvals(b)]);

legend(labels);