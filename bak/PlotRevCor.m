function result = PlotRevCor(Expt, varargin)
%
% Add uncorr to plot
% make plots for AC -> and COR ->
%get default from variance wrt time.
times = -100:2000;
np = 1;
sdfw = 20;
plotsum = 0;
plotbest = 0;
slices = [500 600 700 800 900 1000];
xrange = [];
labela = 'ACRCFiga';
labelb = 'ACRCFigb';
secondorder = 0;
nloops = 1;
legendpos = 0;
result = [];
showplot = 1;
needvar = 0;
autoslice = 1;
iskip = 1;

j = 1; 
while j < nargin
    if strncmpi(varargin{j},'sdfw',4)
        j = j+1;
        sdfw = varargin{j};
    elseif strncmpi(varargin{j},'autoslices',4)
        autoslice = 1;
    elseif strncmpi(varargin{j},'slices',4)
        j = j+1;
        slices = varargin{j};
    elseif strncmpi(varargin{j},'add',3)
        plotsum = 1;
    elseif strncmpi(varargin{j},'range',3)
        j = j+1;
        xrange = varargin{j};
    elseif strncmpi(varargin{j},'best',4)
        plotbest = 1;
    elseif strncmpi(varargin{j},'figa',4)
        j = j+1;
        labela = varargin{j};
    elseif strncmpi(varargin{j},'figb',4)
        j = j+1;
        labelb = varargin{j}; 
    elseif strncmpi(varargin{j},'noplot',5)
        showplot = 0;
    elseif strmatch(varargin{j},'legendpos')
      legendpos = varargin{j+1};
      j = j+1;
    elseif strncmpi(varargin{j},'skip',4)
      iskip = varargin{j+1};
      j = j+1;
    elseif strncmpi(varargin{j},'splitce',6)
        secondorder = 1;
        nloops = 2;
    elseif strncmpi(varargin{j},'splitdx',6)
        secondorder = 2;
        nloops = 2;
    elseif strncmpi(varargin{j},'showvar',6)
        needvar = 1;
    end
    j = j+1;
end

Expt.Header.Name = strrep(Expt.Header.Name,'\','/');
if(autoslice)
    needvar = 1;
end
colors = mycolors;
starts = [Expt.Trials.Start];

xvals = unique([Expt.Trials.dO]);
yvals = unique([Expt.Trials.ce]);

if ~isempty(xrange)
    idx = find(xvals >= xrange(1) & xvals <= xrange(2));
    xvals = xvals(idx);
end
frameperiod = 10000/Expt.Stimvals.fz;
endoffset = 5; %% 0.5ms
for j = 1:length(Expt.Trials)
  tx(:,j) = Expt.Trials(j).dO;
  ty(:,j) = Expt.Trials(j).ce;
  durs(:,j) = [diff(Expt.Trials(j).Start); endoffset + Expt.Trials(j).End(end) - Expt.Trials(j).Start(end)];
end

stimdur = mode(durs);
result.nbad = length(find(abs(durs(2:end,:) - stimdur) > frameperiod/3));
badtimes = find(abs(durs - stimdur) > frameperiod/3);
good(badtimes) = 0;
good(1,:) = 1;
good = ones(size(durs));
if(iskip)
    good(1:iskip,:) = 0;
end
GetFigure(labela);
hold off;
h = [];

if size(xvals,1) > size(xvals,2)
    xvals = xvals';
end


for loop = 1: nloops;
for x = xvals;
    for y = [1 -1 0];
        tidx = [];
        for j = 1:length(Expt.Trials)
            if secondorder == 1
                if loop == 1
                    fidx = find(ty(:,j) == 1);
                else
                    fidx = find(ty(:,j) == -1);
                end
                idx = find(tx(:,j) == x & ty(:,j) == y);
                idx = intersect(idx,fidx+1);
            elseif secondorder == 2
                if loop == 1
                    fidx = find(tx(:,j)  < 0 & ty(:,j) ==1);
                else
                    fidx = find(tx(:,j) > 0 & ty(:,j) ==1);
                end
                idx = find(tx(:,j) == x & ty(:,j) == y);
                idx = intersect(idx,fidx+1);
            else
                idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j));
            end
            if ~isempty(idx)
                Expt.Trials(j).Trigger = Expt.Trials(j).Start(idx)' - Expt.Trials(j).TrialStart;
                tidx = [tidx j];
            end
        end
        if ~isempty(tidx)
            [sdf, n] = trigsdf(Expt.Trials(tidx),sdfw,times,'halfgauss');
            sdfs(np).x = x;
            if y > 0
                sdfs(np).corr = sdf;
                if( ~plotsum)
                    h(np) = plot(times/10,sdf,'color',colors{np});
                    hold on;
                end
            elseif y < 0
                sdfs(np).ac = sdf;
                if ~plotsum
                    plot(times/10,sdf,':','color',colors{np});
                    hold on;
                end
            elseif y == 0 & ~isempty(tidx)
                sdfs(1).uc = sdf;
                if ~plotsum
                    plot(times/10,sdf,'--','color','k');
                    hold on;
                end
            end
            if(plotsum)
                h(np) = plot(times/10,sdfs(np).corr +sdfs(np).ac,'color',colors{np});
                hold on;
            end        
            labels{np} = sprintf('%.2f n = %d',x,n);
        end
    end
    if ~isempty(sdfs(np).corr)
        np = np+1;
    end
end
end
if(legendpos < 6)
    legend(h,labels,legendpos);
end

if ~isfield(sdfs,'uc')
    sdfs(1).uc = [];
end
result.sdfs = sdfs;
GetFigure(labelb);
hold off;

h = [];
labels = {};
np = 1;

bestvar = 0;

if plotbest | needvar
    nv = 1;
    for j = 400:20:length(sdfs(1).corr)
        x = [];
        y = [];
        z = [];
        for k = 1:length(sdfs)
            x = [x sdfs(k).x];
            y = [y sdfs(k).corr(j)];
            z = [z sdfs(k).ac(j)];
        end
        dvar = var(y) + var(z);
        vars(nv) = dvar;
        delays(nv) = times(j);
        result.x(:,nv) = x;
        result.y(:,nv) = y;
        result.z(:,nv) = z;
        if  ~isempty(sdfs(1).uc)
            result.uc(nv) = sdfs(1).uc(j);
        else
            result.uc = [];
        end
        if(dvar > bestvar)
            bestvar = dvar;
            bestj = j;
            bestx = x;
            besty = y;
            bestz=  z;
            bestnv = nv;
        end
        nv = nv + 1;
    end
    result.vars = vars;
    result.delays = delays;
    result.bestdelay = bestnv;
    result.timeoff = times(1);
    if showplot & plotbest
        result.h(1) = plot(bestx,besty,'-','color',colors{np});
        hold on;
        result.h(2) = plot(bestx,bestz,':','color',colors{np});
        title(sprintf('%s at %d',Expt.Header.Name,bestj));
    end
    GetFigure(labela);
    plot(delays./10,sqrt(vars).*4,'k','linewidth',2);
end

if ~plotbest
    GetFigure(labelb);
if autoslice
    [mv, mj] = max(result.vars);
    minv = min(result.vars);
    th = minv + (mv-minv)/10;
    id = find(result.vars > th);
    step = range(result.delays(id))/6;
    slices = round(result.delays(id(1)):step:result.delays(id(end)));
end
    
for j = slices;
    x = [];
    y = [];
    z = [];
    for k = 1:length(sdfs)
        x = [x sdfs(k).x];
        y = [y sdfs(k).corr(j)];
        z = [z sdfs(k).ac(j)];
        u = sdfs(1).uc(j);
    end
       h(np) = plot(x,y,'color',colors{np});
       hold on;
       plot(x,z,':','color',colors{np});
       plot([min(x) max(x)],[u u],'--','color',colors{np});
       labels{np} = sprintf('dT = %.0fms',j/10);
       np = np+1;
end
if legendpos < 6
    legend(h,labels);
end
title(sprintf('%s',Expt.Header.Name));
end
