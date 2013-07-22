function result = simplerc(varargin)

% simulate timcourse of simple reverse correlation sequence.
%
%
% simplerc('expon',3,'plotloop'); vs simplerc('expon',3,'plotloop','noise',10); 
%    shows that variation in addition to the stim reducess the effect of OPNL
%simplerc('threshold',5,'plotloop');
%    shows that a simpel threshold looks very like an exponent =
%    continuous bending. Also flattened by noise just like the exponend
%    is.
% simplerc('threshold',5,'log'); plots resp as log(respi/resp0), to test
% the linear system eqn from Ringach. 
colors = mycolors;
sd = 2;
thr = NaN;
plotresp = 0;
plotloop= 0;
plotlog = 0;
plotpeak = 0;
nstim = 100000;
sigma = 0;
gamma = 1;
framerepeat = 1;
result = [];
done = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'thresh',4)
        j = j+1;
        thr = varargin{j};
    elseif strncmpi(varargin{j},'expons',4)
        j = j+1;
        gamma = varargin{j};
    elseif strncmpi(varargin{j},'nstim',5)
        j = j+1;
        nstim = varargin{j};
    elseif strncmpi(varargin{j},'noise',5)
        j = j+1;
        sigma = varargin{j};
    elseif strncmpi(varargin{j},'figure',3)
        j = j+1;
        MakeFigure(varargin{j});
        done = done+1;
    elseif strncmpi(varargin{j},'frpt',4)
        j = j+1;
        framerepeat = varargin{j};
    elseif strncmpi(varargin{j},'logresp',3)
        plotlog = 1;
     elseif strncmpi(varargin{j},'plotloop',5)
        plotloop = 1;
    elseif strncmpi(varargin{j},'tsd',3)
        j = j+1;
        sd = varargin{j};
    elseif strncmpi(varargin{j},'plotr',5)
        plotresp = 1;
    elseif strncmpi(varargin{j},'plotpeak',5)
        plotpeak = 1;
    end
    j = j+1;
end
if done
    return;
end
xvals = [0:10];
x = floor(rand(nstim,1) * length(xvals)); %random integers 0-10
if framerepeat > 1
    for j = 1:framerepeat
        y(j:framerepeat:nstim*framerepeat) = x;
    end
    x = y;
end

scale = 1./(mean(xvals).^gamma);
noise = randn(size(x)) * sigma;
k = Gauss(sd,-10:10);
sum(k);
resp = conv(x+noise,k);
if ~isnan(thr)
    resp = resp - thr;
    resp(resp<0) = 0;
end
mresp = mean(resp);
resp =  scale .* resp.^gamma;
y = resp(10:end-10);
%y is Gaussian convoled with 1-D random time series

np = 1;
%Now find temporal response to three stimuli, treating this
%as a reverse correlation sequence.
k = 0;
if plotpeak
    stims = xvals;
else
    stims = [min(x) median(x) max(x)];
end

for stim = stims
    k = k+1;
    starts = 1:framerepeat:length(x);
idx = find(x(starts) == stim);
all = zeros(21,length(idx));
for j = 1:length(idx)
    id = starts(idx(j));
    all(:,j)= resp(id:id+20);
end
   allresp(k,:) = mean(all,2);
end
if plotresp
    plot(y);
elseif plotloop
    b = mean(y);
    hold off;
    plot(allresp(1,:),allresp(end,:));
    m = b-min(allresp(1,:));
    hold on;
    plot([b b-m],[b b+m],'k');
elseif plotlog
    for j = 2:size(allresp,1)
    plot(log(allresp(j,:)./allresp(1,:)),'color',colors{j});
    hold on;
    end
    plot([0 size(all,1)],[mean(y) mean(y)],'color',colors{j+1});
elseif plotpeak
    [a,t] = max(var(allresp));
    plot(xvals,allresp(:,t));
    result.peakresp = allresp(:,t);
    result.xvals = xvals;
    result.allresp = allresp;
else
    for j = 1:size(allresp,1)
    plot(allresp(j,:),'color',colors{j});
    hold on;
    end
    plot([0 size(all,1)],[mean(y) mean(y)],'color',colors{j+1});
end



function MakeFigure(fig)
colors = 'rgbmcyk';
if fig == 1  %effect of exponent, integration time, and noise, on resp asymmetry
   a{1} = simplerc('tsd',0.5,'plotp');
   labels{1} = '\tau = 0.5';
   a{2} = simplerc('tsd',2,'plotp');
   labels{2} = '\tau = 4';
   a{3} = simplerc('tsd',0.5,'expon',2,'plotp');
   labels{3} = '\tau = 0.5,\gamma=2';
   a{4} = simplerc('tsd',2,'plotp','expon',2);
   labels{4} = '\tau = 4,\gamma=2';
   a{5} = simplerc('tsd',0.50,'plotp','expon',2,'noise',8);
   labels{5} = '\tau = 4,\gamma=2,noise =2';
   a{6} = simplerc('tsd',2,'plotp','expon',2,'frpt',4);
   labels{6} = '\tau = 4,\gamma=2,40ms';
   a{7} = simplerc('tsd',2,'plotp','expon',1,'frpt',4,'thresh',4);
   labels{7} = '\tau = 4,\gamma=1,threshold,40ms';
   GetFigure('Peak Resps');
   hold off;
   for j = 1:length(a)
        plot(a{j}.xvals,a{j}.peakresp,colors(j));
        a{j}.ratio = (a{j}.peakresp(11)-a{j}.peakresp(6))./(a{j}.peakresp(6)-a{j}.peakresp(1));
        hold on;
        text(9,max(a{j}.peakresp),sprintf('%.2f',a{j}.ratio));
   end
    legend(labels,'position',[0 0.7 0.3 0.3]);
   GetFigure('Loop Resps');
   hold off; 
   for j = 1:length(a)
       scale = 1./mean(a{j}.allresp(:,1));
       plot(a{j}.allresp(1,:).*scale,a{j}.allresp(end,:).*scale,colors(j))
       hold on;
   end
   title('Normalized by mean sdf rate');
elseif fig == 2 %
end
