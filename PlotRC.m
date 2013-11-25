function [x,y,details] = PlotRC(RC, varargin)
%
%PlotRC takes a result file from PlotRevCorAny, and replots the data
%aviods recalculating everything in order to view another plot.
%
% plottypes are:
% default tuning curves for x, one for each y (linestyle) and
% timeslice(color)
% 'ty'  plot sdfs for y(1) vs sdfs for y(2) (AC or psych)
%5 Caluclate Psych kernel CP, and predicted rate difference.
%
% PlotRC(res,'acresp') to build SDF comparisons for Corr, AC.
% PlotRC(res,'acresp','reftype',2) Estimates UC from UC+disps with no resp.
%
%  PlotRC(res, 'sdfs','labela', tag)  uses figure tagged wiht "tag" to plot

smoothw = 50;
smoothsd = 0;
colors = mycolors(0);
times = -1000:100:4000;
step = 1;
nplot = 1;
plotoffset = 0;
showonset = 0;
showvar = 1;
zcheck = 0;
linestyles = {'-', '--', ':','-.'};
symbols = 'nnnnnnnnn';
Expt = [];
RCres = [];
RCs = {};
yvals = [1 2];
timerange = [300 1500];
sumy = 0;
NETSPK = 6;
ACLOOP = 7;
BLANKSDF = 8;
ACRESP = 9;
SDFDIFF = 10;  %- mean rate, to remove framerate
ACDIFF = 11;
DXRESP = 12;
DXDIFF = 13;
BLANKDIFF = 14;
BLANKRATIO = 15;
DIAGONAL = 16;  %for second order plots, look just at main diagonal 
OCULARSDF = 17;  %look at binocular summation
plottype = 0;
resptype = 0;
xid = [];
zid = [];
sumz = 0;
holdon = 0; 
details = [];
coff = 0;
normalize = 0;
framerate = 166.7;
yid = [];
y = [];
details = [];
x = [];
subset = [];
starttime = 0;
showblank = 0;
calcvar = 0;
plotratio = 0;
nexp = 1;
resp = [];

if ~isfield(RC,'delaysamples') %no spike data, just LFP
    x = NaN; y = NaN; h  =0;
    return;
end
if isfield(RC,'subres') 
  subres = [0 1:length(RC.subres)];
else
  subres = [];
end
      
      
if isfield(RC,'slices')
    slices = RC.slices;
else
    slices = RC.delaysamples;  %% ? or times(delaywamples) - check
end
if isfield(RC,'frameperiod')
    framerate = RC.frameperiod;
end
legendpos = 'Best';
plotname = '';
delay = 0;
labela = 'RCplot';
labelb = 'RCplotb';
setcolor = {};
j = 1;
while j  <= nargin -1
    if isstruct(varargin{j}) & isfield(varargin{j},'Trials')
        Expt = varargin{j};
    elseif isstruct(varargin{j}) & isfield(varargin{j},'subsdf') %Plot RC result with multiple subsdfs
        RCres = varargin{j};
    elseif isstruct(varargin{j}) & isfield(varargin{j},'sdfs') %Plot RC result 
        if nexp == 1
            RCs{1} = RC;
        end
        nexp = nexp+1;
        RCs{nexp} = varargin{j};
    elseif strncmpi('+blank',varargin{j},4)
        showblank = 1;
    elseif strncmpi('acdiff',varargin{j},4)
        plottype = ACDIFF;
    elseif strncmpi('acresp',varargin{j},4)
        plottype = ACRESP;
    elseif strcmpi(varargin{j},'cumsum')
        plotname = varargin{j};
    elseif strcmpi(varargin{j},'cmpsd')
        plotname = varargin{j};
    elseif strncmpi('dxresp',varargin{j},4)
        plottype = DXRESP;
    elseif strncmpi('dxdiff',varargin{j},4)
        plottype = DXDIFF;
    elseif strncmpi('ac',varargin{j},2)
        plottype = ACLOOP;
    elseif strncmpi('blankdiff',varargin{j},8)
        plottype = BLANKDIFF;
        resptype = BLANKDIFF;
    elseif strncmpi('blankratio',varargin{j},8)
        plottype = BLANKRATIO;    
        resptype = BLANKRATIO;
    elseif strncmpi('diagonal',varargin{j},4)
        plottype = DIAGONAL;    
    elseif strncmpi('step',varargin{j},3)
     step = varargin{j+1};
     j = j+1;
    elseif strncmpi('best',varargin{j},3)
        [a,b] = max(RC.vars);
        slices = RC.delays(b);
  elseif strncmpi('blanksdf',varargin{j},7)
      plottype = BLANKSDF;
  elseif strncmpi(varargin{j},'coloroff',7)
        j = j+1;
        coff = varargin{j};
  elseif strncmpi('delay',varargin{j},3)
     delay = varargin{j+1};
     j = j+1;
  elseif strncmpi('hold',varargin{j},3)
      holdon = 1;
  elseif strncmpi('labelb',varargin{j},6)
      j = j+1;
      labelb = varargin{j};
  elseif strncmpi('label',varargin{j},3)
      j = j+1;
      labela = varargin{j};
  elseif strncmpi('normalize',varargin{j},4)
      normalize = 1;
  elseif strncmpi('nplot',varargin{j},3)
     plotoffset = varargin{j+1};
     j = j+1;
  elseif strncmpi('netspk',varargin{j},3)
      plottype = NETSPK;
  elseif strncmpi('novar',varargin{j},3)
      showvar = 0;
  elseif strncmpi('noline',varargin{j},3)
      linestyles = {'none' 'none' 'none' 'none' };
      symbols = 'oooooo';
     
  elseif strncmpi('plot',varargin{j},3)
     plottype = varargin{j+1};
     j = j+1;
  elseif strncmpi('pk',varargin{j},2)
     x = PsychKernel(Expt);
     id = find(x.dispvals > -999);
     plot(x.dispvals(id),x.nearcounts(id),'o-');
     hold on;
     plot(x.dispvals(id),x.farcounts(id),'ro-');
     return;
  elseif strncmpi('sdfdiff',varargin{j},5)
      plottype = SDFDIFF;
  elseif strncmpi('sdfslice',varargin{j},3)
      plottype = 5;
      slices = 0;
  elseif strncmpi('sdfs',varargin{j},3)
      plottype = 5;
  elseif strncmpi('slices',varargin{j},3)
     slices = varargin{j+1};
     j = j+1;
  elseif strncmpi('smooth',varargin{j},3)
     smoothw = varargin{j+1};
     j = j+1;
  elseif strncmpi('sdsmooth',varargin{j},5)
     j = j+1;
     smoothsd = varargin{j};
  elseif strncmpi('subres',varargin{j},4)
     j = j+1;
     if isfield(RC,'subres')
         subres = varargin{j};
     else
         subres = [];
     end
  elseif strncmpi('sumy',varargin{j},4)
      sumy = 1;
  elseif strncmpi('sumz',varargin{j},4)
      sumz = 1;
    elseif strncmpi('ty',varargin{j},2)
      plottype = 2;
      if length(varargin) > j & isnumeric(varargin{j+1})
          j = j+1;
          yvals = varargin{j};
      end
  elseif strncmpi('xy',varargin{j},2)
      plottype = 1;
  elseif strncmpi('ocularsdf',varargin{j},2)
      plottype = OCULARSDF;
  elseif strncmpi('psych',varargin{j},2)
      plottype = 3;
  elseif strncmpi('timerange',varargin{j},6)
      j = j+1;
      timerange = varargin{j};
  elseif strncmpi('var',varargin{j},3)
      if strncmpi('varonly',varargin{j},6)
          calcvar = 2;
      else
          calcvar = 1;
      end
  elseif strncmpi('xid',varargin{j},3)
      j = j+1;
      xid = varargin{j};
      showvar =0;
  elseif strncmpi('zid',varargin{j},3)
      j = j+1;
      zid = varargin{j};
    elseif strncmpi('zcheck',varargin{j},3)
% check for rows with missing elements, and eliminate them
% For Dc exps this removes artifacts where the signal values are not
% in the Dc=0 stimulus
       zcheck = 1;
  end
j = j+1;
end

if isfield(RC,'types') && ~isfield(RC,'type')  %old file
    RC.type = RC.types;
end

if length(RCs) > 1
    CompareRCa(RCs, plotname);
    return;
end
if ~isempty(RCres)
    CompareRCs(RC, RCres,plotname);
    return;
end

if ismember(plottype,[OCULARSDF])
    PlotMonocSDF(RC, varargin{:});
    return;
end

if ismember(plottype,[ACLOOP, ACRESP, ACDIFF, DXRESP, DXDIFF])      
    GetFigure(labela);
    details.figa = gcf;
    if ~holdon
        hold off;
    end
      [x, details] = PlotACResp(RC, varargin{:});

      h = details.h;
      labels = details.labels;
end


if calcvar
        S = cat(2,RC.sdfs.s{:});
    for k = 1:size(RC.sdfs.extras)
        S = cat(2,S,RC.sdfs.extras{k}.sdf);
    end
    details.var = var(S,[],2);
   [details.maxvar, details.maxvart] = max(details.var);
   if calcvar ==2
       x = details.var;
       y = details.maxvart;
       return;
   end
end

sumrate = zeros(size(RC.sdfs.s{1}));
sumn = 0;
for k = 1:prod(size(RC.sdfs.s))
    if RC.sdfs.n(k) > 0 & length(RC.sdfs.s{k})
    sumrate = sumrate + (RC.sdfs.s{k} .* RC.sdfs.n(k));
    sumn = sumn +  RC.sdfs.n(k);
    end
end
for k = 1:size(RC.sdfs.extras)

    sumrate = sumrate + RC.sdfs.extras{k}.sdf .* RC.sdfs.extras{k}.n;
    sumn = sumn + RC.sdfs.extras{k}.n;
end
sumrate = sumrate/sumn;
    state.linestyles = linestyles;
    state.colors = colors(1+coff:end);
    state.symbols = symbols;
    state.sumrate = sumrate;
    if strmatch(RC.type{1},{'dx' 'dO' 'dP'})
        state.showuc = 1;
    end

    
if ismember(plottype,[0 1])
    GetFigure(labela);
    if ~holdon
        hold off;
    end
    [x,y, details] = PlotSlices(RC, slices, plottype, state);
end
details.sumrate = sumrate;
if plottype == 2  %plot sdf for y(1) vs y(2), for each x
    for j = 1:size(RC.sdfs.x,1)
        x(:,j) = RC.sdfs.s{j,yvals(1)}(10:end);
        y(:,j) = RC.sdfs.s{j,yvals(2)}(10:end);
        h(j) = plot(x(:,j),y(:,j),'color',colors{j});
        labels{j} = sprintf('%.2f %d,%d',RC.sdfs.x(j,1),RC.sdfs.n(j,yvals(1)),RC.sdfs.n(j,yvals(2)));
        hold on;
    end
    MarkIdentity(gca);
    if isfield(RC,'cp');
        ylabel(sprintf('Y = %.2f, Choice %.2f',RC.sdfs.y(1,yvals(2)),RC.cp.upstim));
        xlabel(sprintf('Y = %.2f, Choice %.2f',RC.sdfs.y(1,yvals(1)),RC.cp.dnstim));
    end
end
if plottype == 3  %psych, averaged over whole impulse
    [a,b,c] = ImpulseSamples(RC,0.5);
    ids = a:b;
    if delay == 0
        delay = RC.bestdelay;
    end
    ds = RC.delaysamples(ids);
    br = [mean(RC.sdfs.extras{2}.sdf(ds)) mean(RC.sdfs.extras{4}.sdf(ds))];
    nullresp = mean(RC.y(:,1,ids),3)-br(1);
    prefresp = mean(RC.y(:,2,ids),3)-br(2);
    h(1) = plot(nullresp,prefresp,'o');
    hold on;
    labels{1} = sprintf('%.1f - %.1f ms',RC.times(ds(1))./10,RC.times(ds(end))./10);
    plot(mean(RC.sdfs.extras{2}.sdf(ds))-br(1),mean(RC.sdfs.extras{4}.sdf(ds))-br(2),'ro');
    plot(0,0,'ro');
    plot(mean(RC.sdfs.extras{1}.sdf(ds))-br(1),mean(RC.sdfs.extras{3}.sdf(ds))-br(2),'go');
    MarkIdentity(gca);
    legendpos = 'NorthWest';
    diffrate = (mean(prefresp) - mean(nullresp)) * 2;
% would be just mean of meanrates for 1 sec trial
    fprintf('Pred diff: %.2f\n', diffrate );
elseif plottype == 4  %psych
    if delay == 0
        delay = RC.bestdelay;
    end
    ds = RC.delaysamples(delay);
    h(1) = plot(RC.y(:,1,delay),RC.y(:,2,delay),'o');
    hold on;
    labels{1} = sprintf('%.1f ms',RC.times(RC.delaysamples(delay))./10);
    plot(RC.sdfs.extras{2}.sdf(ds),RC.sdfs.extras{4}.sdf(ds),'ro');
    plot(RC.sdfs.extras{1}.sdf(ds),RC.sdfs.extras{3}.sdf(ds),'go');
    MarkIdentity(gca);
    legendpos = 'NorthWest';
elseif plottype == 5 & isfield(RC.sdfs,'z') & length(zid) %second order RC 
%if zid is specified, select resps where stimulus at t=0 is in zid, and
%plot resp for each preceding stim. 
    [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    GetFigure(labela);
    if ~holdon
        hold off;
    end
    if sumz
% first calculate first order kernels
        for j = 1:size(RC.sdfs.x,1)
            [x(:,j) xn(j)] = WeightedSum(squeeze(RC.sdfs.s(j,1,:)),squeeze(RC.sdfs.n(j,1,:)));
        end
        id = find(RC.times > RC.times(1)+166.67);
        zmean =  mean(x(:,zid),2); %First order kernel for the stim condition
        sdfall = WeightedSum(x,xn);
        for j = 1:size(RC.sdfs.x,3)
            z(:,j)  = (x(id:end,j) + zmean(1:length(id)));
            z(:,j) = z(:,j) - sdfall(id);
            [y(:,j), n(j)] = WeightedSum(squeeze(RC.sdfs.s(zid,1,j)),squeeze(RC.sdfs.n(zid,1,j)));
            h(j) = plot(RC.times./10,y(:,j),'color',colors{j});
            labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
            hold on;
            plot(RC.times(1:length(id))./10,z(:,j),':','color',colors{j});
            spikediff(j,k) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
        end
        plot(RC.times./10,zmean,'k');
    else
    for j = 1:size(RC.sdfs.x,1)
    for k = 1:length(xid)
        ls = 1+mod(k-1,length(linestyles));
        y(:,j,k) = RC.sdfs.s{j,1,xid(k)};
        h(j) = plot(RC.times./10,y(:,j,k),'color',colors{j},'linestyle',linestyles{ls});
        labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
        hold on;
        spikediff(j,k) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
    end
    end
    end
    nskip = size(y,1) - length(id) +1;
    details.nskip = nskip;
    details.firstorder = x;
    details.sdfall = sdfall;
    details.zmean = zmean;
    details.ymean = WeightedSum(y,n);
    x = y(1:end-nskip+1,:);
    y = z;
elseif plottype == DIAGONAL
    [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    GetFigure(labela);
    tid = find(RC.times>starttime);
    for j = 1:size(RC.sdfs.x,1)
        y(:,j) = RC.sdfs.s{j,1,j}(tid);
        if resptype == BLANKDIFF
            resp(:,j) = y(:,j)-y(:,1);
            h(j) = plot(RC.times(tid)./10,resp(:,j),'color',colors{j});
        elseif resptype == BLANKRATIO
            resp(:,j) = log(y(:,j)./y(:,1));
            h(j) = plot(RC.times(tid)./10,resp(:,j),'color',colors{j});
        else
            resp(:,j) = y(:,j);
            h(j) = plot(RC.times(tid)./10,y(:,j),'color',colors{j});
        end
        labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
        hold on;
        spikediff(j) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
    end
    x = resp;
elseif ismember(plottype,[5 BLANKDIFF]) & isfield(RC.sdfs,'z') %second order RC 
    [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    GetFigure(labela);
        for j = 1:size(RC.sdfs.x,1)
            [x(:,j) xn(j)] = WeightedSum(squeeze(RC.sdfs.s(j,1,:)),squeeze(RC.sdfs.n(j,1,:)));
        end
        id = find(RC.times > RC.times(1)+166.67);
        zmean =  mean(x(:,xid),2); %First order kernel for the stim condition
        sdfall = WeightedSum(x,xn);
        nskip = id(1);
    
    if ~holdon
        hold off;
    end
    if isempty(xid)
        xid = 1:size(RC.sdfs.z,3);
    end
    tid = find(RC.times > starttime);
    if sumz
        nskip = 1;
        for j = 1:size(RC.sdfs.x,1)
            [y(:,j), n(j)] = WeightedSum(squeeze(RC.sdfs.s(j,1,xid)),squeeze(RC.sdfs.n(j,1,xid)));
            if resptype == BLANKDIFF
                resp(:,j) = y(tid,j)-y(tid,1);
            elseif resptype == BLANKRATIO
                resp(:,j) = log(y(tid,j)./y(tid,1));
            else
                resp(:,j) = y(tid,j);
            end
            h(j) = plot(RC.times(tid)./10,resp(:,j),'color',colors{j});
            z(:,j)  = (x(1:length(id),j) + zmean(id));
            z(:,j) = z(:,j) - sdfall(id);
            labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),n(j));
            hold on;
            plot(RC.times(1:length(id))./10,z(:,j),':','color',colors{j});
            spikediff(j,k) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
        end
        y = y(tid,:);
    else
    for j = 1:size(RC.sdfs.x,1)
    for k = 1:length(xid)
        ls = 1+mod(k-1,length(linestyles));
        y(:,j,k) = RC.sdfs.s{j,1,xid(k)};
        if resptype == BLANKDIFF
        resp(:,j,k) = y(tid,j,k)-y(tid,1,k);
        elseif resptype == BLANKRATIO
        resp(:,j,k) = log(y(tid,j,k)./y(tid,1,k));
        else
            resp(:,j,k) = y(tid,j,k);
        end
        h(j) = plot(RC.times(tid)./10,resp(:,j,k),'color',colors{j},'linestyle',linestyles{ls});
        labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
        hold on;
        spikediff(j,k) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
    end
    end
    end
    x = y(1:end-nskip+1,:);
    details.n = xn;
    if sumz
    y = z;
    end
elseif ismember(plottype,[5 SDFDIFF BLANKDIFF]) && (isempty(subres) || sum(subres == 0)) %%Basic RC plot - sdfs for each condition 
   legendpos = [0.7 0.6 0.3 0.4];
   [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    tid = find(RC.times > starttime);
    details.tid = tid;
    details.figa = GetFigure(labela);
    if ~holdon
        hold off;
    end
    msdf = WeightedSum(RC.sdfs.s,RC.sdfs.n);
    if plottype == BLANKDIFF
        msdf = RC.sdfs.extras{1}.sdf;
    end
    if smoothsd
        msdf = smooth(msdf,smoothsd,'gauss');
    end
    if ndims(RC.sdfs.s) == 3  %second order plot
        yid = 1:size(RC.sdfs.s,3);
        secondorder = 1;
    else
        yid = 1:size(RC.sdfs.s,2);
        secondorder = 0;
    end
    if isempty(xid)
        xid = 1:size(RC.sdfs.s,1);
    end
        for k = yid;
            for j = xid;
                if ~isempty(RC.sdfs.s{j,k})
                    if smoothsd
                        y(:,j,k) = smooth(RC.sdfs.s{j,k}(tid),smoothsd,'gauss');
                    else
                        if secondorder
                            y(:,j,k) = RC.sdfs.s{k,1,j}(tid);
                        else
                            y(:,j,k) = RC.sdfs.s{j,k}(tid);
                        end
                    end
                    if plottype == SDFDIFF || plottype == BLANKDIFF
                        if secondorder & plottype == BLANKDIFF;
                            resp(:,j,k) = y(:,j,k)-y(:,j,1);
                        else
                            resp(:,j,k) = y(:,j,k)-msdf(tid);
                        end
                        h(j) = plot(RC.times(tid)./10,resp(:,j,k),'color',colors{j},'linestyle',linestyles{k});
                    elseif plottype == 5
                        h(j) = plot(RC.times(tid)./10,y(:,j,k),'color',colors{j},'linestyle',linestyles{k});
                    end
                    labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
                    hold on;
                    spikediff(j,k) = mean(RC.sdfs.s{j,k}(atime:btime) - sumrate(atime:btime));
                end
                details.n(j,k) = RC.sdfs.n(j,k);
            end
            nx = j;
        end
    ny = k;
    for j = 1:length(RC.sdfs.extras)
        h(j+nx) = plot(RC.times(tid)./10,RC.sdfs.extras{j}.sdf(tid),':','color',colors{j});
        spikediff(j+nx) = mean(RC.sdfs.extras{j}.sdf(atime:btime) - sumrate(atime:btime));
        labels{j+nx} = sprintf('%s n=%d',RC.sdfs.extras{j}.label,RC.sdfs.extras{j}.n);
        if strmatch(RC.type{1},{'dx' 'dO'}) & RC.sdfs.extraval(j) == -1005 %uncorr
            set(h(j+nx),'color','k','linewidth',2,'linestyle','-');
            state.showuc = 1;
        end
    end
    if sumz & isfield(RC,'z')
        hold off;
        resp = [];
        for j = 1:size(RC.sdfs.x,1)
            [x(:,j) xn(j)] = WeightedSum(squeeze(RC.sdfs.s(j,1,:)),squeeze(RC.sdfs.n(j,1,:)));
            resp(:,j) = x(tid,j) - msdf(tid);
            h(j) = plot(RC.times(tid)./10,resp(:,j),'color',colors{j});
            hold on;
            labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),xn(j));
        end
    end    
    details.h = h;
    spikediff = spikediff * (RC.times(btime) - RC.times(atime))/10000;
    if plottype ~= SDFDIFF
        plot(RC.times(tid)./10,sumrate(tid),'k--');
    end
    yvar = squeeze(var(y,[],2));
    if smoothsd == 0
        yvar = smooth(yvar,5,'gauss');
    end
    if showvar
        plot(RC.times(tid)./10,yvar .* max(y(:)./max(yvar(:))),'k','linewidth',2);
    else
        [a,b] = max(mean(yvar,2));
        plot([RC.times(tid(b))./10 RC.times(tid(b))./10],[min(y(:)) max(y(:))],':');
    end
    details.smoothsd = 5;
    details.varsmooth = yvar;
    if ~isempty(Expt)
     x = PsychKernel(Expt,'vals',RC.sdfs.x);
     x.spikediff = spikediff;
     id = find(x.dispvals > -999);
     diffspk(1) = sum(x.diff(id) .* spikediff(1:nx));
     diffspk(2) = mean([Expt.Trials(x.rid{1}).Count]) - mean([Expt.Trials(x.rid{2}).Count]);
     GetFigure(labelb);
     [y,z] = smhist(x.dxresp(x.rid{1}));
     plot(z,y);
     hold on;
     [y,z] = smhist(x.dxresp(x.rid{2}));
     plot(z,y,'r');
     title(sprintf('CP from kernel %.3f',x.cp));
     hold off;
     y = diffspk;
    end
    if slices(1) == 0
        details.figb = GetFigure(labelb);
        if ~holdon
        hold off; %need hold on if doing this for a subres 
        end
        [a,b] = max(mean(yvar,2));
        slices = RC.times(tid(b));
        details.bestslice = tid(b);
        PlotSlices(RC, slices, 0, state);
    end

elseif plottype == BLANKSDF
    [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    GetFigure(labela);
    k = 1;
    for j = 1:length(RC.sdfs.extras)
        if strcmp(RC.sdfs.extras{j}.label,'Blank')
            if normalize
                h(k) = plot(RC.times./10,RC.sdfs.extras{j}.sdf./sumrate,'-','color',colors{k+coff});
            else
                h(k) = plot(RC.times./10,RC.sdfs.extras{j}.sdf,'-','color',colors{k+coff});
            end
        hold on;
        spikediff(k) = mean(RC.sdfs.extras{j}.sdf(atime:btime) - sumrate(atime:btime));
        labels{k} = sprintf('%s n=%d',RC.sdfs.extras{j}.label,RC.sdfs.extras{j}.n);
        details.n(k) = RC.sdfs.extras{j}.n;
        k = k+1;
        end
    end
    x = RC.sdfs.extras{1};
    details.h = h;
    details.labels = labels;
    y = [];
    if normalize == 0
    plot(RC.times./10,sumrate,'k--');
    end

elseif plottype == NETSPK
    if zcheck
        id = find(prod(RC.sdfs.n') > 0);
    else
        id = 1:size(RC.sdfs.x,1);
    end
    details.xv = RC.sdfs.x(id,1);
    for k = 1:size(RC.sdfs.y,2)
    for j = 1:length(id)
        if isempty(RC.sdfs.s{id(j),k})
        y(:,j,k) = 0;
        elseif smoothsd
            y(:,j,k) = smooth(RC.sdfs.s{id(j),k},smoothsd,'gauss');
        else
            y(:,j,k) = RC.sdfs.s{id(j),k};
        end
        details.n(j,k) = RC.sdfs.n(id(j),k);
    end 
    end
    meanrate = squeeze(mean(mean(y,3),2));
    if smoothsd
        meanrate = smooth(meanrate,smoothsd,'gauss');
    end
    if isnan(timerange(1))  %determine range automatically
        yvar = var(y,[],2);
        [a,mi] = max(yvar);
        [a,b] = min(yvar(1:mi));
        timerange(1) = RC.times(b);
        [a,b] = min(yvar(mi:end));
        timerange(2) = RC.times(b+mi-1);
    end
    tid = find(RC.times > timerange(1) & RC.times < timerange(2));
    for k = 1:size(RC.sdfs.y,2)
    for j = 1:length(id)
        netspk(j,k) = sum(y(tid,j,k) - meanrate(tid));
    end 
    end
    netspk = netspk * mean(diff(RC.times))./10000;
    netspk = netspk + mean(sumrate(20:end)) * framerate/10000;
    for k = 1:size(RC.sdfs.extras,1)
        details.ex(k) = sum(RC.sdfs.extras{k}.sdf(tid) - meanrate(tid)) * mean(diff(RC.times)./10000);
        details.ex(k) = details.ex(k)+ mean(sumrate(20:end)) * framerate/10000;
        details.exval(k) = RC.sdfs.extraval(k);
    end
    x = netspk;
    y = details;
elseif plottype == BLANKRATIO
      legendpos = [0.7 0.6 0.3 0.4];
   [a,atime] = min(abs(RC.times - 200));
    [a,btime] = min(abs(RC.times - 1000));
    tid = find(RC.times > starttime);
    GetFigure(labela);
    msdf = RC.sdfs.extras{1}.sdf;
    if smoothsd
        msdf = smooth(msdf,smoothsd,'gauss');
    end
    yid = 1:size(RC.sdfs.s,2);
    if isempty(xid)
        xid = 1:size(RC.sdfs.s,1);
    end
    for k = yid;
    for j = xid;
        if ~isempty(RC.sdfs.s{j,k})
        if smoothsd
            y(:,j,k) = smooth(RC.sdfs.s{j,k}(tid),smoothsd,'gauss');
        else
            y(:,j,k) = RC.sdfs.s{j,k}(tid);
        end
        resp(:,j,k) = log(y(:,j,k)./msdf(tid));
        h(j) = plot(RC.times(tid)./10,resp(:,j,k),'color',colors{j},'linestyle',linestyles{k});
        labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),RC.sdfs.n(j,1));
        hold on;
        end
        details.n(j,k) = RC.sdfs.n(j,k);
    end
    nx = j;
    end
    ny = k;
    for j = 2:length(RC.sdfs.extras)
        h(j+nx) = plot(RC.times(tid)./10,RC.sdfs.extras{j}.sdf(tid),':','color',colors{j});
        spikediff(j+nx) = mean(RC.sdfs.extras{j}.sdf(atime:btime) - sumrate(atime:btime));
        labels{j+nx} = sprintf('%s n=%d',RC.sdfs.extras{j}.label,RC.sdfs.extras{j}.n);
        if strmatch(RC.type{1},{'dx' 'dO'}) & RC.sdfs.extraval(j) == -1005 %uncorr
            set(h(j+nx),'color','k','linewidth',2,'linestyle','-');
            state.showuc = 1;
        end
    end
    
    details.h = h;
    yvar = squeeze(var(y,[],2));
    if smoothsd == 0
    yvar = smooth(yvar,5,'gauss');
    end
    if showvar
    plot(RC.times(tid)./10,yvar .* max(y(:)./max(yvar(:))),'k','linewidth',2);
    else
        [a,b] = max(mean(yvar,2));
        plot([RC.times(tid(b))./10 RC.times(tid(b))./10],[min(resp(:)) max(resp(:))],':');
    end
    if sumz
        hold off;
        resp = [];
        for j = 1:size(RC.sdfs.x,1)
            [x(:,j) xn(j)] = WeightedSum(squeeze(RC.sdfs.s(j,1,:)),squeeze(RC.sdfs.n(j,1,:)));
            resp(:,j) = log(x(tid,j)./msdf(tid));
            h(j) = plot(RC.times(tid)./10,resp(:,j),'color',colors{j});
            hold on;
            labels{j} = sprintf('%.2f n=%d',RC.sdfs.x(j,1),xn(j));
        end
    end    


    if ~isempty(Expt)
     x = PsychKernel(Expt,'vals',RC.sdfs.x);
     id = find(x.dispvals > -999);
     diffspk(1) = sum(x.diff(id) .* spikediff(1:nx));
     diffspk(2) = mean([Expt.Trials(x.rid{1}).Count]) - mean([Expt.Trials(x.rid{2}).Count]);
     GetFigure(labelb);
     [y,z] = smhist(x.dxresp(x.rid{1}));
     plot(z,y);
     hold on;
     [y,z] = smhist(x.dxresp(x.rid{2}));
     plot(z,y,'r');
     title(sprintf('CP from kernel %.3f',x.cp));
     hold off;
     y = diffspk;
    end
    if slices(1) == 0
        GetFigure(labelb);
        if ~holdon
        hold off; %need hold on if doing this for a subres 
        end
        [a,b] = max(mean(yvar,2));
        slices = RC.times(tid(b));
        details.bestslice = tid(b);
        PlotSlices(RC, slices, 0, state);
    end 
end

if ismember(resptype,[BLANKDIFF BLANKRATIO])
if size(resp,2) == 3
    plot(RC.times(tid)./10,resp(:,1)+resp(:,3),'k:');
    yl = get(gca,'ylim');
    ratio = (resp(:,1)+resp(:,3))./resp(:,2);
    plot(RC.times(tid)./10,ratio .* yl(2)/2,'k');
    set(gca,'ylim',yl);
elseif size(resp,2) == 4
    yl = get(gca,'ylim');
    for j = 1:size(resp,3)
    plot(RC.times(tid)./10,resp(:,2,j)+resp(:,4,j),'k:');
    if plotratio
    ratio = (resp(:,2,j)+resp(:,4,j))./resp(:,3,j);
    plot(RC.times(tid)./10,ratio .* yl(2)./2,'k');
    end
    end
    set(gca,'ylim',yl);
end
end
if ~isempty(resp)
    details.resp = resp;
end

details.bestms = RC.delays(RC.bestdelay)./10;

if exist('h')
mylegend(h,labels,'Location',legendpos);
end


if ~isempty(subres) & ~isfield(details,'donesub') & isfield(RC,'subres');
    if isfield(details,'h') & ishandle(details.h)
        h = details.h;
    else
    h = [];
    end
    if isfield(details,'dxid')
        args = {'dxid' details.dxid details.refid 'hold'};
    else
        args = {'hold'};
    end
   k =0;
   id = find(subres > 0);
   sdfs{1} = x;
   nsub = 1;
   for j = id
       if sum(subres == 0) %aready plotted top one.
       hold on;
       end
       [a,b,c] = PlotRC(RC.subres{subres(j)},varargin{:},'coloroff',j,args{:});
       hold on;
       if plottype == 0 & length(slices) == 1
           labels{j+k} = sprintf('%.2f %.0fms n=%.0f',RC.subres{subres(j)}.(RC.ctype),c.bestms,mean(c.n(:)));
       elseif isfield(c,'n')
           labels{j+k} = sprintf('%.2f n=%.0f',RC.subres{subres(j)}.(RC.ctype),mean(c.n(:)));
       else
           labels{j+k} = sprintf('%.2f n=0',RC.subres{subres(j)}.(RC.ctype));
           fprintf('No trials in %s\n',RC.name);
       end
       nsub = nsub+1;
       sdfs{nsub} = a;
       
       if isfield(c,'h')
       h(j+k) = c.h(1);
       end
   end
   if isempty(details) || ~isfield(details,'n');
       details = c;
   end
   if length(sdfs) > 1
   details.subsdf = sdfs;
   end
   if plottype == 0 & length(slices) == 1
       labels{1} = sprintf('%s= %.2f %.0fms n=%.0f',RC.ctype,RC.(RC.ctype),details.bestms,mean(details.n(:)));
   elseif ~isfield(details,'n')
       labels{1} = sprintf('%s = %.2f, n = 0',RC.ctype, RC.(RC.ctype));
   else
       labels{1} = sprintf('%s = %.2f, n = %.0f',RC.ctype, RC.(RC.ctype),mean(details.n(:)));
   end
   if length(h) == length(labels);
   mylegend(h,labels);
   end
   title(MakeTitle(RC));
end

   
if plottype == 1 
    MarkIdentity(gca);
end


function s = MakeTitle(rc, varargin)
s = '';
if isfield(rc,'name')
    [a,s] = fileparts(rc.name);
    s = strrep(s,'\','/');
end
if isfield(rc.Header,'cellnumber')
    s = sprintf('%s Cell%.0f (P%.1f)',s,rc.Header.cellnumber,rc.Header.probe);
elseif isfield(rc,'probe')
    s = sprintf('%s P%.0f',s,rc.probe);
end


function [startt, endt, maxt] = ImpulseSamples(res, crit)
%Caclulate time samples that are in the impulse resp
[maxv, maxt] = max(res.vars);
minv = min(res.vars);
t = maxt;
while(res.vars(t) > minv + (maxv-minv)  * crit)
    t = t+1;
end
endt = t;
t = maxt;
while(res.vars(t) > minv + (maxv-minv)  * crit)
    t = t-1;
end
startt = t;

function [msdf, details] = CalcMeanSdf(res, id)

details.n = 0;
msdf =[];
nsdf = 0;

if nargin == 1
        msdf = res.sdfs.s{1};
        nsdf = 1;
        allsdf(1,:) = res.sdfs.s{1};
        for j = 2:prod(size(res.sdfs.s))
            if ~isempty(res.sdfs.s{j})
                msdf = msdf + res.sdfs.s{j};
                allsdf(j,:) = res.sdfs.s{j};
                nsdf = nsdf+1;
            end
        end
        msdf = msdf ./nsdf;
else
    msdf = zeros(size(res.sdfs.s(id(1))));
    for j = 1:length(id)
       msdf = msdf + res.sdfs.s{id(j)};
       nsdf = nsdf+1;
    end    
    msdf = msdf ./nsdf;
end
details.n = nsdf;

function [acres, details] = PlotDXResp(res, varargin)

secondorder = 0;
plottype = 0;
colors = mycolors;
coff = 0;
bsdf = [];
reftype = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'acloop',5)
        plottype = 1;
    elseif strncmpi(varargin{j},'acdiff',5)
        plottype = 2;
    elseif strncmpi(varargin{j},'coloroff',7)
        j = j+1;
        coff = varargin{j};
    elseif strncmpi(varargin{j},'reftype',7)
        j = j+1;
        reftype = varargin{j};
    end
    j = j+1;
end
yvals = res.y(:,:,res.bestdelay);

function rc = CleanRC(rc)



good = ones(size(rc.sdfs.s,1),1);
for j = 1:size(rc.sdfs.s,1)
for k = 1:size(rc.sdfs.s,2)
    if isempty(rc.sdfs.s{j,k})
        good(j) = 0;
    end
end
end
if sum(good(:) == 0)
    good = find(good);
rc.sdfs.s = rc.sdfs.s(good,:);
rc.sdfs.n = rc.sdfs.n(good,:);
rc.sdfs.x = rc.sdfs.x(good,:);
rc.sdfs.y = rc.sdfs.y(good,:);
rc.y = rc.y(good,:,:);
rc.x = rc.x(good,:,:);
end

function [acres, details] = PlotACResp(res, varargin)

secondorder = 0;
plottype = 0;
colors = mycolors;
coff = 0;
bsdf = [];
reftype = 0;
dxid = [];
dxsigns = [];
refid = 0;
showblank = 0;
prettify = 0;
j = 1;

res = CleanRC(res);
while j <= length(varargin)
    if strncmpi(varargin{j},'acloop',5)
        plottype = 1;
    elseif strncmpi(varargin{j},'acdiff',5)
        plottype = 2;
    elseif strncmpi(varargin{j},'+blank',5)
        showblank=1;
    elseif strncmpi(varargin{j},'dxid',4)
        j = j+1;
        dxid = varargin{j}(1,:);
        dxsigns = varargin{j}(2,:);
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            refid = varargin{j};
        end
    elseif strncmpi(varargin{j},'dxdiff',5)
        plottype = 4;
    elseif strncmpi(varargin{j},'dxresp',5)
        plottype = 3;
    elseif strncmpi(varargin{j},'pretty',5)
        prettify = 1;
    elseif strncmpi(varargin{j},'coloroff',7)
        j = j+1;
        coff = varargin{j};
    elseif strncmpi(varargin{j},'reftype',7)
        j = j+1;
        reftype = varargin{j};
    end
    j = j+1;
end
details.snr = 0;

         yvals = res.y(:,:,res.bestdelay);
% take disparities that differ by more than 1 SD from mean of all
        [a, tval] = min(abs(res.times - res.delays(res.bestdelay)));
        endvals = res.y(:,:,end);
        [msdf,a] = CalcMeanSdf(res);
        if ismember(plottype, [3 4]) %just the dx resp
            details.nsdf(3) = a.n;
            pid = find(yvals > prctile(yvals(:),80));
            [xp, a] = CalcMeanSdf(res, pid);
            details.nsdf(1) = a.n;
            nid = find(yvals < prctile(yvals(:),20));
            xn = CalcMeanSdf(res, nid);
            details.nsdf(2) = a.n;
            details.n = res.sdfs.n([pid; nid]);

            acres.sdf(1,:) = xp;
            acres.sdf(2,:) = xn;
            acres.sdf(3,:) = msdf;
            if plottype == 3
            details.h(1) = plot(xp);
            hold on;
            details.h(2) = plot(xn,'r');
            details.labels{1} = sprintf('Pref');
            details.labels{2} = sprintf('Null');
            elseif plottype == 4
                np = 1;
                details.h(np) = plot(xp-xn,'color',colors{1+coff});
                details.labels{np} = sprintf('Pref-Null %.0f',sum(details.n(:)));
                np = np+1;
                if isfield(res,'subres')
                    hold on;
                    [a,b] = PlotACResp(res.subres{1},'dxdiff','coloroff',1);
                    hold on;
                    details.donesub = 1;
                    details.h(np) = b.h(1);
                    details.labels{np} = b.labels{1};
                    details.n = [details.n; b.n];
                    details.subres{np} = a;
                    np = np+1;
                    details.h(np) = plot(msdf-a.sdf(3,:)','k');
                    details.labels{np} = sprintf('High - Low %.0f',sum(details.n(:)));
                end
            end
            return;
        end
 % or use deviation from UC, based on SD at end
 
        if reftype == 1 || (reftype == 2 && length(refid) && sum(refid))
 %if refid is empty, it means this has been forced by an input argument
 %if refid has entries, it means ref disps already chosen.
            if ~isempty(refid) & refid == 0
            s = std(cat(1,res.y(:,1,:),res.y(:,2,:)));
            id = find(diff(yvals,[],2)< s(1));
            else
                id = refid;
            end
            usdf = zeros(size(res.sdfs.s{id(1),1}));
            for j = 1:length(id)
                usdf = usdf + res.sdfs.s{id(j),1};
                usdf = usdf + res.sdfs.s{id(j),2};
            end
            uid = find(res.sdfs.extraval == -1005);
            x.refid = id;
            if length(uid)
                usdf = usdf + res.sdfs.s{uid(1),1};
                usdf = usdf ./(1+2*length(id));
            else
            usdf = usdf ./(2*length(id));
            end
            uc = usdf(tval);
        elseif isempty(res.sdfs.extras) || sum(find(res.sdfs.extraval == -1005)) == 0
            uc = mean(yvals(:));
            usdf = msdf;
        else
            uc = res.sdfs.extras{1}.sdf(tval);
            usdf = res.sdfs.extras{1}.sdf;
            if reftype == 2
                yrange = uc + [-2 2] .* std(endvals(:));
                if ~isempty(refid) & refid == 0
                    id = find(yvals(:,1) > yrange(1) & yvals(:,1) < yrange(2) & yvals(:,end) > yrange(1) & yvals(:,end) < yrange(2));
                    k = 0;
                for j = 1:length(id)
                    for k = 1:size(res.sdfs.s,2)
                    usdf = usdf + res.sdfs.s{id(j),k};
                    end
                end
                usdf = usdf ./(1+k*length(id));
                else
                    id = refid;
                    for j = 1:length(id)
                        usdf = usdf + res.sdfs.s{id(j),1};
                        usdf = usdf + res.sdfs.s{id(j),2};
                    end
                usdf = usdf ./(1+2*length(id));
                end
                details.refid = id;
            else
               details.refid = []; 
            end
        end
        
        id = find(res.sdfs.extraval == -1009);
        if length(id) == 1
            bsdf = res.sdfs.extras{id}.sdf;
        end
        yrange = uc + [-2 2] .* std(endvals(:));
        zscores = (yvals-uc)./std(endvals(:));
        allid = find(abs(zscores) > 2 & sum(res.sdfs.n,3) > res.nmin);
%        pid  = find(yvals > yrange(2) & sum(res.sdfs.n,3) > res.nmin);
%        nid  = find(yvals < yrange(1) & sum(res.sdfs.n,3) > res.nmin);
%        allid = [pid; nid];
        [l,m] = ind2sub(size(yvals),allid);
        if length(unique(l)) > size(yvals,1)*0.7
            crit = max(prctile(abs(zscores),50));
            allid = find(abs(zscores) > crit & sum(res.sdfs.n,3) > res.nmin);
            [l,m] = ind2sub(size(yvals),allid);
        end
        if isempty(dxid)
        dxid = unique(l)';
        else
            allid = sub2ind(size(yvals),[dxid dxid],[ones(size(dxid)) 2 * ones(size(dxid))]);
        end
        acres.did = dxid;
        
        npairs = floor(size(res.sdfs.s,2)/2);
        if npairs == 0 || isempty(dxid)
            details.h = 0;
            details.labels = [];
            return;
        end
        lsdf = zeros(size(res.sdfs.s{allid(1)},1),npairs);
        csdf = lsdf;
        asdf = lsdf;
        na(1:npairs) = 0;
        nc(1:npairs) = 0;
        if isempty(dxsigns)
            for j = 1:length(dxid)
                if yvals(dxid(j),end) <  yrange(1) || (yvals(dxid(j),1)> yrange(2) && yvals(dxid(j),end) < yrange(2))
                    dxsigns(j) = -1;
                else
                    dxsigns(j) = 1;
                end
            end
        end
        for j = 1:length(dxid)
            for k = 1:npairs
            l = dxid(j);
            c = size(res.sdfs.s,2)+1-k;
            a = k;
            if dxsigns(j) == 1
                lsdf(:,k) = lsdf(:,k) + res.sdfs.s{l,c,1}-usdf;
                lsdf(:,k) = lsdf(:,k) - res.sdfs.s{l,a,1}+usdf;
                csdf(:,k) = csdf(:,k) + res.sdfs.s{l,c,1}-usdf;
                asdf(:,k) = asdf(:,k) + res.sdfs.s{l,a,1}-usdf;
                if secondorder
                    lsdf(:,k) = lsdf(:,k) + res.sdfs.s{l,m,2};
                end
            else
                lsdf(:,k) = lsdf(:,k) + res.sdfs.s{l,a,1}-usdf;
                lsdf(:,k) = lsdf(:,k) - res.sdfs.s{l,c,1}+usdf;
                csdf(:,k) = csdf(:,k) - res.sdfs.s{l,c,1}+usdf;
                asdf(:,k) = asdf(:,k) - res.sdfs.s{l,a,1}+usdf;
                if secondorder
                    lsdf(:,k) = lsdf(:,k) - res.sdfs.s{l,m,2};
                end
            end
                na(k) = na(k) + res.sdfs.n(1,a,1);
                nc(k) = nc(k) + res.sdfs.n(1,c,1);
            end
        end
        csdf = csdf./length(dxid);
        asdf = asdf./length(dxid);
        lsdf = lsdf./length(dxid);
        acres.csdf = csdf;
        acres.asdf = asdf;
        acres.lsdf = lsdf;
        acres.msdf = msdf;
        acres.usdf = usdf;
        if exist('bsdf')
            acres.bsdf = bsdf;
        end
        details.n = [na nc];
        details.dxid = [dxid; dxsigns];
        details.refid = refid;
        n = length(dxid);
        if plottype ==1
            tid = find(res.times > 400);
            dxstr = sprintf('%.2f ',res.x(dxid,1));
            for k = 1:size(asdf,2)
            a = plot(csdf(:,k),-asdf(:,k),'color',colors{k+coff});
            h(k) = a(1);
            [a,b] = max(lsdf(tid,k));
            b = b+tid(1)-1;
            id = find(lsdf(1:b,k) < a/2); %id(end) is half max on rising edged
            if isempty(id)
                ht(1) = tid(1);  %40ms
            else
                ht(1) = id(end);
            end
            id = find(lsdf(b:end,k) < a/2);
            ht(2) = b+ht(1)-1;
            hold on;
            plot(csdf(1:b,k),-asdf(1:b,k),'color',colors{k+coff},'linewidth',2);
%            plot(csdf(1:b,k),asdf(1:b,k),'k.');
            if npairs > 1
                labels{k} = sprintf('Corr  %1f',abs(res.yvals(k)));
            else
                labels{k} = ['C' dxstr];
            end
            end
            yl = get(gca,'ylim');
            plot([0 yl(2)],[0 yl(2)],'k');
            
            
        elseif plottype == 2
            h(1) = plot(res.times./10,csdf-asdf);
            hold on; plot([res.times(1)./10 res.times(end)./10],[0 0],'--');
            labels{1} = 'C-A';
        else
            nc = 1;
            if prettify
                nc = 2; %different colors for C,AC
                labels{1} = 'Correlated';
                labels{2} = 'Anticorrelated';
            else
                labels{1} = 'C';
                labels{2} = 'A';
            end
            for k = 1:size(asdf,2)
                a(k) = plot(res.times./10,csdf(:,k),'color',colors{coff+k});
                hold on;
            end
            h(1) = a(1);
            
        hold on;
        for k = 1:size(asdf,2)
            if prettify
                a(k) = plot(res.times./10,-asdf(:,k),'-','color',colors{nc+coff+k-1});
            else
                a(k) = plot(res.times./10,-asdf(:,k),'--','color',colors{nc+coff+k-1});
            end
        end
        h(2) = a(1);
       
%        plot(lsdf,'g');
tid = find(res.times > 100);
        if prettify ~= 1
        h(3) = plot(res.times(tid)./10,usdf(tid)-mean(usdf),':','color',colors{1+coff});
            labels{3} = 'U';
        end
            if length(bsdf) & showblank
                plot(bsdf-usdf);
            end
        end
        sm = smooth(csdf,10);
        details.noise = std(csdf-sm);
        details.snr = max(sm)./details.noise;
        details.h = h;
        details.labels= labels;
        if prettify
            xlabel('Time (ms)');
            ylabel('Rate (spikes/sec)');
        end
function kernel = PsychKernel(Expt, varargin)

dvals = [];
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'vals',3)
        j = j+1;
        dvals = varargin{j};
    end
    j = j+1;
end
details.h = h;
nearid = find([Expt.Trials.RespDir] < 0);
farid = find([Expt.Trials.RespDir] > 0);

if isempty(dvals)
    dvals = unique([Expt.Trials.Pd]);
end
neardisps = [Expt.Trials(nearid).Pd];
fardisps = [Expt.Trials(farid).Pd];
% calculate mean # frames of each val per trial
for j = 1:length(dvals)
    kernel.nearcounts(j) = sum(neardisps(:) == dvals(j))./size(neardisps,2);
    kernel.farcounts(j) = sum(fardisps(:) == dvals(j))./size(fardisps,2);
end
if 0  %% disps < -999 make a mess of this method. Or any disp vals not in dvals.
    kernel.nearcounts = hist(neardisps(:),dvals)./size(neardisps,2);
    kernel.farcounts = hist(fardisps(:),dvals)./size(fardisps,2);
end
kernel.dispvals = dvals;
kernel.diff = kernel.nearcounts - kernel.farcounts;
for j = 1:length(Expt.Trials)
    dxdist(j,:) = hist(Expt.Trials(j).Pd,dvals);
    dxresp(j) = sum(dxdist(j,:) .* kernel.diff);
end
kernel.cp = CalcCP(dxresp(nearid),dxresp(farid));
kernel.dxresp = dxresp;
kernel.rid{1} = nearid; 
kernel.rid{2} = farid;


function [x,y, details] = PlotSlices(RC, slices, plottype, state)
for it = 1:length(slices)
    [a,b] = min(abs(RC.times - slices(it)));
    for k = 1:size(RC.sdfs.x,2)
        for j = 1:size(RC.sdfs.x,1)
            if b <= length(RC.sdfs.s{j,k})
                y(j,k,it) = RC.sdfs.s{j,k}(b);
                x(j,k,it) = RC.sdfs.x(j,k);
                details.n(j,k) = RC.sdfs.n(j,k);
            end
        end
        if plottype == 0
            h(it) = plot(x(:,k,it),y(:,k,it),'linestyle',state.linestyles{k},'color',state.colors{it},'marker',state.symbols(k));
            labels{it} = sprintf('%.0fms',slices(it)/10);
            hold on;
            plot(minmax(x(:,1,1)),[state.sumrate(b) state.sumrate(b)])
            if isfield(state,'showuc')
                id = find(RC.sdfs.extraval == -1005);
                if length(id)
                    plot(minmax(x(:,1,1)),[RC.sdfs.extras{id}.sdf(b) RC.sdfs.extras{id}.sdf(b)],'k-')
                    details.extras.x = -1005;
                    details.extras.y = RC.sdfs.extras{id}.sdf(b);
                end
            end
        end
    end
    details.h = h;

    if plottype ==1
        h(it) = plot(y(:,it,1),y(:,it,2),'o','color',colors{it});
        labels{it} = sprintf('%.0fms',slices(it)/10);
        hold on;
    end
end



function CompareRCa(RCs, plotname)
GetFigure('Compare');
hold off;
t = RCs{1}.times./10;
colors = mycolors;
linestyles = {'-' '--' ':' '-:'};
for exp = 1:length(RCs)
if sum(RCs{exp}.sdfs.extraval == -1009)
    id = find(RCs{exp}.sdfs.extraval == -1009);
    for j = 1:size(RCs{exp}.sdfs.s,1)
        resp(j,exp,:) = RCs{exp}.sdfs.s{j}-RCs{exp}.sdfs.extras{id}.sdf;
        tmp = plot(t,squeeze(resp(j,exp,:)),'-','color',colors{j},'linestyle',linestyles{exp});
        hold on;
        if exp == 1
            h(j) = tmp;
            labels{j} = sprintf('%.1f',RCs{exp}.sdfs.x(j));
        end
    end
else
    for j = 1:size(RCs{exp}.sdfs.s,1)
        resp(j,exp,:) = RCs{exp}.sdfs.s{j};
        tmp = plot(t,squeeze(resp(j,exp,:)),'-','color',colors{j},'linestyle',linestyles{exp});
        if exp == 1
            h(j) = tmp;
            labels{j} = sprintf('%.1f',RCs{exp}.sdfs.x(j));
        end
    end
end
end
for exp = 1:length(RCs)
    if length(plotname) & strmatch(plotname,{'cmpsd' 'cmpSD'})
        h(j+exp) = plot(t,squeeze(std(resp(:,exp,:),1)),'k-','linewidth',2,'linestyle',linestyles{exp});
    else
        h(j+exp) = plot(t,squeeze(mean(resp(:,exp,:),1)),'k-','linewidth',2,'linestyle',linestyles{exp});
    end
    labels{j+exp} = sprintf('%.0fms',RCs{exp}.stimdur./10);
end
    legend(h,labels);
    title(strrep(RCs{1}.name,'_',' '));
    sds = std(resp,1);
    scales = max(sds,[],3);
    sumvar = sum(std(resp,1)./repmat(scales,[1 1 size(sds,3)]),2);
    [vmax, tmax] = max(sumvar);
    plot([t(tmax) t(tmax)], get(gca,'ylim'),'k:');
    GetFigure('RCplotb');
    plot(RCs{1}.x(:,1),squeeze(resp(:,:,tmax)./repmat(scales,[size(resp,1) 1])));
    title(sprintf('RFs at %.1f',t(tmax)));

function CompareRCs(RC, RCres, plotname)
GetFigure('Compare');
hold off; 

ms = RC.times./10;
tdiff = 10;  %need to derive this from data one day...
x = RCres.subsdf{2}.csdf;
y = RCres.subsdf{2}.asdf;
lsum = x(11:end-10) + x(1:end-20) + x(21:end);
asum = y(11:end-10) + y(1:end-20) + y(21:end);

if strcmp(plotname,'cumsum')
plot(ms,cumsum(RCres.subsdf{1}.csdf));
hold on;
plot(ms,cumsum(RCres.subsdf{1}.asdf),':');
plot(ms,cumsum(RCres.subsdf{2}.csdf),'r');
plot(ms,cumsum(RCres.subsdf{2}.asdf),'r:');
tdiff = 10;  %need to derive this from data one day...
plot(ms(21:end),cumsum(lsum),'g');
plot(ms(21:end),cumsum(asum),'g:');
else 
plot(ms,RCres.subsdf{1}.csdf);
hold on;
plot(ms,RCres.subsdf{1}.asdf,':');
plot(ms,RCres.subsdf{2}.asdf,'r:');
plot(ms,RCres.subsdf{2}.csdf,'r');
plot(ms(21:end),lsum,'g');
plot(ms(21:end),asum,'g:');
end    

function PlotMonocSDF(rc, varargin)
%
%might a loop plot be useful, showing summed response(all) (or variance, to
%avoid weighting too much on blank) vs ratio. Showing that fall in ratio
%late is not becuase of rate decline. Perhaps look at ratio at halfmax each
%end of the response?  Smooth heavily to get times first. 

if rc.sdfs.x(1) == -1009
blnk = rc.sdfs.s{1};
xo = 1;
else
blnk = rc.sdfs.extras{1}.sdf;
xo = 0;
end
ClearPlot(gcf);
tm = rc.times./10;
colors = mycolors;
hold off;
for j = 1:size(rc.sdfs.s,1)-xo
    resps(j,:) = rc.sdfs.s{j+xo}-blnk;
    plot(tm,resps(j,:),'color',colors{j});
    hold on;
end

plot(tm,resps(1,:)+resps(3,:),'k');
ratio = resps(2,:)./(resps(1,:)+resps(3,:));
ratio = smooth(ratio,5);
id = find(abs(ratio) < 2);
lax = gca;
rax = AddRPlot(lax, tm(id),ratio(id),'.');
hold on;
plot(minmax(tm(id)),[1 1]);
