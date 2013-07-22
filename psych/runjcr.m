hold off;
if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/anal/psych/dt/data/'];
else
    prefix = '\\lsr-bgc1\bgc/bgc/anal/psych/dt/data/';
end
if ~exist('nresample','var')
    nresample = 0;
end

summ = [];
colors = mycolors;
PSYCH.opt.forceread = 1;
PSYCH.opt.xmin = -0.1;
PSYCH.opt.xmax = 0.2;
labels = {};
h = [];
allblocks = 4:10;
zeroblocks = 2:10;
[Data, blocklist] = readpsychsum([prefix 'jcrdtxp1']);
allblocks = blocklist.blocks(4:end);
xo = max([Data.expno]);

[bdata, blocklist]  = readpsychsum([prefix 'jcrdtxp0']);
zeroblocks = blocklist.blocks(15:end);
for j = 1:length(bdata)
  bdata(j).expno = bdata(j).expno + xo;
end
Data = [Data bdata];

figure(1);

nfit = 1;

if ~isempty(Data)
  for j = 1:length(Data)
      px = findstr(Data(j).name,'sd=');
      if ~isempty(px)
	Data(j).delay = sscanf(Data(j).name(px(end)+3:end),'%d');
      else
	Data(j).delay = NaN;
      end
      px = findstr(Data(j).name,'xo=');
      if ~isempty(px)
	Data(j).xo = sscanf(Data(j).name(px(end)+3:end),'%d');
      else
	Data(j).xo = NaN;
      end
  end
  idx = find([Data.n] >= 10);
  Data = Data(idx);
  for ex = unique([Data.expno]);
    fit = fitpsf(Data,'expno',ex,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax,'nmin',30,'sderr');
    if ~isnan(fit.fit(1))
      h(nfit) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color', ...
		    colors{nfit},'shown');
      labels{nfit} = sprintf('%s %.3g',fit.data(1).name,abs(fit.fit(2)));
      summ(nfit).delay = fit.data(1).delay;
      summ(nfit).xo = fit.data(1).xo;
      summ(nfit).sd = abs(fit.fit(2));
      summ(nfit).pse = fit.fit(1);
      summ(nfit).sdlim = abs(fit.sdlim);
      summ(nfit).data = fit.data;
      summ(nfit).exit = fit.exit;
      nfit = nfit+1;
    end
  end
end

legend(h,labels,3);
figure(2);
nplot = 1;
labels = {};
[a,b] = sort([summ.delay]);
sumdata = summ(b);

    b = find([summ.xo] < 0 & [summ.delay] == 1) ;
    c = find([summ.xo] < 0 & [summ.delay] == -1) ;
    a = find([summ.xo] < 0 & [summ.delay] == 0) ;

    d = find([summ.xo] > 0 & [summ.delay] == 1) ;
    e = find([summ.xo] > 0 & [summ.delay] == -1) ;
    f = find([summ.xo] > 0 & [summ.delay] == 0) ;

    r = (summ(d).sd-summ(f).sd)/(summ(e).sd-summ(f).sd);
    d = (1-r)/(r+1);
    r = (summ(b).sd-summ(a).sd)/(summ(c).sd-summ(a).sd);
    a = (1-r)/(r+1);

    state.freebase = 0;
id = find([sumdata.xo] > 0)
plot([sumdata(id).delay],1./[sumdata(id).sd],'ro-');
state.logy = 1;
fitl = fitgauss([sumdata(id).delay],1./[sumdata(id).sd],state)
labels{nplot} = sprintf('xo = %.1f dt =%.2f',mean([sumdata(id).xo]),fitl.mean);
hold on;
nplot = nplot+1;
id = find([sumdata.xo] < 0)
plot([sumdata(id).delay],1./[sumdata(id).sd],'go-');
fitr = fitgauss([sumdata(id).delay],1./[sumdata(id).sd],state)
labels{nplot} = sprintf('xo = %.1f dt %.2f',mean([sumdata(id).xo]),fitr.mean);
legend(labels,3);

x = -2:0.1:2;
plot(x,gauss([fitl.mean fitl.sd fitl.amp],x),'r');
plot(x,gauss([fitr.mean fitr.sd fitr.amp],x),'g');

drawnow;

[aData, ablocklist] = readpsychsum([prefix 'jcrdtxp1']);
[aData, bblocklist] = readpsychsum([prefix 'jcrdtxp0']);


figure;
idl = find([sumdata.xo] < 0)
idr = find([sumdata.xo] > 0)
hold off;
PlotSDth(sumdata);


newdata = sumdata;
  for j = 1:length(sumdata)
      for k = 1:length(sumdata(j).data)
          newdata(j).data(k).reresp = binornd(sumdata(j).data(k).n,sumdata(j).data(k).p,1,nresample);
      end
  end

  for nrun = 1:nresample
      
      summ = [];

  for j = 1:length(sumdata)
      for k = 1:length(sumdata(j).data)
          newdata(j).data(k).resp = newdata(j).data(k).reresp(nrun);
      end
  end
rData = newdata;
     
    nfit = 1;
    for j = 1:length(rData);
    fit = fitpsf(rData(j).data);
    if ~isnan(fit.fit(1))
      summ(nfit).delay = fit.data(1).delay;
      summ(nfit).xo = fit.data(1).xo;
      summ(nfit).sd = abs(fit.fit(2));
      summ(nfit).pse = fit.fit(1);
      summ(nfit).data = fit.data;
      summ(nfit).exit = fit.exit;
      nfit = nfit+1;
    end
    end
    id = find([summ.xo] < 0);
    state.logy = 1;
    fit = fitgauss([summ(id).delay],1./[summ(id).sd],state);
    a = fit.mean;
    id = find([summ.xo] > 0);
    ;
    b = fit.mean;
    t(nrun,1) = a;
    t(nrun,2) = b;
    sums{nrun} = summ;
    if(nresample < 100)
      idl = find([summ.xo] < 0);
      idr = find([summ.xo] > 0);
      hold on;
      PlotSDth(summ);
    end
end


return;

nrun = 1;
while nrun <= nresample
  newdata = [];
  summ = [];
  blocks = Bresample(allblocks);
  for j = blocks
    id = find([ablocklist.alldata.block] == j);
    newdata = [newdata ablocklist.alldata(id)];
  end
  blocks = Bresample(zeroblocks);
  for j = blocks
    id = find([bblocklist.alldata.block] == j);
    newdata = [newdata bblocklist.alldata(id)];
  end
    rData = consolidate(newdata);
    for j = 1:length(rData)
      px = findstr(rData(j).name,'sd=');
      if ~isempty(px)
	rData(j).delay = sscanf(rData(j).name(px+3:end),'%d');
      else
	rData(j).delay = NaN;
      end
      px = findstr(rData(j).name,'xo=');
      if ~isempty(px)
	rData(j).xo = sscanf(rData(j).name(px+3:end),'%d');
      else
	rData(j).xo = NaN;
      end
    end
    nfit = 1;
    for ex = unique([rData.expno]);
    fit = fitpsf(rData,'expno',ex,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax,'nmin',20);
    if ~isnan(fit.fit(1))
      summ(nfit).delay = fit.data(1).delay;
      summ(nfit).xo = fit.data(1).xo;
      summ(nfit).flag = fit.exit;
      summ(nfit).sd = abs(fit.fit(2));
      summ(nfit).pse = fit.fit(1);
      summ(nfit).data = fit.data;
      nfit = nfit+1;
    end
    end
    id = find([summ.xo] < 0);
    state.logy = 1;
    fit = fitgauss([summ(id).delay],1./[summ(id).sd],state);
    a = fit.mean;
    id = find([summ.xo] > 0);
    fitl = fit;
    fit = fitgauss([summ(id).delay],1./[summ(id).sd],state);
    b = fit.mean;
    if(fit.exit & nfit == 7 & fit.sd < 1000 & fitl.sd < 1000)
      t(nrun,1) = a;
      t(nrun,2) = b;
      sums{nrun} = summ;
      if(nresample < 100)
	idl = find([summ.xo] < 0);
	idr = find([summ.xo] > 0);
	hold on;
	PlotSDth(summ);
      end
      nrun = nrun+1;
    else
      fprintf('Cant fit run %d exit %d, fits %d\n',nrun,fit.exit,nfit-1);
    end
end


return;

% old calcuation
b = find([summ.xo] < 0 & [summ.delay] == 1) ;
    c = find([summ.xo] < 0 & [summ.delay] == -1) ;
    a = find([summ.xo] < 0 & [summ.delay] == 0) ;
    d = find([summ.xo] > 0 & [summ.delay] == 1) ;
    e = find([summ.xo] > 0 & [summ.delay] == -1) ;
    f = find([summ.xo] > 0 & [summ.delay] == 0) ;
    r = (summ(d).sd-summ(f).sd)/(summ(e).sd-summ(f).sd);
    t(nrun,1) = (1-r)/(r+1);
    r = (summ(b).sd-summ(a).sd)/(summ(c).sd-summ(a).sd);
    t(nrun,2) = (1-r)/(r+1);


nplot = 1;
for delay = unique([sumdata.delay]);
    
id = find([sumdata.delay] == delay);
plot([sumdata(id).xo],1./[sumdata(id).sd],'color',colors{nplot});
labels{nplot} = sprintf('Delay %d',delay);
hold on;
nplot = nplot+1;
end
legend(labels,3);
drawnow;

