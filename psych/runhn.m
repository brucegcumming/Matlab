hold off;
if ~exist('nresample','var')
    nresample = 0;
end
summ = [];

if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/anal/psych/dt/data/'];
else
    prefix = '\\lsr-bgc1\bgc/bgc/anal/psych/dt/data/';
end

colors = mycolors;
PSYCH.opt.forceread = 1;
PSYCH.opt.xmin = -0.3;
PSYCH.opt.xmax = 0.3;
PSYCH.opt.nmin = 30;
labels = {};
h = [];
allblocks = 2:17;
zeroblocks = 2:10;
[Data, ablocklist] = readpsychsum([prefix 'hndtxp0'],'useblocks',5:15);
xo = max([Data.expno]);
[bdata, bblocklist] = readpsychsum([prefix 'hndtxp1'],'useblocks',7:14);
for j = 1:length(bdata)
  bdata(j).expno = bdata(j).expno + xo;
end
Data = [Data bdata];

xo = max([Data.expno]);
[bdata, cblocklist] = readpsychsum([prefix 'hndtxp2']);
twoblocks = cblocklist.blocks(1:end);

for j = 1:length(bdata)
  bdata(j).expno = bdata(j).expno + xo;
end
Data = [Data bdata];

figure(1);

nfit = 1;

sdk = 10;
if ~isempty(Data)
  for j = 1:length(Data)
      px = findstr(Data(j).name,'Hz96');
      if ~isempty(px)
          sdk = 10;
      else
          sdk = 14;
      end
      
      px = findstr(Data(j).name,'sd=');
      if ~isempty(px)
          Data(j).delay = sdk * sscanf(Data(j).name(px(end)+3:end),'%d');
      else
          px = regexp(Data(j).name,'sd[0-9] ');
          if ~isempty(px)
              Data(j).delay = sdk * sscanf(Data(j).name(px+2:end),'%d');
          else
              Data(j).delay = NaN;
          end    
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
    fit = fitpsf(Data,'expno',ex,'xmin',PSYCH.opt.xmin,'xmax',PSYCH.opt.xmax,'nmin',PSYCH.opt.nmin,'sderr');
    if ~isnan(fit.fit(1))
      h(nfit) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color', ...
		    colors{nfit});
      labels{nfit} = sprintf('%s %.3g',fit.data(1).name,abs(fit.fit(2)));
      summ(nfit).delay = fit.data(1).delay;
      summ(nfit).xo = fit.data(1).xo;
      summ(nfit).sd = abs(fit.fit(2));
      summ(nfit).pse = fit.fit(1);
      summ(nfit).sdlim = abs(fit.sdlim);
      summ(nfit).data = fit.data;
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

figure;
idl = find([sumdata.xo] < 0);
idr = find([sumdata.xo] > 0);
PlotSDth(sumdata);
drawnow;

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
      nfit = nfit+1;
    end
    end
    id = find([summ.xo] < 0);
    state.logy = 1;
    fit = fitgauss([summ(id).delay],1./[summ(id).sd],state);
    a = fit.mean;
    id = find([summ.xo] > 0);
    fit = fitgauss([summ(id).delay],1./[summ(id).sd],state);
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

