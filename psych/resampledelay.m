sums = {};
Trials = [];

allresumm = summ;
redata = summ;
for j = 1:length(summ)
    for k = 1:length(summ(j).data)
        allresumm(j).data(k).resp = binornd(summ(j).data(k).n,summ(j).data(k).p,nresample,1);
    end
end


for res = 1:nresample;

for d = 1:length(summ)
    for k = 1:length(summ(d).data)
        redata(d).data(k).resp = allresumm(d).data(k).resp(res);    
        redata(d).data(k).p = redata(d).data(k).resp / redata(d).data(k).n;    
    end
end

resumm = [];
trialcounts = [];
nplot = 1;
for sd = unique([summ.delay])
  sddata = redata(find([redata.delay] == sd));
  if length(sddata) >0
    fit = fitpsf(sddata.data,'xmin',xl,'xmax',xu,'nmin',nmin);
    if ~isnan(fit.fit(1))
      resumm(nplot).delay = sd;
      resumm(nplot).sd = fit.fit(2);
      resumm(nplot).pse = fit.fit(1);
      resumm(nplot).data = fit.data;
      resumm(nplot).xo = mean([sddata.xo]);
      resumm(nplot).exit = fit.exit;
      trialcounts(nplot) = sum([resumm(nplot).data.n]);
      nplot = nplot+1;
  else
        fprintf('%d Bad fit for %.2f\n',res,sd);
    end
else
        fprintf('No data for %.2f\n',sd);
 end
  
end
  fit = fitgauss([resumm.delay],1./[resumm.sd],state);
  if(nresample < 100)
      plot([resumm.delay],[resumm.sd],'o');
      hold on;
      plot(x,1./gauss([fit.mean fit.sd fit.amp],x),'g');
      drawnow;
  end
  if mod(res,100) == 0
    fprintf('Resample %d\n',res);
    drawnow;
  end
  counts(res) = sum(trialcounts);
  means(res) = fit.mean;
  sums{res} = resumm;
  exits(res) = fit.exit;
end
