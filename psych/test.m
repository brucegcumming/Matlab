sums = {};
Trials = [];

for j = 1:length(alltrials)
    if isempty(alltrials(j).score)
        alltrials(j).score = 0;
    end
    if isempty(alltrials(j).sd)
        alltrials(j).sd = nan;
    end
end

j = 1;
ndat = 1;

for sd = unique([sumData.sd])
    
  idx = find([alltrials.sd] == sd & [alltrials.score] ~= 0);
  sdtrials{j} = alltrials(idx);
  for x = unique([sdtrials{j}.x])
      idx = find([sdtrials{j}.x] == x);
      tdata(ndat).ups = [sdtrials{j}(idx).score] .* [sdtrials{j}(idx).sign];
      tdata(ndat).x = x;
      tdata(ndat).sd = sd;
      tdata(ndat).expno = j;
       bredata(ndat).n = length([tdata(ndat).ups]);
        bredata(ndat).sd = sd;
        bredata(ndat).x = x;
        bredata(ndat).xo = mean([sdtrials{j}(idx).xo]);
       bredata(ndat).expno = j;
       ndat = ndat +1;
   end
  Trials = [Trials sdtrials{j}];
  j = j+1;
end
ndat = ndat -1;

x = -30:30;
for res = 1:nresample
  nplot = 1;
  state.freebase = 0;
  state.logy = 1;
  Trials = [];
  for j = 1:length(sdtrials);
    id = 1:length(sdtrials{j});
    rid = Bresample(id);
    retrials{j} = sdtrials{j}(rid);
    Trials = [Trials retrials{j}];
end
aredata = CountPsychTrials(Trials,'sortby','sd');

    for j = 1:ndat
        rid = Bresample(1:length(tdata(j).ups));
        ups = [tdata(j).ups];
        bredata(j).resp = length(find(ups(rid) >= 1));
        bredata(j).p = bredata(j).resp / bredata(j).n;
    end
redata = bredata;
trialcounts = [];
for sd = unique([sumData.sd])
  sddata = redata(find([redata.sd] == sd));
  if length(sddata) >1
    fit = fitpsf(sddata,'xmin',xl,'xmax',xu,'nmin',nmin);
    if ~isnan(fit.fit(1))
      summ(nplot).delay = sd;
      summ(nplot).sd = fit.fit(2);
      summ(nplot).pse = fit.fit(1);
      summ(nplot).data = fit.data;
      summ(nplot).xo = mean([sddata.xo]);
      trialcounts(nplot) = sum([summ(nplot).data.n]);
      nplot = nplot+1;
    end
  end
end
  fit = fitgauss([summ.delay],1./[summ.sd],state);
  if(nresample < 100)
      plot([summ.delay],[summ.sd],'o');
      hold on;
      plot(x,1./gauss([fit.mean fit.sd fit.amp],x),'g');
      drawnow;
  end
  if mod(res,10) == 0
    fprintf('Resample %d\n',res);
    drawnow;
  end
  counts(res) = sum(trialcounts);
  means(res) = fit.mean;
  sums{res} = summ;
end
