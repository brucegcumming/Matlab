%Now combine across days.

suml = [];
sumr = [];
sumData = [];
ndat = 1;

j = 0;
for sd = unique([alldata.sd])
  idx = find([alldata.sd] == sd);
  j = j+1;
  for val = unique([alldata.x])
    id = find([alldata.x] == val & [alldata.sd] == sd);
    if ~isempty(id)
      sumData(ndat).x = val;
      sumData(ndat).expno = j;
      sumData(ndat).n = sum([alldata(id).n]);
      sumData(ndat).resp = sum([alldata(id).resp]);
      sumData(ndat).p = sumData(ndat).resp/sumData(ndat).n;
      sumData(ndat).name = alldata(ndat).name;
      sumData(ndat).sd = sd;
      sumData(ndat).xo = mean([alldata(id).xo]);
      ndat = ndat+1;
    end
  end
end

    
summ = [];

hold off;
nplot = 1;
h = [];
labels = {};
for sd = unique([sumData.sd])
  sddata = sumData(find([sumData.sd] == sd));
  if(abs(sd) > 19)
    xl = -1;
    xu = 1;
  else
    xl = -0.1;
    xu = 0.1;
  end
  nmin = setnmin;
  id = find([ranges.sd] == sd);
    if ~isempty(id)
      xl = ranges(id).min;
      xu = ranges(id).max;
      nmin = ranges(id).nmin;
    end
  if length(sddata) >1
    fit = fitpsf(sddata,'xmin',xl,'xmax',xu,'nmin',nmin,'sderr');
    if ~isnan(fit.fit(1))
    h(nplot) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color',colors{nplot});
    labels{nplot} = sprintf('Sd = %.0f %d Trials %.3g',sd, ...
			    sum([fit.data.n]),fit.fit(2));
    drawnow;
    summ(nplot).delay = sd;
    summ(nplot).sd = fit.fit(2);
    summ(nplot).sdlim = fit.sdlim;
    summ(nplot).pse = fit.fit(1);
    summ(nplot).data = fit.data;
    summ(nplot).xo = mean([sddata.xo]);
    nplot = nplot+1;
    else
      fprintf('Error in %d\n',mean([sddata.sd]));
    end;
    hold on;
  end
end
  legend(h,labels,2);

  
figure
h = [];
plotargs{1} = 'linewidth';
plotargs{2} = 2;

hold off;
labels = {};
lims = reshape([summ.sdlim],2,length([summ.sdlim])/2);
if(side == '1')
  plot([summ.delay], [summ.sd],'go');
elseif(side == 'L')
  suml = summ;
  h(1) = plot([summ.delay], [summ.sd],'go');
  hold on;
  e = errorbar([summ.delay], [summ.sd],lims(1,:)-[summ.sd], ...
	       [summ.sd]-lims(2,:),'go');
  labels{1} = sprintf('xo = %.2f',mean([summ.xo]));
  if exist(rightname,'file')
  load(rightname);
  lims = reshape([sumr.sdlim],2,length([sumr.sdlim])/2);  
  h(2) = plot([sumr.delay], [sumr.sd],'ro');
  errorbar([sumr.delay], [sumr.sd],lims(1,:)-[sumr.sd],[sumr.sd]- ...
	   lims(2,:),'ro');
  labels{2} = sprintf('xo = %.2f',mean([sumr.xo]));
  end
  legend(h,labels);
elseif(side == 'R')
  sumr = summ;
  labels{1} = sprintf('xo = %.2f',mean([summ.xo]));
  h(1) = plot([summ.delay], [summ.sd],'ro');
  hold on;
  if exist(leftname,'file')
    load(leftname);
    h(2) = plot([suml.delay], [suml.sd],'go');
    labels{2} = sprintf('xo = %.2f',mean([suml.xo]));
  end
  legend(h,labels);
else
  if exist(rightname,'file')
  load(rightname);
  plot([sumr.delay], [sumr.sd],'ro','MarkerFaceColor','r');
  hold on;
  end
  if exist(leftname,'file')
  load(leftname);
  plot([suml.delay], [suml.sd],'go','MarkerFaceColor','g');
  end
end
set(gca,'Yscale','log');

h = [];
labels = {};
%
%
% showgauss == 1 is broken. need showgauss = 2 to get through this
% (which is preferred anyway

if(showgauss == 1)
  state.freebase = 0;
  hold off;
  plot([suml.delay],1./[suml.sd],'go');
  hold on;
  lims = reshape([suml.sdlim],2,length([suml.sdlim])/2);  
  errs = [];
  errs(1,:) = 1./[suml.sd] - 1./lims(1,:);
  errs(2,:) = 1./lims(2,:) - 1./[suml.sd];
  errorbar([suml.delay], 1./[suml.sd],errs(2,:),errs(1,:),'go');
  plot([suml.delay], [suml.sd],'go','MarkerFaceColor','g');
  x = -30:30;
  fitl = fitgauss([suml.delay],1./[suml.sd],state);
  plot(x,gauss([fitl.mean fitl.sd fitl.amp],x),'g-');
  plot([sumr.delay],1./[sumr.sd],'ro');

  lims = reshape([sumr.sdlim],2,length([sumr.sdlim])/2);  
  errs = [];
  errs(1,:) = 1./[sumr.sd] - 1./lims(1,:);
  errs(2,:) = 1./lims(2,:) - 1./[sumr.sd];
  errorbar([sumr.delay], 1./[sumr.sd],errs(2,:),errs(1,:),'ro');

  fitr = fitgauss([sumr.delay],1./[sumr.sd],state);
  plot(x,gauss([fitr.mean fitr.sd fitr.amp],x),'ro');
elseif(showgauss == 2)
  state.freebase = 0;
  hold off;
  plot([suml.delay],[suml.sd],'go');
  hold on;
  lims = reshape([suml.sdlim],2,length([suml.sdlim])/2);  
  errs = [];
  errs(1,:) = [suml.sd] - lims(1,:);
  errs(2,:) = lims(2,:) - [suml.sd];
  e = errorbar([suml.delay], [suml.sd],errs(2,:),errs(1,:),'go');
  set(e,plotargs{:});
  hold on;
  h(1) = plot([suml.delay], [suml.sd],'go','MarkerFaceColor','g');
  labels{1} = sprintf('xo = %.2f',mean([suml.xo]));
  x = -30:30;
  fitl = fitgauss([suml.delay],1./[suml.sd],state);

  plot(x,1./gauss([fitl.mean fitl.sd fitl.amp],x),'g',plotargs{:});
  plot([sumr.delay],[sumr.sd],'ro',plotargs{:});

  lims = reshape([sumr.sdlim],2,length([sumr.sdlim])/2);  
  errs = [];
  errs(1,:) = [sumr.sd] - lims(1,:);
  errs(2,:) = lims(2,:) - [sumr.sd];
  e = errorbar([sumr.delay], [sumr.sd],errs(2,:),errs(1,:),'ro'); 
  set(e,plotargs{:});
  h(2) = plot([sumr.delay], [sumr.sd],'ro','MarkerFaceColor','r',plotargs{:});
  labels{2} = sprintf('xo = %.2f',mean([sumr.xo]));

  fitr = fitgauss([sumr.delay],1./[sumr.sd],state);
  plot(x,1./gauss([fitr.mean fitr.sd fitr.amp],x),'r',plotargs{:});
  set(gca,'Yscale','log');
end
if exist('ylim','var');
  set(gca,'Ylim',ylim);
  set(gca,'YTickLabel',{'0.01', '0.1'});
  set(gca,'TickLength',[0.025 0.05],'Color','none');
end

legend(h,labels);
set(gca,'linewidth',2,'FontSize',20);

