function PlotSDth(data, varargin)

idl = find([data.xo] < 0);
idr = find([data.xo] > 0);
suml = data(idl);
sumr = data(idr);

base = 0.0000001;
plotargs{1} = 'linewidth';
plotargs{2} = 2;
j = 1;
while j < nargin-2
  if strncmpi(varargin{j},'ylim')
    j = j+1;
    ylim = varargin{j};
  end
  j = 1+1;
end

state.freebase = 0;
state.logy = 1;
  plot([suml.delay],[suml.sd],'go');
  hold on;
  if isfield(suml,'sdlim');
    lims = reshape([suml.sdlim],2,length([suml.sdlim])/2);  
    errs(1,:) = [suml.sd] - lims(1,:);
    errs(2,:) = lims(2,:) - [suml.sd];
    if mean(errs(2,:)) < 0
      errs(2,:) = [suml.sd] - lims(2,:);
      errs(1,:) = lims(1,:) - [suml.sd];
    end
    e = errorbar([suml.delay], [suml.sd],errs(2,:),errs(1,:),'go');
    set(e,plotargs{:});
  end
  hold on;
  h(1) = plot([suml.delay], [suml.sd],'go','MarkerFaceColor','g');
  fitl = fitgauss([suml.delay],1./[suml.sd],state);
  labels{1} = sprintf('xo = %.2f, sd %.2f delay = %.2f',mean([suml.xo]),fitl.sd,fitl.mean);

  clear errs;
  step = (max([suml.delay sumr.delay])-min([suml.delay sumr.delay]))/100;;
  x = min([suml.delay sumr.delay]):step:max([suml.delay sumr.delay]);
  plot(x,1./(gauss([fitl.mean fitl.sd fitl.amp],x)+base),'g',plotargs{:});
  plot([sumr.delay],[sumr.sd],'ro',plotargs{:});

  if isfield(sumr,'sdlim');
  lims = reshape([sumr.sdlim],2,length([sumr.sdlim])/2);  
  errs(1,:) = [sumr.sd] - lims(1,:);
  errs(2,:) = lims(2,:) - [sumr.sd];
  if mean(errs(2,:)) < 0
  errs(2,:) = [sumr.sd] - lims(2,:);
  errs(1,:) = lims(1,:) - [sumr.sd];
  end
  
  e = errorbar([sumr.delay], [sumr.sd],errs(2,:),errs(1,:),'ro'); 
  set(e,plotargs{:});
  end
  h(2) = plot([sumr.delay], [sumr.sd],'ro','MarkerFaceColor','r',plotargs{:});

  fitr = fitgauss([sumr.delay],1./[sumr.sd],state);
  labels{2} = sprintf('xo = %.2f, sd %.2f delay = %.2f',mean([sumr.xo]),fitr.sd,fitr.mean);
  plot(x,1./gauss([fitr.mean fitr.sd fitr.amp],x),'r',plotargs{:});
  set(gca,'Yscale','log');

  
if exist('ylim','var');
  set(gca,'Ylim',ylim);
  set(gca,'YTickLabel',{'0.01', '0.1'});
  set(gca,'TickLength',[0.025 0.05],'Color','none');
  set(gca,'linewidth',2,'FontSize',20);
end
legend(h,labels);