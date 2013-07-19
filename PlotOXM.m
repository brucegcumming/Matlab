function result = PlotOXM(Expt,state)

hold off;
colors = {[0 1 0], [1 0 0], [0 0 1],[1 1 0],[1 0 1]};
if(~exist('Expt') | ~isfield(Expt.Trials,'or'))
  return;
end

Expt = FillExpt(Expt, 'st');

[eye, eyes] = GetEval(Expt,'me','mode');

latency = 500;
duration = min([Expt.Trials.End] - [Expt.Trials.Start]);
for j = 1:length(Expt.Trials);
    counts(j) = length(find([Expt.Trials(j).Spikes] > latency & ...
			    [Expt.Trials(j).Spikes] < duration+latency));
end
ors = sort(unique([Expt.Trials.or]));

ie = 1;
h = [];
x = ors(1):(ors(end)-ors(1))/100:ors(end);
for eye = eyes;
  io = 1;
  for or = ors;
    idx = find([Expt.Trials.or] == or & [Expt.Trials.me] == eye & ...
	       [Expt.Trials.st] > 0);
    if ~isempty(idx)
      allcounts{io,ie} = counts(idx);
      result.ors(io,ie) = or;
      result.means(io,ie) = mean(counts(idx));
      result.n(io,ie) = length(counts(idx));
      result.sd(io,ie) = std(counts(idx));
      result.flags(io,ie) = 0;
      io = io+1;
    end
  end
  h(ie) = plot(result.ors(:,ie),result.means(:,ie),'o', ...
	       'color',colors{ie},'MarkerFaceColor',colors{ie});
  hold on;
  er = errorbar(result.ors(:,ie),result.means(:,ie),result.sd(:,ie) ./ sqrt(result.n(:,ie)),'o');
  set(er,'color',colors{ie});
  labels{ie} = sprintf('%s',val2str(eye,'me'));
  fitg = FitGauss(result.ors(:,ie),result.means(:,ie),0);
  plot(x,Gauss(x,fitg),'-','color',colors{ie},'LineWidth',2);
  ie = ie+1;
end
set(gca,'FontName','Arial','FontWeight','Bold','LineWidth',2, ...
	'Color','none');
if isfield(state,'fontsize')
  set(gca,'FontSize',state.fontsize);
end
xlabel('Orientation');
ylabel('Spikes/sec');
legend(h,labels);


function str = val2str(val, type)
if(strmatch(type,'me'))
  if(val == -1)
    str = 'Left';
  elseif(val == 0)
    str = 'Binoc';
  else
    str = 'Right';
  end
end

    
function Expt = FillExpt(Expt,field)

if ~isfield(Expt.Trials,field)
  val = eval(['Expt.Stimvals.' field]);
  eval(['[Expt.Trials.' field '] = deal(val)']);
end
