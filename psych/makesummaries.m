alltrials = [];
alldata = [];
h = [];
labels = {};
colors = mycolors;
hold off;

for j = 1:alln
  prop(j).color = colors{nplot};
  [Data, Trials, lb, hn] = readplotpsych(prefix, prop(j).file, prop(j));
  if ~isempty(Data)
  labels{nplot} = sprintf('%s %s',prop(j).file,lb);
  h(nplot) = hn;
  for k = 1:length(Data)
    Data(k).sd = Data(k).sd * prop(j).sdk;
    Data(k).file = prop(j).file;
  end
  alldata = [alldata Data];
  len = length(alltrials);
% have to do this the clunky way because the order of the fields is
% not always the same....
  for tr = 1:length(Trials)
    alltrials(len+tr).xo = Trials(tr).xo;
    alltrials(len+tr).yo = Trials(tr).yo;
    alltrials(len+tr).sd = Trials(tr).sd * prop(j).sdk;
    alltrials(len+tr).sign = Trials(tr).sign;
    alltrials(len+tr).score = Trials(tr).score;
    alltrials(len+tr).x = Trials(tr).x;
  end
  
  legend(h, labels,2);
  nplot = nplot+1;
end

if ~isempty(prop(j).title)
    title(prop(j).title);
    nplot = 1;
    h = [];
    labels = {};
    figure;
    hold off;
  end
  drawnow;
end
