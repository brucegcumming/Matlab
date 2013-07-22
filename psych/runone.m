hold off;
nplot = 1;
colors = mycolors;
labels = {};
h = [];

  
  


[Data, labels{nplot}, h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul22'],...
    'first',1,'color',colors{nplot},'xmin',-0.1,'xmax',0.1,'sd',-1);
labels{nplot} = sprintf('Jul22 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

