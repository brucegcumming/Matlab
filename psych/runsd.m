colors = {'r','b','g','c','y','m','k',[0 1 0.5],[0 0.5 1], [0.5 0 1], [0.5 1 0]};
prefix = '\\lsr-bgc1';


leftname = 'rufleft.mat';
rightname = 'rufright.mat';
prop = [];
alldata = [];
nplot = 1;
h = [];
labels = {};

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul04'],...
    'last',0,'first',1,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Jul04 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul08'],...
 'last',0,'first',1,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Jul08 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul09'],...
    'last',0,'first',1,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Jul09 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
title('Delay 0 ms');
drawnow;

figure
nplot = 1;
h = [];
labels = {};
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug19'],...
    'last',978,'first',1,'color',colors{nplot});
labels{nplot} = sprintf('Aug19 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug20'],...
    'last',1424,'first',1,'color',colors{nplot});
labels{nplot} = sprintf('Aug20 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug21'],...
    'last',943,'first',24,'color',colors{nplot});
labels{nplot} = sprintf('Aug21 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug22'],...
    'last',1344,'first',571,'color',colors{nplot});
labels{nplot} = sprintf('Aug9 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug22'],...
    'last',2243,'first',1937,'color',colors{nplot});
labels{nplot} = sprintf('Aug9 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
title('Delay 28 ms');

drawnow;

%SD 14
figure;
hold off;
nplot = 1;
h = [];
labels = {};

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul22'],...
    'last',1752,'first',1,'color',colors{nplot},'xmin',-0.1,'xmax',0.06,'sd',-1,'nmin',10);
labels{nplot} = sprintf('Jul22 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul18_1'],...
    'last',1752,'first',1,'color',colors{nplot},'xmin',-0.1,'xmax',0.1,'sd',-1);
labels{nplot} = sprintf('Jul18 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul29'],...
    'last',1752,'first',1,'color',colors{nplot},'xmin',-0.1,'xmax',0.1);
labels{nplot} = sprintf('Jul29 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;


[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Jul30_1'],...
    'last',143,'first',1,'color',colors{nplot},'xmin',-0.1,'xmax',0.1);
labels{nplot} = sprintf('Jul30 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
title('Delay 14 ms');
drawnow;


%everyting up to here has been 14ms frames
for j = 1:length(alldata)
alldata(j).sd = alldata(j).sd .* 1.4;
end

figure
hold off;
nplot = 1;
h = [];
labels = {};

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug06'],...
    'last',0,'first',100,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug06 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug07'],...
    'last',806,'first',20,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug07 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug08_10'],...
    'last',601,'first',1,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug08 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug11'],...
    'last',0,'first',1215,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug11 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
title('Delay 10 ms');
drawnow;


h = [];
labels = {};

figure
title('Delay 20 ms');
hold off;
nplot = 1;
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug12'],'last',999,'color','r','xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug12 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

%this one (Aug13) seems to break....
all = 0;
if all
[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug13'],'last',2005,'color','b','xmin',-0.05,'xmax',0.05);
alldata = [alldata Data];
labels{nplot} = sprintf('Aug13 %s',labels{nplot});
nplot = nplot+1;
  legend(h,labels,2);
drawnow;
end

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug14'],'last',1913,'color','g','xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug14 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug15'],...
    'last',1604,'first',1102,'color',colors{nplot},'xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug15 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

[Data, labels{nplot} h(nplot)] = readandplotpsych([prefix '/bgc/data/rufus/psych/Aug18'],'last',1789,'color','c','xmin',-0.05,'xmax',0.05);
labels{nplot} = sprintf('Aug18 %s',labels{nplot});
alldata = [alldata Data];
nplot = nplot+1;
  legend(h,labels,2);
drawnow;

for j = 1:length(alldata)
alldata(j).sd = alldata(j).sd .* 10;
end

ranges = [];
ranges(1).sd = 14;
ranges(1).min = -0.02;
ranges(1).max = 0.02;
ranges(1).nmin = 70;

side = 'L';
ylim = [0.005 0.2];
summarize;
suml = summ;

save(leftname,'suml');

return;


%Now combine across days.
Data = [];
ndat = 1;

j = 0;
for sd = unique([alldata.sd])
  idx = find([alldata.sd] == sd);
  j = j+1;
  for val = unique([alldata.x])
    id = find([alldata.x] == val & [alldata.sd] == sd);
    if ~isempty(id)
      Data(ndat).x = val;
      Data(ndat).expno = j;
      Data(ndat).n = sum([alldata(id).n]);
      Data(ndat).resp = sum([alldata(id).resp]);
      Data(ndat).p = Data(ndat).resp/Data(ndat).n;
      Data(ndat).name = alldata(ndat).name;
      Data(ndat).sd = sd;
      ndat = ndat+1;
    end
  end
end

figure
hold off;
nplot = 1;
h = [];
labels = {};
for sd = unique([Data.sd])
  sddata = Data(find([Data.sd] == sd));
  if(abs(sd) > 19)
    xr = 0.1;
  else
    xr = 0.03;
  end
  if length(sddata) >1
    fit = fitpsf(sddata,'xmin',-xr,'xmax',xr,'nmin',50,'sderr');
    h(nplot) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color',colors{nplot});
    labels{nplot} = sprintf('Sd = %.0f %.3g',sd,fit.fit(2));
    summ(nplot).delay = sd;
    summ(nplot).sd = fit.fit(2);
    summ(nplot).sdlim = fit.sdlim;
    summ(nplot).pse = fit.fit(1);
    nplot = nplot+1;
    hold on;
  end
end
  legend(h,labels,2);


return;
