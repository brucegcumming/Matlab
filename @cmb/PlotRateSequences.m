function PlotRateSequences(DATA, varargin)

colors = mycolors;
offset = 0;
normalize = 0;

j = 1;
while j <= length(varargin) 
if strncmpi(varargin{j},'normalize',5)
normalize = 1;
elseif strncmpi(varargin{j},'offset',5)
offset = 1;
end
j = j+1;
end
hold off; 
AllTrials = [];
AllBlocks = [];
AllIds = [];
for j = 1:length(DATA.AllExpts)
E = DATA.AllExpts{j};
AllTrials = cat(2,AllTrials,[E.Trials.Trial]);
AllIds = cat(2,AllIds,[E.Trials.id]);
AllBlocks = cat(2,[E.Header.BlockStart], AllBlocks);
end
[AllTrials, id] = unique(AllTrials);
AllIds = AllIds(id);
AllBlocks = unique(AllBlocks);
dy = 0;
for j = 1:length(DATA.AllExpts)
E = DATA.AllExpts{j};
if normalize
scale = 1./mean([E.Trials.count]);
else
scale = 1;
end
h(j)  = plot(find(ismember(AllTrials,[E.Trials.Trial])), scale.*[E.Trials.count]+dy,'o','color',colors{j});
if normalize
dy = dy+offset;
else
dy = dy + mean([E.Trials.count]) .* offset;
end
labels{j} = sprintf('C%d',DATA.AllExpts{j}.Header.cellnumber);
hold on;
end
legend(h,labels);
yl = get(gca,'ylim');
for j = 2:length(AllBlocks)
[dt, t] = min(abs(AllTrials-AllBlocks(j)));
line([t t],yl,'linestyle','--');
text(t,yl(2),sprintf('Id%d',AllIds(t)),'rotation',90,'verticalalignment','top','horizontalalignment','right');
end


