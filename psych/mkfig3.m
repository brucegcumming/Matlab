
idx = find([suml.delay] == -14);
plotpsych(suml(idx).data,suml(idx).pse,suml(idx).sd,'color','r');
hold on;
idx = find([suml.delay] == 14);
plotpsych(suml(idx).data,suml(idx).pse,suml(idx).sd,'color','b');
set(gca,'Xlim',[-0.1 0.1],'color','none');
set(gca,'YTickLabel',{},'XTickLabel',{},'TickLength',[0 0]);
