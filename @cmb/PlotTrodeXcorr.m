function PlotTrodeXcorr(a,b)
DATA = GetDataFromFig(a);
allspks = DATA.spklist;
for j = 1:length(allspks);
xc = corrcoef(reshape(DATA.AllData.Spikes.values(allspks(j),:),32,4));
xcs(j,:) = [xc(1,2) xc(2,3) xc(3,4)];
end
GetFigure('SpkieCorr');
hold off;
plot(xcs(:,2),xcs(:,3),'.','markersize',DATA.ptsize);
hold on; plot(median(xcs(:,2)),median(xcs(:,3)),'r+');


