function PlotPCAs(DATA, eid)
for e = 1:length(eid)
[DATA, allspks] = SetExptSpikes(DATA,eid(e),'setrange');
end

GetFigure('PCA');
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;];
pw=0.33;
ph = 0.5;
cid = unique(DATA.AllData.Spikes.codes(allspks,2));
for j = 1:size(pairs,1);
x = mod(j-1,3) * 0.33;
y = floor((j-1)/3) * 0.5;
subplot('Position' ,[x y pw ph]);
hold off;
for k = 1:length(cid)
id = find(DATA.AllData.Spikes.codes(allspks,2) == cid(k));
plot(DATA.AllData.pcs(allspks(id),pairs(j,1)),DATA.AllData.pcs(allspks(id),pairs(j,2)),'.','markersize',1,'color',DATA.spkcolor{k});
hold on;
end
set(gca,'Xtick',[],'YTick',[]);
end

