function CalcTetrodeCorrs(DATA)


for j = 1:length(allspks);
    xc = corrcoef(reshape(DATA.AllData.Spikes.values(allspks(j),:),32,4));
    xcs(j,:) = [xc(1,2) xc(2,3) xc(3,4)];
end
