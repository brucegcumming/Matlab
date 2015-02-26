function DATA = PlotXcorr(DATA, probes, type)
nc = length(probes);
GetFigure('CrossCorrelation');
if nc > 2
for j = 1:nc
subplot(nc,nc,((j-1)*nc)+1);
for k = 1:j
subplot(nc,nc,((j-1)*nc) + k);
if type == 3
cmb.CalcXcorr(DATA, DATA.currentexpt(1),probes(j),probes(k));
elseif type == 1
cmb.CalcXcorrDC(DATA, DATA.currentexpt(1),probes(j),probes(k));
else
cmb.CalcXcorrV(DATA, DATA.currentexpt(1),probes(j),probes(k));
end
set(gca,'ytick',[],'xtick',[]);
if k == 1
ylabel(sprintf('P%d',probes(j)));
end
if j == nc
xlabel(sprintf('P%d',probes(k)));
end
end
end
else
if type == 3
cmb.CalcXcorr(DATA, DATA.currentexpt(1),probes(1),probes(2));
elseif type == 1
cmb.CalcXcorrDC(DATA, DATA.currentexpt(1),probes(1),probes(2));
else
cmb.CalcXcorrV(DATA, DATA.currentexpt(1),probes(1),probes(2));
end
end


