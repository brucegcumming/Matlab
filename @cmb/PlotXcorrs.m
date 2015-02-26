function PlotXcorrs(DATA, plottype)

if isfield(DATA,'xcorrs') 
xc = DATA.xcorrs;
else
name = get(DATA.toplevel,'Name');
set(DATA.toplevel,'Name','Building Xcorrs')
xc = ExptListCorrs(DATA.AllExpts);
DATA.xcorrs = xc;
set(DATA.toplevel,'UserData',DATA,'Nmae',name);
end
hold off;
for j = 1:length(xc)
Ea = DATA.AllExpts{xc(j).cells(1)};
Eb = DATA.AllExpts{xc(j).cells(2)};
xcs(j) = xc(j).rsc(1);
if strcmp(plottype,'xcorr-distance')
plot(abs(xc(j).probesep),xc(j).rsc(1),'o','buttondownfcn',{@cmb.HitXcorr, j});
elseif strcmp(plottype,'xcorr-cp')
plot((0.5 - Eb.cp)*(0.5-Ea.cp),xc(j).rsc(1),'o','buttondownfcn',{@cmb.HitXcorr, j});
elseif strcmp(plottype,'xcorr-test')
plot(abs(Eb.cp-Ea.cp),xc(j).rsc(1),'o','buttondownfcn',{@cmb.HitXcorr, j});
cmb.PredictChoices(DATA);
end
hold on;
end
title(sprintf('Mean %.3f',mean(xcs(~isnan(xcs)))));
set(gcf,'UserData',DATA.toplevel);


