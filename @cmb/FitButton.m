function FitButton(a,b)

DATA = GetDataFromFig(a);
[DATA.Expt, plotres] = cmb.PlotCombined(DATA, DATA.Expt);
fit = FitExpt(plotres,'plotfit');
DATA = cmb.AddFitToData(DATA, plotres, fit);
set(DATA.toplevel,'UserData',DATA);

