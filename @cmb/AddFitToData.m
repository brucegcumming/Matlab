function DATA = AddFitToData(DATA, plotres, fit)

type = strmatch(plotres.type{1},{'Op' 'Pp' 'or'});
h = get(gca,'title');
s = get(h,'string');
if ~isempty(type)
if type(1) == 1
DATA.fitvals.Op = fit.mean;
DATA.fitvals.Opw = fit.sd;
DATA.fitvals.OpRo = GetEval(DATA.Expt,'Ro');
DATA.Expt.fits.Op = fit;
set(h,'string',[s sprintf(' FitPeak = %.1f',fit.mean)]);
elseif type(1) == 2
DATA.fitvals.Pp = fit.mean;
DATA.fitvals.PpRo = GetEval(DATA.Expt,'Ro');
DATA.fitvals.Ppw = fit.sd;
DATA.Expt.fits.Pp = fit;
set(h,'string',[s sprintf(' FitPeak = %.1f',fit.mean)]);
elseif type(1) == 3
set(h,'string',[s sprintf(' FitPeak = %.1f',fit.mean)]);
end
DATA.Expt.fit = fit;
end


