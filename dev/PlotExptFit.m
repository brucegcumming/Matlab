function PlotExptFit(fit, varargin)
%Plots result from FitExpt

if iscell(fit)
    for j = 1:length(fit)
        PlotExptFit(fit{j});
    end
    return;
end

plot(fit.xv,fit.fitcurve);
if isfield(fit,'resp')
    hold on;
    plot(fit.x, fit.resp,'o');
end
    