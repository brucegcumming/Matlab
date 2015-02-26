function PlotExptFit(fit, varargin)
%Plots result from FitExpt
color = [];
matchstr = [];
showlabel = '';

if iscell(fit)
    for j = 1:length(fit)
        PlotExptFit(fit{j}, varargin{:});
    end
    return;
end

if length(fit) > 1
    colors = mycolors;
    for j = 1:length(fit)
        PlotExptFit(fit(j), 'color', colors{j},varargin{:});
    end
    return;
end

j = 1;
for j = 1:length(varargin)
    if strncmpi(varargin{j},'color',5)
        j = j+1;
        color = varargin{j};
    elseif strncmpi(varargin{j},'match',5)
        j = j+1;
        matchstr = varargin{j};
    elseif strncmpi(varargin{j},'label',5)
        showlabel = 'name';
    end
    j = j+1;
end

if ~isempty(matchstr)
    if isfield(fit,'expt') & regexp(fit.expt,'matchstr')
        go = 1;
    else
        return;
    end
end
       
h(1) = plot(fit.xv,fit.fitcurve);
if isfield(fit,'resp')
    hold on;
    if isfield(fit.state,'sds')
        if isfield(fit.state,'nreps')
            sem = fit.state.sds./sqrt(fit.state.nreps);
        else
            sem = fit.state.sds;
        end
        h(2) = errorbar(fit.x, fit.resp,sem,'o');
    else
        h(2) = plot(fit.x, fit.resp,'o');
    end
end
if ~isempty(color)
    set(h,'color',color);
end
if ~isempty(showlabel) && isfield(fit,'name')
    [a,b] = max(fit.fitcurve);
    [c,filename] = fileparts(fit.name);
    text(fit.xv(b),a,sprintf('%s',filename));
end