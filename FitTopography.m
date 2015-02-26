function [fit, fity] = FitTopograpy(x,y,px,py, varargin)
%FitTopograpy(x,y,px,py) Fit RF topography for penetraions


fittype = 'plane';
guess = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'cubic',4)
        fittype = 'cubic';
    elseif strncmpi(varargin{j},'guess',4)
        j = j+1;
        guess = varargin{j};
    elseif strncmpi(varargin{j},'planar',4)
        fittype = 'plane';
    elseif strncmpi(varargin{j},'quad',4)
        fittype = 'quadratic';
    elseif strncmpi(varargin{j},'rotate',4)
        fittype = 'rotate';
    end
    j = j+1;
end
if isstruct(x) %evaluate a fit
    if nargin > 3 && strcmp(py,'invert')
        [fit, fity] = InvertFit(x,y,px);
    else
        [fit, fity] = EvaluateFit(x,y,px);
    end
    return;
end

id = find(~isnan(x) & ~isnan(y) & ~isnan(px) & ~isnan(py));
Map.x = x(id);
Map.px = px(id);
Map.py = py(id);
Map.y = y(id);

params(1) = 0; %rotation;
params(2) = 1; %scaling mm -> degrees
params(3) = 0; %offset degree VA at0
angles = 0:0.2:2*pi;
for j = 1:length(angles)
    [a,b,xcs(j)] =TryMap([angles(j) params(2:end)], Map);
end
[a,b] = max(xcs);
params(1) = angles(b);
    
if strcmp(fittype,'rotate')
options = optimset('MaxFunEvals',100000,'display','off');
[fittedparams,fval,exitflag, output] = fminsearch(@TryMap,params, options,Map);
[ssq(1), fitted] = TryMap(fittedparams, Map);
fit.xfit = fittedparams;
fit.xfitted = fitted;
angles = fitted(1) + [0 pi];
clear xcs;
Map.x = y(id);
for j = 1:length(angles)
    [a,b,xcs(j)] =TryMap([angles(j) params(2:3)], Map);
end
[a,b] = max(xcs);
params(1) = angles(b);
params(3) = 0;
guess = params;
[fittedparams,fval,exitflag, output] = fminsearch(@TryMap,params, options,Map);
[ssq(2), fitted] = TryMap(fittedparams, Map);

fit.yfit = fittedparams;
fit.yfitted = fitted;
else
    if ~isempty(guess)
        xguess = guess(1,:);
        yguess = guess(2,:);
    elseif strcmp(fittype,'quadratic')
        xguess = [5 0.8 0.8 0 0];
        yguess = [5 0.8 0.8 0 0];
    elseif strcmp(fittype,'cubic')
        xguess = [5 0.8 0.8 0 0 0 0];
        yguess = [5 0.8 0.8 0 0 0 0];
    else
        xguess = [2 0.5 0.5];
        yguess = [5 0.8 0.8 ];
    end
      
    [guessq(1), fitz] = TryPlane(xguess,Map.px,Map.py,Map.x);
    options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
    fit.xfit = fminsearch(@TryPlane, xguess, options, Map.px,Map.py,Map.x);
    [ssq(1), fit.xfitted] = TryPlane(fit.xfit,Map.px,Map.py,Map.x);
    [guessq(1), fitz] = TryPlane(yguess,Map.px,Map.py,Map.y);
    fit.yfit = fminsearch(@TryPlane, yguess, options, Map.px,Map.py,Map.y);
    [ssq(2), fit.yfitted] = TryPlane(fit.yfit,Map.px,Map.py,Map.y);
end
fit.fittype = fittype;
fit.px = Map.px;
fit.py = Map.py;
fit.guess = guess;
fit.ssq = ssq;

function [ssq, fitx, xc] = TryMap(x, Map)

cosa = cos(x(1));
sina = sin(x(1));
fitx = x(3) + (x(2) .* (Map.px .* cosa + Map.py .* sina));
xc = corrcoef(Map.x, fitx);
xc = xc(1,2);
ssq = sum((fitx-Map.x) .^2);

function [ssq, p] = TryPlane(x, X,Y,Z)

if length(x) ==7
    p = x(1) + x(2) .* X + x(3) .* Y + x(4) .* X.^2 + x(5) .* Y.^2 + x(6) .* X.^3 + x(7) .* Y.^3;
elseif length(x) ==5
    p = x(1) + x(2) .* X + x(3) .* Y + x(4) .* X.^2 + x(5) .* Y.^2;
else
    p = x(1) + x(2) .* X + x(3) .* Y;
end
if nargin > 3
    id = find(~isnan(Z));
    ssq = sum((Z(id) - p(id)).^2);
else
    ssq = NaN;
end

function [fitx, fity] = InvertFit(fit, x,y, varargin)
%z,y are positions in space.  Fitx fitx are penetration co-ordninates
fx = linspace(min(fit.px(:)),max(fit.px(:)));
fy = linspace(min(fit.py(:)),max(fit.py(:)));
[X,Y] = meshgrid(fx,fy);
[xfit, yfit] = EvaluateFit(fit,X,Y);
for j = 1:length(x(:))
    [d(j), id(j)] = min(abs((x(j)+i.*y(j)) - (xfit(:) + i.* yfit(:))));
end
fitx = X(id);
fity = Y(id);


function [fx, fy] = EvaluateFit(fit, x,y, varargin)

if strcmp(fit.fittype,'rotate')
    cosa = cos(fit.xfit(1));
    sina = sin(fit.xfit(1));
    fx = fit.xfit(3) + (fit.xfit(2) .*(x.*cosa + y.* sina));
    cosa = cos(fit.yfit(1));
    sina = sin(fit.yfit(1));
    fy = fit.yfit(3) + (fit.yfit(2) .*(x.*cosa + y.* sina));
else
    [~, fx] = TryPlane(fit.xfit,x,y);
    [~, fy] = TryPlane(fit.yfit,x,y);
end
fitted = cat(3,fx,fy);

