function [FitGauss,fval] = FitPhaseGaussAll(x,y, ipy, varargin)
%
% FitPhaseGuass(x,y) fits the sum of a Gaussian and a Cumulative Gaussian
% to resposnes. Used for fitting sqcorrug responses as a function of phase.

BASE = 1;
CUMAMP = 2;
AMP = 3;
SD = 4;
MEAN = 5;
CUMAMP = 6;



p(1:size(x,2)) = deal(-1000);
q(1:size(x,2)) = deal(1000);
px = [p;x; q];
p(1:size(x,2)) = deal(ipy(1));
q(1:size(x,2)) = deal(ipy(2));
py = [p;y; q];



%
% if lengths match, this is fitting so need to make a guess. Otherwise,
% just evalutating a fit, when y is fit parameters
if length(x) == length(y)
    for j = 1:size(px,2)
        guess(MEAN,j) = 0;
%        guess(j,SD) = std(x(j,:) .* y(j,:))/mean(y(j,:));
        guess(SD,j) = 20;
        guess(CUMAMP,j) = ipy(2) - ipy(1);
        guess(BASE,j) = ipy(1);
        guess(AMP,j) = range(y(:,j))/2;
    end
end



state.freebase = 1;
state.guessmean = NaN;
state.guesssd = NaN;
state.sdmax = 0;
state.logy = 0;
state.minoffset = NaN;
state.sdmin = mean(mean(diff(x))/2);


j = 1;
while j < nargin -2
    if strncmpi(varargin{j},'guess',5)
        j = j+1;
        guess = varargin{j};
    elseif strncmpi(varargin{j},'state',5)
        j = j+1;
        state = varargin{j};
    elseif strncmpi(varargin{j},'eval',4)
        j = j+1;
        FitGauss = PhaseGauss(x,y);
        fval = 0;
        return;
    end
  j = j+1;
end


if(state.freebase)
  guess(AMP,:) = (max(py) -min(py))./1.2;
  guess(BASE,:) = deal(min(ipy));
else
  guess(AMP,:) = max(y)/1.2;
end

if ~isnan(state.guessmean)
    guess(MEAN,:) = state.guessmean;
end
if ~isnan(state.guesssd)
    guess(MEAN,:) = state.guesssd;
end


 sdmax = state.sdmax;

  minoffset = state.minoffset;

FitGauss.guessrss = Minimise(guess,px,py,0,NaN,state);
FitGauss.guess = guess;

options = optimset('MaxFunEvals',100000,'maxiter',5000);

for j = 1:size(px,2)
    fit = FitPhaseGauss(px(:,j),py(:,j));
    params(1,j) = fit.params(4);
    params(2,j) = fit.params(5);
    params(3,j) = fit.params(1);
    params(4,j) = fit.params(2);
    params(5,j) = fit.params(3);
end
FitGauss.params = params;
return;

[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,px,py,sdmax,minoffset,state);
FitGauss.amp = fittedparams(AMP);
FitGauss.sd = fittedparams(SD);
if state.freebase
  FitGauss.base = fittedparams(BASE);
else
  FitGauss.base = 0;
end
FitGauss.mean = fittedparams(MEAN);
FitGauss.rss = fval;
FitGauss.exit = exitflag;
FitGauss.cumamp = fittedparams(CUMAMP);
FitGauss.params = fittedparams;
FitGauss.x = px;
FitGauss.y = y;


function SSD = Minimise(params,x,y,sdmax,minoffset,state)

%BASE=4 AMPGAUSS=1, SD=2, XCEN=3, AMPCUMGAUSS = 5; 
%BASE = 1;
%CUMAMP = 2;
%AMP = 3;
%SD = 4;
%MEAN = 5;

params(4,:) = abs(params(4,:));

for j = 1:size(x,2)
    if(state.sdmin > 0 & abs(params(4,j)) < state.sdmin)
        params(4,j) = state.sdmin;
    end
    fity(:,j) = params(1,j) + params(3,j) .* exp(-(x(:,j)-params(5,j)).^2/(2 * ...
        params(4,j)^2));
    fity(:,j) = fity(:,j) + params(2,j) .* erf(((x(:,j)-params(5,j))/(sqrt(2) * params(4,j))) +1)/2;
end

if(state.logy)
  e = min(y(find(y > 0)))/1000000;
  diffs = (log(fity+e) - log(y+e)).^2;
else
  diffs = (fity - y).^2;
end

SSD = sum(sum(diffs));

if(sdmax > 0 & abs(params(2)) > sdmax)
  SSD = NaN;
end
if(~isnan(minoffset) & params(3) < minoffset)
  SSD = NaN;
end


function fity = PhaseGauss(x, params)

params(4,:) = abs(params(4,:));
for j = 1:size(x,2)
    fity(:,j) = params(1,j) + params(3,j) .* exp(-(x(:,j)-params(5,j)).^2/(2 * ...
        params(4,j)^2));
    fity(:,j) = fity(:,j) + params(2,j) .* erf(((x(:,j)-params(5,j))/(sqrt(2) * params(4,j))) +1)/2;
end
