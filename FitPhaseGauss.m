function [FitGauss,fval] = FitPhaseGauss(x,y, varargin)
%
% FitPhaseGuass(x,y) fits the sum of a Gaussian and a Cumulative Gaussian
% to resposnes. Used for fitting sqcorrug responses as a function of phase.

BASE = 4;
AMP = 1;
SD = 2;
MEAN = 3;
CUMAMP = 5;


if length(x) == length(y)
    guess(MEAN) = mean(x .* y)/mean(y);
    guess(SD) = std(x .* y)/mean(y);
    guess(SD) = 25;
    guess(MEAN) = 0;
    guess(CUMAMP) = 0;
end

state.freebase = 1;
state.guessmean = NaN;
state.guesssd = NaN;
state.sdmax = 0;
state.logy = 0;
state.minoffset = 0;
state.sdmin = mean(diff(x))/2;
state.sdmin = 10;
state.minu = -40;
state.maxu = 40;
state.maxy = max(y);

j = 1;
while j < nargin -1
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


guess(CUMAMP) = 0;
if(state.freebase)
  guess(BASE) = min(y);
  guess(AMP) = (max(y) - min(y))/1.2;
else
  guess(AMP) = max(y)/1.2;
end

if ~isnan(state.guessmean)
    guess(MEAN) = state.guessmean;
else
    guess(MEAN) = mean(x);
    guess(MEAN) = 20;
end
if ~isnan(state.guesssd)
    guess(MEAN) = state.guesssd;
end


 sdmax = state.sdmax;

  minoffset = state.minoffset;

FitGauss.guessrss = Minimise(guess,x,y,0,NaN,state);
FitGauss.guess = guess;

options = optimset('MaxFunEvals',100000,'maxiter',5000);
[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,x,y,sdmax,minoffset,state);
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



function SSD = Minimise(params,x,y,sdmax,minoffset,state)

%BASE=4 AMPGAUSS=1, SD=2, XCEN=3, AMPCUMGAUSS = 5; 


fity = params(4) + params(1) .* exp(-(x-params(3)).^2/(2 * ...
						  params(2)^2));


                      fity = fity + params(5) .* erf(((x-params(3))/(sqrt(2) * abs(params(2))) +1))/2;

if(state.logy)
  e = min1(y(find(y > 0)))/1000000;
  diffs = (log(fity+e) - log(y+e)).^2;
else
  diffs = (fity - y).^2;
end

SSD = sum(diffs);

if(sdmax > 0 & abs(params(2)) > sdmax)
  SSD = NaN;
end
if(state.sdmin > 0 & abs(params(2)) < state.sdmin)
  SSD = NaN;
end
if(~isnan(minoffset) & params(3) < minoffset)
  SSD = NaN;
end

if params(4) < 0 | params(4) > state.maxy
    SSD = NaN;
end

if params(3) > state.maxu | params(3) < state.minu
    SSD = NaN;
end

function fity = PhaseGauss(x, params)

fity = params(4) + params(1) .* exp(-(x-params(3)).^2/(2 * ...
						  params(2)^2));

fity = fity + params(5) .* erf(((x-params(3))/(sqrt(2) * params(2)) +1))/2;

