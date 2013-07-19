function [Fit,fval,exitflag] = FitExp(x,y, varargin)

%Fit = FitExp(x,y,....)  fits an exponential y = exp(-x/a)
%FitExp(x,Fit.params,'eval') evaluates and existing fit at x
%Fit = FitExp(x,y,'freebase')  fits an exponential y = c + exp(-x/a)
TAU= 1;
PEAK = 2;
BASE = 3;
state.freebase = 0;

guess = [];


j = 1;
while j <= nargin -2
  if strncmpi(varargin{j},'guess',5)
    j = j+1;
    guess = varargin{j};
  elseif strncmpi(varargin{j},'freebase',7)
      state.freebase = 1;
  elseif strncmpi(varargin{j},'eval',4)
    Fit = EvalExp(x,y);
    fval = 0;
    exitflag = 0;
    return;
  end
  j = j+1;
end

if isempty(guess)
    guess(PEAK) = mean(y(1:3));


    if(isfield(state,'freebase') & state.freebase)
        guess(BASE) = mean(y(end-10:end));
        guess(TAU) = trapz(x,y) - guess(BASE) * (x(end) - x(1));
    else
        guess(TAU) = trapz(x,y);
    end
end

if ~isfield(state,'logy')
  state.logy = 0;
end

if(isfield(state,'sdmax'))
  sdmax = state.sdmax;
else
  sdmax = 0;
end
if(isfield(state,'minoffset'))
  minoffset = state.minoffset;
else
  minoffset = NaN;
end
Fit.guessrss = Minimise(guess,x,y,0,NaN,state);
Fit.guess = guess;

options = optimset('MaxFunEvals',100000,'maxiter',5000);
[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,x,y,sdmax,minoffset,state);
Fit.peak = fittedparams(PEAK);
Fit.tau = fittedparams(TAU);
if(isfield(state,'freebase') & state.freebase)
  Fit.base = fittedparams(BASE);
else
  FitGauss.base = 0;
end
  Fit.params = fittedparams;
FitGauss.rss = fval;
FitGauss.exit = exitflag;


function fity = EvalExp(x, params)

if(length(params) == 3)
  fity = params(3) + params(2) .* exp(-x/params(1));
else
  fity = params(2) .* exp(-x/params(1));
end  

function SSD = Minimise(params,x,y,sdmax,minoffset,state)

%BASE=4 AMP=1, SD=2, XCEN=3,

fity = EvalExp(x, params);
if(state.logy)
  e = min(y(find(y > 0)))/1000000;
  diffs = (log(fity+e) - log(y+e)).^2;
else
  diffs = (fity - y).^2;
end

SSD = sum(diffs);
if(sdmax > 0 & abs(params(2)) > sdmax)
  SSD = NaN;
end
if(~isnan(minoffset) & params(3) < minoffset)
  SSD = NaN;
end
