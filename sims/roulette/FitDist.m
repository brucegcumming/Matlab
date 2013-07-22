function [FitGauss,fval,exitflag] = FitDist(x,y, type, varargin)
%[FitGauss,fval,exitflag] = FitGauss(x,y, type, varargin)
% type = 0 = one Gaussians
% type = 1 = sum of two Gaussians
%        2 Gauss + periodic Gauss
% type = 3 Gauss + Skewed Gauss
% type = 4 = skewed Gaussian

MEAN = 1;
SD = 2;
AMP = 3;

MEAN2 = 4;
SD2=5;
AMP2 = 6;
SKEW2 = 7;
SKEW = 4;
period = pi/2;


autocenter = 1;  %% adjust parameters so that peak is at X = 0;

guessed = 0;
if length(y) < length(x)
    params = y;
    y = zeros(size(x));
end
state.logy = 0;

guess(MEAN) = mean(x .* y)/mean(y);
guess(SD) = std(x .* y)/mean(y);
guess(AMP) = max(y)/1.2;
if type == 1
    guess(MEAN2) = guess(MEAN);
    guess(AMP2) = guess(AMP)/2;
    guess(AMP) = guess(AMP)/2;
    guess(SD2) = guess(SD)/2;
elseif type == 4
    guess(SKEW) = 1.1;
elseif type == 2
    guess(MEAN2) = guess(MEAN);
    guess(AMP2) = guess(AMP)/4;
    guess(SD2) = pi/16;
elseif type == 3
    guess(MEAN2) = guess(MEAN);
    guess(AMP2) = guess(AMP)/2;
    guess(AMP) = guess(AMP)/2;
    guess(SD2) = guess(SD)/2;
    guess(SKEW2) = 1.0;
end

state.freebase = 0;

j = 1;
while j <= nargin -3
  if strncmpi(varargin{j},'freebase',5)
      guess(BASE) = min(y);
      guess(AMP) = (max(y) - min(y))/1.2;
      state.freebase = 1;
  elseif strncmpi(varargin{j},'eval',4)
      [a, FitGauss] = Minimise(params,x,y,type,NaN,state);
      return;
  elseif strncmpi(varargin{j},'posbase',5)
      guess(BASE) = min(abs(y));
      guess(AMP) = (max(y) - min(abs(y)))/1.2;
      state.freebase = 2;
  elseif strncmpi(varargin{j},'mean',7)
    j = j+1;
    guess(MEAN) = varargin{j};
  elseif strncmpi(varargin{j},'sd',5)
    j = j+1;
    guess(SD) = varargin{j};
  elseif strncmpi(varargin{j},'guess',5)
    j = j+1;
    guess = varargin{j};
    guessed = 1;
  end
  j = j+1;
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

[FitGauss.guessrss, fitted] = Minimise(guess,x,y,type,NaN,state);

if guessed
    FitGauss.guessfit = fitted;
end
FitGauss.guess = guess;

options = optimset('MaxFunEvals',100000,'maxiter',5000);
[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,x,y,type,minoffset,state);
FitGauss.amp = fittedparams(AMP);
FitGauss.sd = fittedparams(SD);
if(isfield(state,'freebase') & state.freebase)
    if state.freebase == 2 %% baseline > 0
        fittedparams(BASE) = abs(fittedparams(BASE));
    end
        FitGauss.base = fittedparams(BASE);
else
  FitGauss.base = 0;
end
FitGauss.mean = fittedparams(MEAN);
FitGauss.rss = fval;
FitGauss.exit = exitflag;
[a, FitGauss.fitted] = Minimise(fittedparams,x,y,type,NaN,state);
FitGauss.params = fittedparams; %% in order needed for Gauss(..)

if autocenter
    [a,b] = max(FitGauss.fitted);
    if ismember(type,[0 4])
        FitGauss.params(1) = FitGauss.params(1) - x(b);
    elseif ismember(type,[1 3]) 
        FitGauss.params(1) = FitGauss.params(1) - x(b);
        FitGauss.params(4) = FitGauss.params(4) - x(b);
    elseif type == 2 %%+periodic Gaussian
        FitGauss.params(1) = FitGauss.params(1) - x(b);
        FitGauss.params(MEAN2) = mod(FitGauss.params(MEAN2) - x(b),period);
    end
end
FitGauss.type = type;


function [SSD, fity] = Minimise(params,x,y,type,minoffset,state)

%BASE=4 AMP=1, SD=2, XCEN=3,
MEAN = 1;
SD = 2;
AMP = 3;

MEAN2 = 4;
SD2=5;
AMP2 = 6;
SKEW2 = 7;

SKEW = 4;
period = pi/2;
if type == 0 
fity = params(AMP) .* exp(-(x-params(MEAN)).^2/(2 * ...
						  params(SD)^2));
elseif type == 1 
fita = params(AMP) .* exp(-(x-params(MEAN)).^2/(2 * ...
						  params(SD)^2));
fitb = params(AMP2) .* exp(-(x-params(MEAN2)).^2/(2 * ...
						  params(SD2)^2));
                      fity = fita+fitb;
elseif type == 4 
    xs = (x - min(x))./range(x);
    xs = min(x) + xs.^params(SKEW) .*range(x);
    fity = params(AMP) .* exp(-(xs-params(MEAN)).^2/(2 * ...
						  params(SD)^2));
elseif type == 2 %Gauss + periodic gauss;
    xs = mod(x-params(MEAN2),period)-period/2;
fita = params(AMP) .* exp(-(x-params(MEAN)).^2/(2 * ...
						  params(SD)^2));
fitb = fita .* params(AMP2) .* exp(-(xs).^2/(2 * ...
						  params(SD2)^2));
                      fity = fita+fitb;
elseif type == 3 
    xs = (x - min(x))./range(x);
    xs = min(x) + (xs.^params(SKEW2)) .*range(x);
fita = params(AMP) .* exp(-(xs-params(MEAN)).^2/(2 * ...
						  params(SD)^2));
fitb = params(AMP2) .* exp(-(x-params(MEAN2)).^2/(2 * ...
						  params(SD2)^2));
                      fity = fita+fitb;
end

if(state.logy)
  e = min(y(find(y > 0)))/1000000;
  diffs = (log(fity+e) - log(y+e)).^2;
else
  diffs = (fity - y).^2;
end

SSD = sum(diffs);
if(~isnan(minoffset) & params(MEAN) < minoffset)
  SSD = NaN;
end
