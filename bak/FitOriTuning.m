function [FitGauss,fval,exitflag] = FitOriTuning(x,y, varargin)
%[FitGauss,fval,exitflag] = FitOriTuning(x,y, varargin)
%Fits a double von mieses distribution to orientation tuning data that span
%360 degreees
BASE = 5;
AMP = 3;
SD = 2;
MEAN = 1;
AMP2 = 4;
maxiter = 5000;

if length(x)  < 2 || length(y) < 2
    fprintf('FitGauss:Insufficient Data\n');
    FitGauss = [];
    fval = 0;
    exitflag = -1;
    return;
end

if length(x) == length(y) && length(x) > 1
    id = find(~isinf(x));
    r = MeanVector(y,x,'double');
    guess(MEAN) = angle(r).*180/pi;
    guess(SD) = cv2sd(abs(r));
    guess(AMP) = std(y) * 2;
    guess(AMP2) = guess(AMP);
end

state.freebase = 0;
state.period = 0;
state.wrap = 0;
state.logy = 0;
state.maxamp = 0;
state.posamp = 0;
state.period = 360;
state.wrap = 1;

j = 1;
while j <= nargin -2
    if isstruct(varargin{j}) && isfield(varargin{j},'params')  % a fit struct
        infit = varargin{j};
        if isfield(infit,'state')
            state = infit.state;
        end
    elseif strncmpi(varargin{j},'freebase',5)
       state.freebase = 1;
       if length(x) == length(y) %fitting, not eval'ing, so need new guess
      id = find(~isinf(x));
      guess(BASE) = min(y(id));
      guess(AMP) = (max(y(id)) - min(y(id)))/1.2;
      %guess(MEAN) = mean(x(id))+ mean((x(id)-mean(x(id))) .* (y(id)-min(y(id))))/mean(y(id)-min(y(id)));
      %guess(SD) = std(((x(id)-mean(x(id))) .* (y(id)-guess(BASE))))/mean(y(id)-guess(BASE));
      nguess(BASE) = max(y(id));
      nguess(AMP) = -guess(AMP);
      %nguess(MEAN) = mean(x(id)) + mean((x(id)-mean(x(id))) .* (y(id)-max(y(id))))/mean(y(id)-max(y));
      nguess(SD) = guess(SD);
      guess(AMP2) = guess(AMP);
       end
  elseif strncmpi(varargin{j},'eval',4)
      [sse,FitGauss] = Minimise(y,x,zeros(size(x)),0,NaN,state);
      return;
  elseif strncmpi(varargin{j},'maxiter',5)
      j = j+1;
      maxiter = varargin{j};
  elseif strncmpi(varargin{j},'maxamp',5)
      state.maxamp = range(y) * 2;
  elseif strncmpi(varargin{j},'posamp',5)
      state.posamp = 1;
  elseif strncmpi(varargin{j},'posbase',5)
      guess(BASE) = min(abs(y));
      guess(AMP) = (max(y) - min(abs(y)))/1.2;
      state.freebase = 2;
  elseif strncmpi(varargin{j},'mean',7)
    j = j+1;
    guess(MEAN) = varargin{j};
  elseif strncmpi(varargin{j},'pos',3)
    state.posamp = 1;
  elseif strncmpi(varargin{j},'period',5)
    if strncmpi(varargin{j},'periodwrap',10)
        state.wrap  = 1;
    end
    j = j+1;
    state.period = varargin{j};
    if state.period == 360
        r = MeanVector(y,x);
        guess(MEAN) = angle(r);
        guess(SD) = cv2sd(abs(r));
        guess(AMP) = std(y) * 2;
    elseif state.period == 180;
        r = MeanVector(y,x,'double');
        guess(MEAN) = angle(r);
        guess(SD) = cv2sd(abs(r)) .* 1.5;
        guess(AMP) = std(y) * 2;
    end
  elseif strncmpi(varargin{j},'sd',2)
    j = j+1;
    guess(SD) = varargin{j};
  elseif strncmpi(varargin{j},'nreps',5)
    j = j+1;
    state.nreps = varargin{j};
  elseif strncmpi(varargin{j},'guess',5)
    j = j+1;
    guess = varargin{j};
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
[FitGauss.guessrss, guessfit] = Minimise(guess,x,y,0,NaN,state);
FitGauss.guess = guess;

options = optimset('MaxFunEvals',100000,'maxiter',maxiter,'display','off');
[fittedparams,fval,exitflag, output] = fminsearch(@Minimise,guess,options,x,y,sdmax,minoffset,state);
FitGauss.exitmsg = output.message;
if state.freebase == 1
 [nfittedparams,nfval,nexitflag, noutput] = fminsearch(@Minimise,nguess,options,x,y,sdmax,minoffset,state);
 if nfval < fval
     fittedparams = nfittedparams;
     fval = nfval;
     exitflag = nexitflag;
    FitGauss.exitmsg = noutput.message;
 end
end
    
if state.posamp
fittedparams(AMP) = abs(fittedparams(AMP));
end
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
fittedparams(MEAN) = mod(fittedparams(MEAN),state.period);
   [ssd, fitted] = Minimise(fittedparams, x,y, 0, NaN, state);

   FitGauss.mean = fittedparams(MEAN);
FitGauss.rss = fval;
if isfield(state, 'nreps')
    ymean = WeightedSum(y, state.nreps);
    yvar = sum((y-ymean).^2 .* state.nreps);
    FitGauss.pvar = 1 - fval/yvar;
else
    FitGauss.pvar = 1 - fval/var(y);
end
FitGauss.exit = exitflag;
FitGauss.fitted = fitted;
FitGauss.params = fittedparams; %% in order needed for Gauss(..)
FitGauss.state = state;


function [SSD, fity] = Minimise(params,x,y,sdmax,minoffset,state)

%BASE=4 AMP=1, SD=2, XCEN=3,
BASE = 5;
AMP = 3;
SD = 2;
MEAN = 1;
AMP2 = 4;

id = find(~isnan(x));
if state.period
    if state.wrap
        axd = x(id)-params(MEAN);
    xd(1,:) = axd;
    xd(2,:) = axd+state.period;
    xd(3,:) = axd-state.period;
    else
        xd = mod(state.period/2+x(id)-params(MEAN),state.period)-state.period/2;
    end
else
    xd = x(id)-params(MEAN);
end
if state.posamp
    params(AMP) = abs(params(AMP));
end
y = y(id);

if(length(params) == 5)
    if state.freebase == 2 %% baseline > 0
        params(BASE) = abs(params(BASE));
    end
  fity = params(BASE) + abs(params(AMP)) .* exp(-(xd).^2/(2 * ...
						  params(SD)^2)) +  abs(params(AMP2)).* exp(-(xd+180).^2/(2 * ...
						  params(SD)^2)); 
else
  fity = abs(params(AMP)) .* exp(-(xd).^2/(2 * ...
						  params(SD)^2)) + abs(params(AMP2)).* exp(-(xd+180).^2/(2 * ...
						  params(SD)^2));
end  

if state.wrap
    if length(params) == 5
        fity = sum(fity,1) - 2 * params(BASE);
    else
        fity = sum(fity,1);
    end
    if diff(size(fity)) == -diff(size(y))
        fity = fity';
    end
end
if(state.logy)
  e = min(y(find(y > 0)))/1000000;
  diffs = (log(fity+e) - log(y+e)).^2;
else
  diffs = (fity - y).^2;
end

if isfield(state,'nreps') && length(state.nreps) == length(diffs)
    SSD = sum(diffs .* state.nreps(id));
else
    SSD = sum(diffs);
end
if(state.maxamp > 0 & abs(params(AMP)) > state.maxamp)
  SSD = NaN;
end
if(sdmax > 0 & abs(params(SD)) > sdmax)
  SSD = NaN;
end
if(~isnan(minoffset) & params(MEAN) < minoffset)
  SSD = NaN;
end
