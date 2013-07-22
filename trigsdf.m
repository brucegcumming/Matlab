function [sdf, nsac, nspikes, spikes, times, details] = trigsdf(Trials, width, ...
						times, varargin)
% [sdf, nsac, nspikes, spikes, times] = trigsdf(Trials, width, times, ...)
%
% trigsdf takes a vector list of Trials, spikes
% and calculates a spike density function smoohted with 
% a Gaussian of S.D. 'width'.
% 'times' is a vector of timestamp values at which the density
% function is to be evaluated.
%
% sdf returns the smoothed function.
% nsac returns the number of sweeps making up this average.
% nspikes is the total number of spikes used.
% spikes is an unsmoothed function (a PSTH constructed with a bin
% width of 0.1ms) from which a tradional PSTH can be built.
%
% trigsdf(Trials, width, times, 'exp')
% Uses exponential of time constant width (if flag = 'exp') for smoothing.
%
% In each Trial structure, a vector of Trigger times may be
% specified. In this case, the spikes are aligned relative to the
% trigger time. If Trials.Trigger does not exist, all the trials
% are aligned at t=0. If more than one trigger time is specified,
% then than trial will contribute N times to the average (if there
% are N Trigger times.
%
% trigsdf(Trials, width, 'period', x)
% builds a psth of length period, to show average periodic
% modulation over period x.
%
% trigsdf(Trials, width, 'kernel', x)
% uses the kernel x for smoothing


%
% trigsdf(Trials, width)

% build up a PSTH with 0.1ms bins. Making the bins shorter might
% make this ultrafast.

period = 0;
flag = 'gauss';
clipping = 0; %% if 1, remove a full kernel width at each end.
nvar = nargin - 3;
j = 1;

while j  <= nvar
  str = varargin{j};
  if strmatch('kernel',str)
    j = j+1;
    kl = varargin{j};
  elseif strncmpi('clip',str,3)
      clipping = 1;
  elseif strmatch('period',str)
    j = j+1;
    period = round(varargin{j});
  elseif strncmpi('trig',str,3)
    j = j+1;
    triggers = varargin{j};
    if length(triggers) == 1
        [Trials.Trigger] = deal(triggers);
    else
        for k = 1:length(Trials)
            Trials(k).Trigger = triggers(k).trigger;
        end
    end
  elseif strmatch(str,{'gauss','box','raw','halfg','exp'})
      flag = str;
  end
  j = j+1;
end

if length(Trials) == 0
    sdf = 0;
    nsac = 0;
    nspikes = 0;
    spikes = [];
    return;
end

nbins = 1+times(end) - times(1);
if period
    spikes = zeros(period*2,1);
else
    spikes = zeros(nbins,1);
end

if(~isfield(Trials,'Trigger'))
  [Trials.Trigger] = deal(0);
end

nsac = 0;
allspks = [];
for j = 1:length([Trials.Trial])
    spkcount(j) = 0;
%bin is (spk - trigger) - times(1), = spk - (times(1)+trigger)
    if ~isempty(Trials(j).Trigger)
        for k  = length(Trials(j).Trigger):-1:1;
            t = Trials(j).Trigger(k) + times(1);
            spk = fix([Trials(j).Spikes'-t]);
            idx = find(spk > 0 & spk <= nbins);
            if period
                spk = fix(mod(spk(idx)-t,period))+1;
                idx = 1:length(spk);
            end
% This would become nbins * spk(idx)/duration
%    bins = ceil(nbins .* spk(idx) ./ nbins);
%
            if(~isempty(idx))
                spikes(spk(idx)) = spikes(spk(idx)) +1;
  %              allspks = [allspks spk(idx)];
               if(period)
                    spikes(spk(idx)+period) = spikes(spk(idx)+period) +1;
   %                 allspks = [allspks spk(idx)+period];
               end
            end
            if(period)
                nsac = nsac + (times(end)-times(1))/period;
            else
                nsac = nsac + 1;
            end
 %           spkcount(j) = spkcount(j) + length(idx);
        end
    end
end

details.spkc = spkcount;
details.allspks = allspks;

if(period)
%  spikes = [spikes(1:period); spikes(1:period)+period];
end


if(strmatch('exp',flag))
  x = 0:width*3;
  kl = exp(-x/width);
  rates = conv(spikes,kl);
  fullsdf = [rates(1:end-width*3) .* 10000/(width * nsac)];
elseif(strmatch('gauss',flag))
  x = -width*3:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  if clipping
      fullsdf = [rates(6*width:end-6*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
      clipping = 6 * width;
      tclip(1) = 3 * width;
  else
      fullsdf = [rates(3*width:end-3*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
  end
elseif(strmatch('halfg',flag))
  x = -0:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  
  fullsdf = [rates(1:end-3*width) .* 2 * 10000/(sqrt(2 * pi) * ...
						  width * nsac)];
elseif(strmatch('box',flag))
  kl = ones(width,1);
  rates = conv(spikes,kl);
  if clipping
      fullsdf = [rates(width:end-width+1) .* 10000/(sqrt(2 * pi) * ...
						  width * nsac)];
                      clipping = width;
      tclip(1) = width/2;
  else
      fullsdf = [rates(1:end-width+1) .* 10000/(width * nsac)];
  end
elseif(strmatch('kernel',flag))
  kl = varargin{1};
  rates = conv(spikes,kl);
  w = floor(length(kl)/2);
  fullsdf = [rates(w:end-w) .* 10000/nsac];
elseif(strmatch('raw',flag))
  fullsdf = spikes;  
end
if clipping
    idx = (times - times(1));
    idx = idx(find(idx > clipping)) -clipping;
    times = times(1) + tclip(1) + idx;
else
    idx = 1 + (times - times(1));
end
    if period
    w = round(period/2);
    idx = [period:period+w period+1-w:period-1];
end

idx = idx(find(idx <= length(fullsdf)));
sdf = fullsdf(idx);
nspikes = sum(spikes);


    
  
  