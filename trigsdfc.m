function [sdf, nsac, nspikes, spikes, times, details] = trigsdf(Trials, width, ...
						times, varargin)
% [sdf, nsac, nspikes, spikes, times, details] = trigsdf(Trials, width, times, ...)
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
% modulation over period x. Here, 'times' specifies the time range from
% each trial that will be included. The sdf returned is evaluated at  0.1ms
% intervals
%
% trigsdf(Trials, width, 'kernel', x)
% uses the kernel x for smoothing


%
% trigsdf(Trials, width)

% build up a PSTH with 0.1ms bins. Making the bins shorter might
% make this ultrafast.

skipevent = 0;
nskip = 0;
period = 0;
flag = 'gauss';
clipping = 0; %% if 1, remove a full kernel width at each end.
nvar = nargin - 3;
j = 1;
lfpmin = [];
timebinned = 1;
getallspks = 0;
countwin = [];

while j  <= nvar
  str = varargin{j};
  if strmatch('kernel',str)
    j = j+1;
    kl = varargin{j};
  elseif strncmpi('allspks',str,5)
      getallspks = 1;
  elseif strncmpi('events',str,3)
      for k = 1:length(Trials)
          Trials(k).Trigger = Trials(k).Events{:,2};
      end
  elseif strncmpi('clip',str,3)
      clipping = 1;
  elseif strncmpi('countwin',str,7)
      j = j+1;
      countwin = varargin{j};
  elseif strncmpi('freetimes',str,7)
      timebinned = 0;
  elseif strmatch('lfpmin',str)
    j = j+1;
    lfpmin = varargin{j};
    j = j+1;
    samplerate = varargin{j};
  elseif strmatch('period',str)
    j = j+1;
    period = round(varargin{j});
  elseif strncmpi('skipev',str,5)
      skipevent= 1;
  elseif strncmpi('trig',str,3)
    j = j+1;
    triggers = varargin{j};
    if iscell(triggers) & length(triggers) 
    for k = 1:length(Trials)
        Trials(k).Trigger = round(triggers{k});
    end
    else
    for k = 1:length(Trials)
        Trials(k).Trigger = round(triggers(k).trigger);
    end
    end
  elseif ischar(str) & strmatch(str,{'gauss','box','raw','halfg','exp'})
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

if IsAllExpt(Trials)
    details.times = times;
    for j = 1:length(Trials.Header)
        Expt = All2Expt(Trials,Trials.Header(j).cellnumber);
        [sdf(:,j), b,c,d,e,f] = trigsdfa(Expt.Trials,width, times, varargin{:});
        details.nsac = b;
        details.nspikes(j) = c;
        details = CopyFields(details,f);
    end
    nsac = details;
    return;
end


ntimepts = length(times);
times = [times(1)-width times];
nbins = round(1+times(end) - times(1));

if(~isfield(Trials,'Trigger'))
  [Trials.Trigger] = deal(0);
end

nsac = 0;
allspks = [];
spkcount = [];
%bin is (spk - trigger) - times(1), = spk - (times(1)+trigger)
if ~isfield(Trials,'Trial')
    for j = 1:length(Trials)
        Trials(j).Trial = j;
    end
end

counts = {};

allspk = [];
for j = 1:length([Trials.Trial])
% if skipevent is set, then spikes occuring with the first stimulus period
% after and event are NOT included. - allows the mean PSTH uncontaminated
% by skips to be calculated
if ~isempty(lfpmin)
          lfptimes = Trials(j).lfpo + round(Trials(j).Spikes .* samplerate);
      end
        if skipevent & ~isempty(Trials(j).Events)
            skiper = Trials(j).Events{1,2};
            if skiper >= times(1) & skiper <= times(end)
                nskip = 1;
            else
                nskip = 0;
            end
            idx = find(Trials(j).Spikes < skiper | Trials(j).Spikes > skiper+period);
            txs = repmat(round(times(1) + Trials(j).Trigger),length(idx),1);
            sxs = repmat(Trials(j).Spikes(idx),1,length(Trials(j).Trigger));
        elseif isempty(Trials(j).Spikes)
            spk = [];
        else
            spk = bsxfun(@minus,Trials(j).Spikes,Trials(j).Trigger);
        end
        details.triggers{j} = Trials(j).Trigger;
%used to use fix(sxs-txs), but now 1bin = 1 clock tick, no need.
%timbinned can be zero if using free array of times
%using histc here instead of find/for seems to be slower, prob becuase its
%sparse
%       counts = histc(spk(:),1:nbins);
%        allspk = [allspk spk(:)'];
        %idx = find(spk > 0 & spk <= nbins);
%        spkcount(j) = length(idx);
        if period
%      idx = find(spk >= 0 & spk <= (times(end) - times(1)));
%leave fix in here - period may not be an integer
            spk = fix(mod(spk(idx),period))+1;
            spk = [spk; spk + period];
            idx = 1:length(spk);
        end
% This would become nbins * spk(idx)/duration
%    bins = ceil(nbins .* spk(idx) ./ nbins);
% 

% The for loop is faster than spikes(spk(idx)) = spikes(spk(idx))+1
% And
       allspks{j} = spk(:);
%        allspks = [allspks; spk(idx)];
        if(period)
            nsac = nsac + length(Trials(j).Trigger) * (times(end)-times(1))/period;
            if nskip
                nsac = nsac - nskip;
            end
        else
        end
        for k = 1:size(countwin,1)
            counts{j,k} = sum(spk > countwin(k,1) & spk <= countwin(k,2),1);              
        end
        
        if ~isempty(countwin) && length(counts{j,k}) < length(Trials(j).Trigger)
            if isempty(Trials(j).Spikes)
                counts{j,k} = zeros(1,length(Trials(j).Trigger));
            end
        end
%Do this AFTER test for Trigger, Spikes. If spikes is empty, still count
%trigger
        nsac = nsac + length(Trials(j).Trigger);
end

 spkt = times(1):times(end); %time of bins in spikes

spikes = histc(cat(1,allspks{:}),spkt);
details.spkc = spkcount;
details.allspks = allspks;
if ~isempty(counts)
    details.counts = counts;
end
if(period)
%  spikes = [spikes(1:period); spikes(1:period)+period];
end

%spikes has dimensionns N,1, need to make sure
% kernel is also this way, to be sure the convolution is also
%Otherwise if spiesk si short than kernel, they swap

if(strmatch('exp',flag))
  x = 0:width*3;
  kl = exp(-x/width);
  rates = conv(spikes,kl,'same');
  fullsdf = [rates(1:end-width*3) .* 10000/(width * nsac)];
  spkt = spkt(1:end-3*width);

elseif(strmatch('gauss',flag))
  x = -width*3:width*3;
  kl = exp(-(x'.^2)/(2 * width^2));
  rates = conv(spikes,kl,'same');
  if clipping
      fullsdf = [rates(6*width:end-6*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
      clipping = 6 * width;
      tclip(1) = 3 * width;
      spkt = spkt(6*width:end-6*width);
  else
      fullsdf = [rates(3*width:end-3*width) .* 10000/(sqrt(2 * pi) * ...
          width * nsac)];
      spkt = spkt(3*width:end-3*width);
  end
elseif(strmatch('halfg',flag))
  x = -0:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  
  fullsdf = [rates(1:end-3*width) .* 2 * 10000/(sqrt(2 * pi) * ...
						  width * nsac)];
elseif(strmatch('box',flag))
  width = round(width);
  kl = ones(width,1);
  rates = conv(spikes,kl,'same');
  if clipping
      tclip(1) = width/2;
      clipping = width;
  end
  if nsac == 0
      fullsdf = [rates(width:end)];
      spkt = spkt(width:end);
  elseif clipping
      fullsdf = [rates(width:end)] .* 10000/(width * nsac);
      clipping = width;
      spkt = spkt(width:end);
  else
      fullsdf = [rates .* 10000/(width * nsac)];
  end
elseif(strmatch('kernel',flag))
  kl = varargin{1};
  rates = conv(spikes,kl);
  w = floor(length(kl)/2);
  fullsdf = [rates(w:end-w) .* 10000/nsac];
  spkt = spkt(2:end-w);
elseif(strmatch('raw',flag))
  fullsdf = spikes;  
end
%times(1) is an extra element added at the start, to account for
%smoothing with. So length of sdf should be length(times)-1 at this
%point
if clipping
    idx = 1+(times - times(1));
    idx = idx(find(idx > clipping)) -clipping;
    times = times(1) + tclip(1) + idx;
else
    idx = 1+(times - times(1));
    idx = idx(2:end);
end


if period
    w = round(period/2);
    idx = [period:period+w period+1-w:period-1];
end

idx = floor(idx(find(floor(idx) <= length(fullsdf))));
%if length(idx) < ntimepts
%    idx = [idx(1)-1 idx];
%end
sdf = fullsdf(idx);
times = spkt(idx);
nspikes = sum(spikes);
if getallspks
    spikes = allspks;
end

    
  
  