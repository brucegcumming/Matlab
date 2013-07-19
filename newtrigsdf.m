function [sdf, nsac, nspikes] = newtrigsdf(Trials, width, times,flag)
% trigsdf takes a vector list of spike times, spikes
% and a set of trigger times,
% and calculates a spike density function with either 
% a Gaussian of S.D. width  
% or an exponential of time constant width (if flag = 'exp')
% in steps of step (timestamp units).
%

% build up a PSTH with 0.1ms bins. Making the bins shorter might
% make this ultrafast.

duration = times(end) - times(1);
nbins = 1+(duration);
spikes = zeros(nbins,1);
  if(~isfield(Trials,'Trigger'))
    [Trials.Trigger] = deal(0);
  end

nsac = 0;
for j = 1:length([Trials.Trial])
  for t = [Trials(j).Trigger]
    spk = [Trials(j).Spikes'-t-times(1)];
    idx = find(spk > 0 & spk < nbins);
% This would become nbins * spk(idx)/duration
%    bins = ceil(nbins .* spk(idx) ./ nbins);
%
    if(~isempty(idx))
      spikes(spk(idx)) = spikes(spk(idx)) +1;
    end
    nsac = nsac + 1;
  end
end

k = 1;
for j = 1:10:length(spikes)-10
  subspk(k) = sum(spikes(j:j+10));
  k = k+1;
end
spikes = subspk;

binw = width/10;
if(strmatch('exp',flag))
  x = 0:10:width*3;
  kl = exp(-x ./ width);
  rates = conv(spikes,kl);
  fullsdf = [rates(1:end-binw*3) .* 10000/(width * nsac)];
elseif(strmatch('gauss',flag))
  x = -width*3:10:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  fullsdf = [rates(3*binw:end-3*binw) .* 10000/(sqrt(2 * pi) * width * nsac)];
end
idx = 1 + (times - times(1))/10;
%sdf = fullsdf(idx);
sdf = fullsdf;
nspikes = sum(spikes);



    
  
  