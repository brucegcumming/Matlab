function [sdf, nsac, nspikes] = oldtrigsdf(Trials, width, times,flag)
% trigsdf takes a vector list of spike times, spikes
% and a set of trigger times,
% and calculates a spike density function with either 
% a Gaussian of S.D. width  
% or an exponential of time constant width (if flag = 'exp')
% in steps of step (timestamp units).
%

% build up a PSTH with 0.1ms bins. Making the bins shorter might
% make this ultrafast.

nbins = 1+times(end) - times(1);
spikes = zeros(nbins,1);

nsac = 0;
for j = 1:length([Trials.Trial])
  if(~isfield(Trials(j),'Trigger'))
    Trials(j).Trigger = 0;
  end
  for t = [Trials(j).Trigger]
    spk = floor([Trials(j).Spikes'-t-times(1)]);
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

if(strmatch('exp',flag))
  x = 0:width*3;
  kl = exp(-x/width);
  rates = conv(spikes,kl);
  fullsdf = [rates(1:end-width*3) .* 10000/(width * nsac)];
elseif(strmatch('gauss',flag))
  x = -width*3:width*3;
  kl = exp(-(x.^2)/(2 * width^2));
  rates = conv(spikes,kl);
  fullsdf = [rates(3*width:end-3*width) .* 10000/(sqrt(2 * pi) * width * nsac)];
end
idx = 1 + (times - times(1));
sdf = fullsdf(idx);
nspikes = sum(spikes);


    
  
  