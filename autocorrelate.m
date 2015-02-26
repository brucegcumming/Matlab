function [ac, aceach, nspikes, psth] = autocorrelate(Trials, latency, ...
					       duration, varargin)

% [ACF, ACFeach, nspikes] = autocorrelate(Trials, latency, duration)
%
% autocorrelate calculates the autocorrelation function for
% Spike trains in a set of Trials. Spikes are counted from time
% 'latency' from the beginning of  a trial until duration+latency.
%
% Returns ac              Sum of ACs across trials.
%         aceach          array with one ACF for each trial
%         nspikes         number of spikes in each trial
%         psth            matrix with spikes for each trial      
%                         as a PSTH with a binwidth of 0.1ms. 
%
%
% autocorrelate(Trials, latency, duration, 'Scale')
% Scales autocorrelogram  to give approx mean spike rate
% (spikes/sec) by dividing by the number of spikes, then * 10000
% useful for comparing with PSTH. 

nvar = nargin - 3;
j = 1;
state.scaling = 0;
skiplen = 0;

while j <= nvar
    if isstr(varargin{j})
        if strmatch(varargin{j},'Noscale')
            state.scaling = 0;
        elseif strmatch(varargin{j},'Scale')
            state.scaling = 1;
        elseif strncmpi(varargin{j},'skip',4)
            j = j+1;
            skiplen = varargin{j};
        end
    end
    j = j+1;
end


k = 1;
for trial = 1:length(Trials)
  spikes = Trials(trial).Spikes(find(Trials(trial).Spikes > ...
		 latency+skiplen & Trials(trial).Spikes < duration+latency));
  array = zeros(1,ceil(duration));
  array(ceil(spikes-latency)) = 1;
  z = fft(array);
  aceach(k,:) = real(ifft(z .* conj(z)));
  nspikes(k) = length(spikes);
  psths(k,:) = array;
  k = k+1;
end

%array has units spikes/bin. So ACF of array has spikes^2/bin^2,
%summed over n bins, giving spikes^2/bin. Need to divide this by 
%n bins (length(Trials) * duration) to get ACF in spikes^2/bin^2;
% * 1e8 gives spikes^2/sec^2;

if(state.scaling)
  ac = 10000 .* sum(aceach,1) ./ sum(nspikes);
else
  ac = 1e8 .* sum(aceach,1) ./ (length(Trials) .* duration);
end

psth = mean(psths,1) .* 10000;

aceach = 1e8 .* aceach/duration;