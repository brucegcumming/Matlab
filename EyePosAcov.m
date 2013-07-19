function [ac, aceach, avg] = EyePosAcov(Trials, eye, varargin)

% [ACF, ACFeach, avg] = EyePosAcov(Trials, ...)
%calculates the autocorrelation function for eye position signals
%averaged across a set of trials.
%
nvar = nargin - 3;
Ro = NaN;
state.scaling = 0;

j = 1;
while j <= nvar
    if isstr(varargin{j})
        if strmatch(varargin{j},'Noscale')
            state.scaling = 0;
        elseif strmatch(varargin{j},'Scale')
            state.scaling = 1;
        elseif strmatch(varargin{j},'Ro')
            j = j+1;
            Ro= pi * varargin{j}/180;
        elseif strncmpi(varargin{j},'skip',4)
            j = j+1;
            skiplen = varargin{j};
        end
    end
    j = j+1;
end

buflen = 1e10;
for trial = 1:length(Trials)
    buflen = min([buflen length(Trials(trial).Eyevals.lv)]);
end

k = 1;
for trial = 1:length(Trials)
  if eye == -1
      if isnan(Ro)
          array = Trials(trial).Eyevals.lv;
      else
          array = -Trials(trial).Eyevals.lv .* cos(Ro) + Trials(trial).Eyevals.lh .* sin(Ro);
      end
  else
      array = Trials(trial).Eyevals.rh - Trials(trial).Eyevals.lh;
  end
  array = array(1:buflen)';
  z = fft(array - mean(array));
  aceach(k,:) = real(ifft(z .* conj(z)));
  psths(k,:) = array;
  k = k+1;
end

avg = mean(psths,1);
ac = mean(aceach,1);

