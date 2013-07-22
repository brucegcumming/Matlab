function Filters(varargin)

%plot impulse response of a Butterworth Filter

fc = 1000; % Cut-off frequency (Hz)
fs = 8192; % Sampling rate (Hz)
order = 4; % Filter order
[B,A] = butter(order,2*fc/fs,'low'); % [0:pi] maps to [0:1] here
[B,A] = besself(order,2*fc/fs); % [0:pi] maps to [0:1] here
[sos,g] = tf2sos(B,A)
% sos =
%  1.00000  2.00080   1.00080  1.00000  -0.92223  0.28087
%  1.00000  1.99791   0.99791  1.00000  -1.18573  0.64684
%  1.00000  1.00129  -0.00000  1.00000  -0.42504  0.00000
% 
% g = 0.0029714
%
% Compute and display the amplitude response
Bs = sos(:,1:3); % Section numerator polynomials
As = sos(:,4:6); % Section denominator polynomials
[nsec,temp] = size(sos);
nsamps = 256; % Number of impulse-response samples
% Note use of input scale-factor g here:
x = g*[1,zeros(1,nsamps-1)]; % SCALED impulse signal
for i=1:nsec
  x = filter(Bs(i,:),As(i,:),x); % Series sections
end
%
plot(x); % Plot impulse response to make sure 
          % it has decayed to zero (numerically)
