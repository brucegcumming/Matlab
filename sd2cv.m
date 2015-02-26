function cv = sd2cv(sd)
%
% converts sd in degrees to vector length (1-circular variance)
% Turns out this mapping is a Gaussian with sd 28.64 degrees

cv = exp(-(sd.^2)./(2 .* 28.6462^2));
cv(sd < 0) = cv(sd < 0) .* -1; %not sign in case sd = 0
%specific to current implementation of stimuli
cv(abs(sd) > 120) = 0;
    