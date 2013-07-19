function cv = sd2cv(sd)
%
% converts sd in degrees to vector length
% Turns out this mapping is a Gaussian with sd 28.64 degrees

cv = exp(-(sd.^2)./(2 .* 28.6462^2));
if sd < 0
    cv = cv.* -1;
end
%specific to current implementation of stimuli
if abs(sd) > 120
    cv = 0;
end