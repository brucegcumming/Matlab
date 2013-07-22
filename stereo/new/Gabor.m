function y = Gabor(x, sd, cycles, phase)
% 
% 	y = Gabor(x, sd, cycles, phase, dc)
%

filtersize = size(x,2);
f = Cycles2const(cycles,filtersize);
y = exp( - (x.^2/sd.^2)) .* (cos(f .* x + phase));
