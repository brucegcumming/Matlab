function resp = circgauss(params,x)

% circgauss(params,x) evalutates a circular Gaussian function of an input x, in radians.
% params is a vector defining the Gaussian whose elements are
% [mean  sd  amplitude [verticaloffset]]
%
% if params has only 1 element, this defines the SD, mean is taken
% to be zero and amplitude is 1/SD, so the area is 1. Useful for
% generating kernels for convolution
%
if(length(params) == 1)
  base = 0;
  mean = 0;
  sd = params(1);
  amp = 1/(sd * (2 * pi)^0.5);
else
  mean = params(1);
  sd = params(2);
  amp = params(3);
  if(length(params) > 3)
    base = params(4);
  else
    base = 0;
  end
end



for j = 1:5
    xp = x + ((j-2) * pi * 2) - mean;
    iresp(:,j) = base + amp .* exp(-(xp.^2)/(2 * sd .^2));
end
resp = sum(iresp,2);
