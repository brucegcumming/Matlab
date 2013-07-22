function resp = gauss2(params,x)

% gauss2(params,x) evalutates a Gaussian function of an input x
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

X = meshgrid(x,x);
Y = X';
r = sqrt(X.^2 + Y.^2);


xp = r - mean;
resp = base + amp .* exp(-(xp.^2)/(2 * sd .^2));


