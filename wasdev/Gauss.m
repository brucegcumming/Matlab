function resp = gauss(params,x,varargin)

% gauss(params,x) evalutates a Gaussian function of an input x
% params is a vector defining the Gaussian whose elements are
% [mean  sd  amplitude [verticaloffset]]
%
% if params has only 1 element, this defines the SD, mean is taken
% to be zero and amplitude is 1/SD, so the area is 1. Useful for
% generating kernels for convolution. !!Beware undersampling
%

period = 0;
j = 1;
while j  <= nargin-2
    if strncmpi(varargin{j},'period',4)
        j = j+1;
        period = varargin{j};
    end
    j = j+1;
end
if(length(params) == 1)
  base = 0;
  mean = 0;
  sd = params(1);
  amp = 1/(sd * (2 * pi)^0.5);
%  amp = 1/(sd);
else
  mean = params(1);
  sd = params(2);
  if length(params) > 2
      amp = params(3);
  else
      amp = 1/(sd * (2 * pi)^0.5);
  end
  if(length(params) > 3)
    base = params(4);
  else
    base = 0;
  end
end




xp = x - mean;
if period
    xp = mod(xp+period/2,period) - period/2;
end
resp = base + amp .* exp(-(xp.^2)/(2 * sd .^2));


