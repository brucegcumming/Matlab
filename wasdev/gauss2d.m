function [X,Y,Z] = gauss2(params,x, varargin)

% [X,Y,Z] = gauss2d(params,x) evalutates a Gaussian function of an input x
% params is a vector defining the Gaussian whose elements are
% [mean  sd  amplitude [verticaloffset]]
%
% if params has only 1 element, this defines the SD, mean is taken
% to be zero and amplitude is 1/SD, so the area is 1. Useful for
% generating kernels for convolution
%

j = 1;
mean(2) = params(1);
setmean = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'mean',2)
        j = j+1;
        mean = varargin{j};
        setmean = 1;
    elseif strncmpi(varargin{j},'ymean',2)
        j = j+1;
        mean(2) = varargin{j};
        setmean = 1;
    end
    j = j+1;
end

if(length(params) == 1)
  base = 0;
  if setmean == 0 
      mean(1) = 0; mean(2) = 0;
  end
  sd(1) = params(1);
  sd(2) = params(1);
  amp = 1/(sd(1) * sd(2) * (2 * pi)^0.5);    
elseif(length(params) == 2)
  base = 0;
  if setmean == 0 
      mean(1) = 0; mean(2) = 0;
  end
  sd = param
  amp = 1/(sd(1) * sd(2) * (2 * pi)^0.5);    
else
  mean(1) = params(1);
  sd(1) = params(1);
  sd(2) = params(2);
  amp = params(3);
  if(length(params) > 3)
    base = params(4);
  else
    base = 0;
  end
end

[X, Y] = meshgrid(x,x);
xp = X - mean(1);
yp = Y - mean(2);
r = sqrt(xp.^2 + xp.^2);


Z = base + amp .* exp(-(xp.^2)/(2 * sd(1) .^2)) .* exp(-(yp.^2)/(2 * sd(2) .^2));


