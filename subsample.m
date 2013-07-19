function subsmp = subsample(x, ratio)

% subsmp = subsample(x, ratio)
% returns a shorter vector, with length length(x)/ratio
% by taking the mean of groups of values of x.

nsub = floor(length(x)/ratio);
maxl = nsub * ratio;
smp = reshape(1:maxl,ratio,nsub);
subsmp = mean(x(smp));