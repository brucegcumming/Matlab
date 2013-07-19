function [r, pval, theta] = Rayleigh(angles, varargin)

nresample = 0;
pval = NaN;
j = 1;
while j < nargin
    if strncmpi('resample',varargin{j},3)
        j = j+1;
        nresample = varargin{j};
    end
    j = j+1;
end

sina = mean(sin(angles));
cosa = mean(cos(angles));
r = abs(cosa + i * sina);
theta = atan2(sina,cosa);

if nresample
   rnd = rand(length(angles),nresample) .* 2 * pi;
   rs = abs(i * mean(sin(rnd)) + mean(cos(rnd)));
   pval = sum(rs >= r)./nresample;
end
   