function [pwrs,apwrs] = getori(x,y,z)
%[pwrs,apwrs] = getori(x,y,z)  calulate power at combinations of ori/freq
%for arbitray x,y, locations. x,y, should be in the range +- 1 for the
%frequecny range to be sensible.
%set envtype = 0 to calculate fourier power at each, 
%set envtype = 1 to use Gabors (i.e. include Gaussian spatial weighting
%the choice of sd = 1./(4*f) is arbitrary, and may need tuning, especially
%for very low F (i.e. just a smooth gradient)/

j= 1;
plotenv = 0;
envtype = 0;
xm = 0;
ym = 0;

for f = 0.1:0.1:2
k = 1;
sd = 1./(4*f);
if plotenv
    [xi,yi] = meshgrid([-1:0.1:1],[-1:0.1:1]);
    env = exp(-((xi-xm).^2./sd + (yi-ym).^2./sd));
    %Z = interpf(x,y,z,xi,yi,1,0.2);
    Z = cos(2 * pi .* xi .* f) .* env;
    imagesc(Z);
end
if envtype == 1
    env = exp(-((x-xm).^2./sd + (y-ym).^2./sd));
else
    env = ones(size(x));
end
    
for ori = 0:pi/50:pi
   X  = (x-xm) .* cos(ori) + (y-ym) .* sin(ori);
   Y = (y-ym) .*cos(ori) - (x-xm) .* sin(ori);
   c = sum(env(:) .* z(:) .* exp(2*pi*i*f.*X(:)));
   ac = sum(exp(2*pi*i*f.*X(:)));
%  pwr = sum(cos(2 * pi * f .* X(:).*z(:).*env)).^2 + sum(sin(2 * pi * f .* X(:).*z(:))).^2;
   pwrs (j,k) = abs(c);
   apwrs (j,k) = abs(ac);
   k = k+1;
end
   j = j+1;
end