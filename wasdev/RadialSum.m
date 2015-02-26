function [sor, details] = RadialSum(kernel, fixr)
%sor = RadialSum(kernel, fixr)
%
%Calculate the radial sum of an image, with inteporation.
%Default image size 256x256
%default resolution for sum: 1 deg
%if fixr = 0 or fix is omitted, then sum is caclulated up to radius of
%image width/2

if nargin == 1
    fixr = 0;
end

angles = [0:pi/180:pi];
coss = cos(angles);
sins = sin(angles);
%pixel 129,129 is the center of a 256x256 fourier transform;
%unless its been flipped..
linrvals = [0.5:0.5:size(kernel,2)/4];
cc = findmid(kernel);

for j = 1:length(linrvals)
    r = linrvals(j);
    xi = r .* coss;
    yi = r .* sins;
    sir(j,:) = interp2(kernel, cc(1)+xi, cc(2)+yi);
    rwts(j) = r;
end
r= fixr;
details.maxr = max(linrvals);

if fixr > 0
    rwts(r > fixr) = 0;
    details.maxr = fixr;
end
sor = rwts * sir;
details.sfk = sum(sir,2);
details.ork = sor;
details.weighted =1;
details.angles = angles;
details.center = cc;
details.sir = sir;

function sor = OldRadialSum(kernel, fixr)

coss = cos([0:pi/180:pi]);
sins = sin([0:pi/180:pi]);
%pixel 129,129 is the center of a 256x256 fourier transform;
linrvals = [0.5:0.5:45.5];

for j = 1:length(linrvals)
    r = linrvals(j);
    xi = r .* coss;
    yi = r .* sins;
    sir(j,:) = interp2(kernel, 129+xi, 129+yi);
    rwts(j) = r^2;
end
r= fixr;
sor = sum(sir);