function [theta, c, details] = BestAngle(x,y, test, varargin)
%
% Find best angle to cut in a 2D space. With very skewed distributison
% (Energy, ADC value where it is clipped), the bimodality coeffeicient is
% misleading. But otherwise itts smoother and more reliable. 

a = 0:pi/36:pi;

for j = 1:length(a)
xy = xyrotate(x,y,a(j));
if bitand(test,2)
dip(j) = HartigansDipTest(sort(xy(:,1)));
end
skews(j) = skewness(xy(:,1));
kurts(j) = kurtosis(xy(:,1));
coeff(j) = cmb.BimodalCoeff(xy(:,1),1.5);
end

details.bmc = coeff;
if bitand(test,1)
details.coeff = coeff;
if bitand(test,2)
details.dip = dip;
details.hdip = max(dip);
end
elseif bitand(test,2)
details.hdip = dip;
details.coeff = dip;
end
[c, j] = max(details.coeff);
details.besti = j;
theta = a(j);
details.angles = a;

%if using bimodality coeff to find best angle, calc Hartigan for best
%angle so can compare with other measures.
if test == 1
xy = xyrotate(x,y,theta);
details.dip = HartigansDipTest(sort(xy(:,1)));
end

