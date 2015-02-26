function [theta, c, details] = BestAngles(x,y, varargin)
%[theta, c, details] = BestAngles(x,y, varargin)
%
j = 1;
test = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'hartigan',3)
        test = test+2;
    end
    j = j+1;
end
if test == 0
    test = 1;
end
a = 0:pi/36:pi;
  
  for j = 1:length(a)
      xy = xyrotate(x,y,a(j));
      if bitand(test,2)
      dip(j) = HartigansDipTest(sort(xy(:,1)));
      end
      if bitand(test,4)
      [aa,bb] = FindDip(xy(:,1));
      mydip(j) = bb.dipsize(1);
      end
      if bitand(test,1)
      coeff(j) = (1+skewness(xy(:,1)).^2)./(kurtosis(xy(:,1).^1.3)+3);
      end
      stds(j) = std(xy(:,1));
      ystds(j) = std(xy(:,2));
  end
  if test == 2
      coeff = dip;
  end
  [c, j] = max(coeff);
  theta = a(j);
  details.coeff = coeff;
  details.angles = a;
  details.xy = xyrotate(x,y,a(j));
details.stds = stds;
details.ystds = ystds;
