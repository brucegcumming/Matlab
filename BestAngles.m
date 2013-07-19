function [theta, c, details] = BestAngles(x,y, varargin)

  a = 0:pi/36:pi;
  test = 1;
  
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
  end
  [c, j] = max(coeff);
  theta = a(j);
  details.coeff = coeff;
  details.angles = a;
