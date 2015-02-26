function y = Gauss(x, gauss)

if(isfield(gauss,'logfit') & gauss.logfit)
  y = gauss.base + gauss.amp .* exp(-(log(x)-gauss.mean).^2./(2 * ...
						  gauss.sd^2));
else
  y = gauss.base + gauss.amp .* exp(-(x-gauss.mean).^2./(2 * ...
						  gauss.sd^2));
end

