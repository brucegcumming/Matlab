function c = BimodalCoeff(x, e)
e = 1.3;
c = (1+skewness(x).^2)./((kurtosis(x).^e)+3);


