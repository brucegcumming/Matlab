function dp = CalcDprime(x, y)
%dp = CalcDprime(x, y) Calc dprime for two distributions x and y
dp = (mean(x)-mean(y))./sqrt(mean([var(x) var(y)]));
