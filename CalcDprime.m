function dp = CalcDprime(x, y)

dp = (mean(x)-mean(y))./sqrt(mean([var(x) var(y)]));
