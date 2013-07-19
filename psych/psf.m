function psf = getpsf(probit)
     x(0) = 0;
     x(1) = 0.1;
   psf = fminsearch(@llike, x, probit);

function result = llike(params, probit);

mean = param[1];
sd = param[2];
y = ([probit.x] - mean)./sd;
p = erf(y);
q = 1-p;
result = sum(([probit.n] - [probit.resp]) .* log(q) - ([probit.n] .* log(p)));
