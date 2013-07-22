nreps = 1000;
ratio = 2;

x = randn(1,nreps);
y = randn(1,nreps) .* ratio;

chi = x.^2 + y.^2;
invchi = 1./(chi + 0.1);

hist(invchi);