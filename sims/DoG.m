function DoG(varargin)


x = -1000:1000;

sd = 50;
a = Gauss(sd, x);
b = Gauss(sd*5, x).*0.5;
y = a-b;
ft = abs(fft(y));
[a,b] = max(ft);
plot(ft);
period = 2000./(b-1);
ratio = sd./period

