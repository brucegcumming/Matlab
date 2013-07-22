function MFfft(varargin)

x = 1:1024;
amps = zeros(512,1);
A = 1;
phases = pi/2 + [1:512] .* pi/(2 * 512);
y(1:512) = -1;
y(513:1024) = 1;
ft = fft(y(1:1024));
for c = 1:2:511;
    amps(c) = A;
    A = A/3;
    fti(c,:) = cos(phases(c) + (c * 2 * pi * (x)./1024)) * amps(c);
    ffti(c,:) = cos(phases(c) + (c * 2 * pi * (x)./1024)) * abs(ft(c+1));
end
plot(sum(fti));



for c = 1:2:511;
    fti(c,:) = cos(angle(ft(c+1)) + (c * 2 * pi * (x-1)./1024)) * abs(ft(c+1));
end
plot(sum(fti));
zf = zeros(size(ft));
zf(2) = ft(2);
zf(1024) = ft(2);
ft(2) = 0;
ft(1024) = 0;
z = ifft(ft);
f = ifft(zf);
plot(z);
hold on; 
plot(x,f,'r');

plot(sum(fti));