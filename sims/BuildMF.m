function BuildMF(npts,varargin)
%BuildMF(npixels,varargin)
%Build a square wave/MF by Fouriere synthesis
x = 1:npts;
halfpts = floor(npts/2);
amps = zeros(halfpts,1);
A = 1;
%there is a simple rule for phase
phases = pi/2 + [1:halfpts] .* pi/(2 * halfpts);
%the rule for amplitude is tricky. The standard rule from continuous
%transform (reduce by a factor of 3 for each new odd harmonic) fails
%for frequencies about half the nyquist frequency. So use an FFT to
%determine correct amplitudes.
y(1:halfpts) = -1;
y(1+npts-halfpts:npts) = 1;
ft = fft(y);

for c = 1:2:floor((npts+1)/2);
    ffti(c,:) = cos(phases(c) + (c * 2 * pi * (x)./npts)) * abs(ft(c+1));
%also build one with usual amplitude rule for comparison.
    amps(c) = 1/c;
    fti(c,:) = cos(phases(c) + (c * 2 * pi * (x)./npts)) * amps(c);
end
plot(sum(ffti));

