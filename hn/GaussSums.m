function [areas, peaks] = GaussSums(sds)

if isempty(sds)
    sds = [1 2 3 4 5 6 7 8 9 10 12 18 24 32 40 50];
end

for j = 1:length(sds)
    x = -500:500;
    y = gauss([0 sds(j) 1],x);
    y = y.^2;
    fy= abs(fftshift(fft(y)));
    peaks(j) = max(fy);
    areas(j) = trapz(fy);
end
plot(areas,peaks);