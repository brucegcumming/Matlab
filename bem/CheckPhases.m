function f = CheckPhases(phase, freq)

forward = 1;
if forward
    r = BuildGabors(phase, freq);
    return;
end
colors= 'rgb';

hold off;

for j =1:length(phase)
g = Gabor([freq 16 phase(j)],'npts',1024);
f = fft(fftshift(g));
id = find(abs(f(1:128)) > max(abs(f))/10);
[y,x] = hist(angle(f(id)));
[a,b] = unique(x)
subplot(2,1,1);
h = bar(x,y);
set(h,'facecolor',colors(j));
hold on;
subplot(2,1,2);
plot(g,colors(j));
hold on;


end

function f = BuildGabors(phase, freqs)

if isempty(freqs)
freqs = [0.002 0.0015 0.001 0.0005];
end
[nr,nc] = Nsubplots(length(freqs));

for j = 1:length(freqs)
    freq = freqs(j);
subplot(nr,nc,j); 
hold off; 
x = [1:1024] - 1025/2;
sx = 150;
g = exp(-(x).^2/(2.*sx.^2)) .* cos(2 * pi * freq .* (x)+phase(1));
fg  =g;
%g = Gabor([freq 16 0],'npts',1024);
%fg = exp(-(x).^2/(2.*sx.^2)) .* cos(2 * pi * freq .* (x) + phase(1));
f = fft(fftshift(g));
ffg = fft(fftshift(fg));
n = 512;
id = 1:512;
f(id) = abs(f(id)) .* cos(phase) +abs(f(id)) * i * sin(phase);
id = 513:1024;
f(id) = abs(f(id)) .* cos(-phase) + abs(f(id)) * i * sin(-phase);
g = fftshift(real(ifft(f)));
plot(abs(f),abs(ffg));
plot(g);
hold on;
plot(fg,'r');
set(gca,'xlim',[1 1024]);
end