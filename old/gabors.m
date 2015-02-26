sd=20;
f=0.015;
amp = 10;
x = -1000:1000;
y = amp * cos(2 * x *f * pi) .* exp(-(x.^2)/(2 * sd *sd));
z = fft(y);
sum = x .*0;

hold off;
subplot(2,1,1);
plot(x,y);
hold on


n = 0;
for offset = -20:5:20;
p = x+offset;
y = amp * cos(2 * p *f * pi) .* exp(-(p.^2)/(2 * sd *sd));
plot(x,y);
sum = sum+y;
n = n+1;
end
sum = sum/n;
plot(x,sum,'r');
hold off
subplot(2,1,2);
plot(abs(fftshift(z)));
hold on
z = fft(sum);
plot(abs(fftshift(z)),'r');
hold off

