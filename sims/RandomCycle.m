function RandomCycle(p, loops)


for loop = 1:loops;
rnd = randn(1000,1);
rnd = rnd - mean(rnd);
pon = (rand(1000,1) > (1-p));
rndonoff = rnd .* pon .* sqrt(1/p);
res(loop,:) = abs((fft(rnd)));
ores(loop,:) = abs((fft(rndonoff)));
end
subplot(2,1,1);
hold off;
plot(rnd);
hold on;
plot(rndonoff,'g');
subplot(2,1,2);
hold off;
plot(mean(res,1));
hold on;
plot(mean(ores,1),'g');
