function xcorrtest(sd)

scale = 5;
z = rand(100,1) .* scale;
x = z+rand(100,1);
y = z+rand(100,1);
hold off;
xc = xcorr(x,y,50,'unbiased');
plot(xc);
hold on;
xs = smooth(x,5,'gauss');
ys = smooth(y,5,'gauss');
xc = xcorr(xs,ys,50,'unbiased');
plot(xc,'r');
xs = smooth(x,1,'gauss');
ys = smooth(y,1,'gauss');
xc = xcorr(xs,ys,50,'unbiased');
plot(xc,'g');
