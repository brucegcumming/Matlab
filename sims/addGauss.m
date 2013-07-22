function addGauss(sep)

xv = -4:0.1:4;
x = gauss([-sep/2 1 1],xv);
y = gauss([sep/2 1 1],xv);
figure(1);
hold off;
plot(xv,x);
hold on;
plot(xv,y,':');
plot(xv,x+y,'r'),