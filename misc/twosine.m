x = 0:0.01:1;
y = sin(x * 2 * pi);
z = sin(2 * x * 2 * pi);
plot(x,y);
hold on;
plot(x,z,'r');
plot(x,z-y,'g');