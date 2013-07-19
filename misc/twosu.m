function twosu(varargin)

x = -100:100;
c(1,:) = sin(x .* 2 .* pi .* 0.01);
a(1,:) = c(1,:) * -0.5;
c(2,:) = 0.5 + 0.5 * sin(x .* 2 .* pi .* 0.01 + pi/2);
a(2,:) = 1 - c(2,:);

hold off;
plot(c(1,:));
hold on; 
plot(c(2,:),':');
plot(a(1,:),'r');
plot(a(2,:),'r:');
plot(sum(a),'m');
plot(sum(c),'k');