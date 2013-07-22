subplot(1,1,1);

alpha = [1:25] .* pi/180;
z = 2 * pi/180;

scales = [0.25 0.5 1 2 4];
colors = mycolors(0);
for j = 1:length(scales);
x = z*scales(j);
theta = atan(z ./ alpha);
btheta = atan(x ./ alpha);
plot(theta,btheta,'color',colors{j});
hold on;
end

