sd=2;
sdy=2;
f=0.025;
amp = 10;


f = 0.00;
for crit = 1:9
x = 0;
%%first get width at y = 0;
val = 100;
while(val > crit)
val = amp * cos(2 * x * f * pi) * exp(-(x.^2)/(2 * sd *sd));
x = x + 0.1;
end
width = x;
val = 100;

y = 0;
while(val > crit)
val = amp *  exp(-(y.^2)/(2 * sdy *sdy));
y = y + 0.1;
end
height = y;
fprintf('C %.1f ratio %.1f\n',crit,height/width);
end

%kval = 0.022/sd; good for power of 1.3
%kval = 0.036/sd; good for power of 1.2
kval = 0.063;
crit = 5;
X = [];
Y = [];
Z = [];
M = [];
O = [];

for f = 0.0:0.002:0.04
%for f = 0.0:0.005:0.04
%for f = 0.0:0.05:0.4
x = 0;
%%first get width at y = 0;
val = 100;
while(val > crit)
val = amp * cos(2 * x * f * pi) * exp(-(x.^2)/(2 * sd *sd));
x = x + 0.001;
end
width = x; %with at criterion value
val = 100;

x = f/4;
val = 100;
while(val > crit)
val = amp * cos(-pi/2 + 2 * x * f * pi) * exp(-(x.^2)/(2 * sd *sd));
x = x + 0.001;
end
w = x;
val = 100;

x = 0;
while(val > crit)
val = amp *  exp(-(x.^2)/(2 * sd *sd));
x = x + 0.001;
end
height = x;
r = height/width;
%a = (1.179 * sd)/((sd * 6.68  * f) +(f+1)/(f+1));

%with k = 0.038/sd this works well when f is small, but is too
%big for intermediate f.
sf = sd * f;
%a = (1.190 * sd)/((sf * 5.95) +(kval)/(f^(1.2)+kval));
a = (1.1776 * sd)/((sf * 7.04) +(kval)/(sf^1.2+kval));
%b = (1.190 * sd)/((sd * 5.95  * f^0.95));
fprintf('F %.3f ratio %.4f %.2f %.4f %.3f\n',f,height/width, height,width,a);
Y = [Y width/height];
X = [X f];
Z = [Z a];
O = [O w];
end
plot(X, Y);
      hold on;
      plot(X,Z,'k');
      plot(X,O,'b');
      hold off;
