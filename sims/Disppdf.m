function DispPDF(plottype, N)
if nargin < 2
N = 1400;
end
rn = randn(N,1);
d = 1.8+exp(-rn) * 10;
d = d(find(d > 0));
x = [0:3:200];
b = hist(d,x);
subplot(2,1,1);
bar(x,b);
angles = 180/pi * tan(0.065./d);
subplot(2,1,2);
hold off;
median(d)

if plottype == 1
hist(angles(find(d > 0)),100);
    return;
elseif plottype == 4
  a = median(angles);
  alla = hist(angles,100);
  na = hist(angles(find(angles > a)),100);
  fa = hist(angles(find(angles < a)),100);
  plot(conv(alla,na));
  hold on;
  plot(conv(alla,fa),'r');
end
dd = repmat(angles,1,N);
aa= repmat(sort(angles)',N,1);

if plottype == 2
x = [-2:0.02:2];
disp = dd-aa;
hist(disp(:),100);

elseif plottype == 3
disp = dd(:,1:N/2)-aa(:,1:N/2);
x = [-2:0.02:2];
near = hist(disp(:),x);
disp = dd(:,N/2:N)-aa(:,N/2:N);
far = hist(disp(:),x);
hold off;
plot(x,near);
hold on;
plot(x,far,'r');
    
end
