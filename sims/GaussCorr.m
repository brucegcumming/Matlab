function GaussCorr(varargin)

twod = 0;
sd =10;
offsets = [0:10:200];
ths = [0 0.0000001 0.0001 0.001];
th = 0.01;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'th',2)
        j = j+1;
        ths = varargin{j};
    elseif strncmpi(varargin{j},'oned',4)
        twod = 0;
        ths = [0 0.000000001 0.0000001 0.00001 0.001];
    elseif strncmpi(varargin{j},'twod',4)
        twod = 1;
    end
    j = j+1;
end
if twod
[x,y] = meshgrid(-1000:1000,-1000:1000);
for j = 1:length(offsets)
x = -10000:10000;
[X,Y,g1] = Gauss2d([sd sd],x,'mean',[-offsets(j) 0]);
[X,Y,g2] = Gauss2d([sd sd],x,'mean',[offsets(j) 0]);

for k = 1:length(ths)
id = find(g1> ths(k) & g2 > ths(k));
xc = corrcoef(g1(id),g2(id));
c(j,k) = xc(1,2);
hold off;
plot(g1(:),g2(:),'o');
hold on;
plot(g1(id),g2(id),'ro');
drawnow;
end
end
else
x = -10000:10000;
for j = 1:length(offsets)
g1 = Gauss([-offsets(j) sd],x);
g2 = Gauss([offsets(j) sd],x);
g1 = normpdf(x, -offsets(j), sd);
g2 = normpdf(x, offsets(j), sd);
for k = 1:length(ths)
id = find(g1> ths(k) & g2 > ths(k));
if length(id) > 1
xc = corrcoef(g1(id),g2(id));
c(j,k) = xc(1,2);
else
    c(j,k) = NaN;
end
hold off;
plot(g1(:),g2(:),'o');
hold on;
plot(g1(id),g2(id),'ro');
end
drawnow;
end
plot(offsets,c);
legend(num2str(ths'));
end
hold off;
plot(offsets,c);
legend(num2str(ths'));
