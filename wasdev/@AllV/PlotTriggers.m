function result = PlotTriggers(Vall, varargin)
%plot continuous props from FullV to look at poss
%Triggering Spaces
p = [];
result = [];

if isempty(p)
    p = 12;
end
sd = 0;



T(1,:) = [0 0 0 1 0 0 -2 0 0 1 0 0 0 0];
T(2,:) = [0 0 0 0 0 0 -2 0 0 2 0 0 0 0];
T(3,:) = [0 0 0 2 0 0 -2 0 0 0 0 0 0 0];
T(5,:) = [0 0 2 0 0 0 -2 0 0 0 0 0 0 0];
T(4,:) = [0 0 0 0 0 0 -2 0 0 2 0 0 0 0];
T(6,:) = [0 0 0 0 0 0 -2 0 0 0 2 0 0 0];
T(7,:) = [0 0 0 0 0 0 -2 0 0 1 0 1 0 0];
T(8,:) = [0 0 1 0 0 0 -2 0 0 0 1 0 0 0];
T(9,:) = [0 0 0 0 0 0 -2 0 0 0 0 2 0 0];
T(10,:) = Gauss(2,[-7:6]).*6;
T(11,:) = Gauss(1.5,[-7:6]).*6;
T(12,:) = Gauss(1,[-7:6]).*6;
rV = Vall.V(p,:);
if sd > 0
G = Gauss(2,-10:10);
rV = conv(rV,G,'same');
end
X = conv(rV,T(1,:),'same');

for j = 1:size(rV,1)
    sgn(j,:) = diff(sign(diff(X(j,:),1,2)),1,2);
end
    id = find(sgn < 0)+1; %minima
   vds = -10:1:20;
   id = id(10:end-10);
fprintf('%d events\n',size(sgn,2));
for j = 1:size(T,1)
    X = conv(rV,T(j,:),'same');
    dV = X(id);
    AllScores(:,j) = dV;
    sds(j) = std(dV);
    means(j) = mean(dV);
    ks(j) = kurtosis(dV);
    skews(j) = skewness(dV);
end
[a,b] = sort(ks,'descend');
X = conv(rV,T(b(1),:),'same');
Y = X(id);
X = conv(rV,T(b(2),:),'same');
Z = X(id);
DensityPlot(AllScores(b(1),:),AllScores(b(2),:));
[E, Ev] = eig(cov(AllScores));
pcs = AllScores * E;
PlotND(pcs(:,1:6),[]);

for j = 1:length(vds)
    dV = rV(id+vds(j))-rV(id);
    sds(j) = std(dV);
    means(j) = mean(dV);
    ks(j) = kurtosis(dV);
    skews(j) = skewness(dV);
end
        
hold off;
plot(sds,'o-');
hold on;
plot(means,'ro-');
plot(ks,'go-');
plot(skews,'ko');
    