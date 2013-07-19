clear, clc

% neurons in each pool
nn = 100;

% total neurons
Nn = 8 * nn;

% trials count
Tn = 1000;

%initial avergage firing rate for population
PopmFr = 1;

% target correlation matrix
CorrMat = zeros(Nn,Nn)+0.1;
for g = 0: 7
    CorrMat(g*nn+1:g*nn+nn, g*nn+1:g*nn+nn) = 0.2;
%    CorrMat(Nn/2+1:end, Nn/2+1:end) = 0.5;
end
for i = 1: Nn, CorrMat(i,i) = 1; end

% actual trial by trial firing rates
FrsInit = normrnd(PopmFr, PopmFr, [Tn, Nn]);

SQCorr = sqrtm(CorrMat);

Frs = FrsInit * SQCorr;


figure(10), clf, 
imagesc(FrsInit); 
title ('initial population firing rate structre');

figure(11), clf, 
imagesc(Frs); 
title ('Generated firing rate structure');

cg = corr(Frs);
figure(12), clf,  imagesc(cg); 
title ('Generated Correlation structure');

figure(13), clf,
hist(FrsInit, 100, 'k'), hold on,
hist(Frs, 100, 'r')
title ('Freq distribution of initial (black) and final (red) firing rates ');

%%

a = find(mean(Frs(:, 1:Nn/2),2)>mean(Frs(:, Nn/2+1:Nn),2));
b = find(mean(Frs(:, 1:Nn/2),2)<mean(Frs(:, Nn/2+1:Nn),2));

for n = 1:Nn
%    CP(n) = ROCAUC(Frs(a,n), Frs(b,n));
end

%figure(14), clf, 
%plot(CP);
