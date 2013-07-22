load gabgauss.mat

hold off;
maxi = 60;

S = F(1:maxi);
subplot(1,2,1);
plot(S,OWR(1:maxi));
hold on;
plot(S,EWR(1:maxi),'r');
%plot(S,AWR(1:maxi),'r');
%plot(S,S * 1.1776 * 6);

ofit = polyfit(S,OWR(1:maxi),3)
fit = polyval(ofit, S);
%plot(S, fit, 'g');
efit = polyfit(S,EWR(1:maxi),3)
fit = polyval(efit, S);
%plot(S, fit, 'g');
kval = 34;
myfit = (S .* 1.1776 .* 6) + 1 ./ (1+S .^ 1.5 * kval);
%plot(S, myfit, 'g');
kval = 12

subplot(1,2,2);
K = (1 ./ (EWR(1:maxi) - (S .*1.1776 .* 6)) -1) ./S;
myfit = (S .* 1.1776 .* 6) + 1 ./ (1+(S .* K));
plot(S,K);
hold on;
efit = polyfit(S,K,2)
fit = polyval(efit, S);
plot(S, fit, 'g');

K = (1.4695 ./ (OWR(1:maxi) - (S .*1.1776 .* 6)) -1) ./S;
myfit = (S .* 1.1776 .* 6) + 1 ./ (1+(S .* K));
plot(S,K,'r');
hold on;
ofit = polyfit(S,K,2)
fit = polyval(ofit, S);
plot(S, fit, 'g');
%plot(S, myfit, 'c');
%Odd fit

K = (1.4695 ./ (AWR(1:maxi) - (S .*1.1776 .* 6)) -1) ./S;
myfit = (S .* 1.1776 .* 6) + 1 ./ (1+(S .* K));
%plot(S,K,'r');
afit = polyfit(S,K,2)
fit = polyval(afit, S);
%plot(S, fit, 'g');
%plot(S, myfit, 'c');
%Odd fit

kval = 24;
  myfit = (S .* 1.1776 .* 6) + 1.4695 ./ (1+ S .^ 1.5  * kval);
%  plot(S, myfit,'g');
kval = 9;
  myfit = (S .* 1.1776 .* 6) + 1.4695 ./ (1+ S .^ 1.1  * kval);
%  plot(S, myfit,'g');

%plot(OWR(1:maxi),EWR(1:maxi))