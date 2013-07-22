sd=1;
f=0.025;
amp = 10;

EW = [];
OW = [];
EWR = [];
OWR = [];
AWR = [];
F = [];

for f = 0.001:0.01:1.001;

%kval = 0.022/sd; good for power of 1.3
%kval = 0.036/sd; good for power of 1.2
kval = 0.063;
crit = 5;
if f < 0.2
  off = 2;
  X = 0:0.0001:5;
else
  off = 0.1/f;
  X = 0:0.0001:1/(2*f);
end

G = amp .* cos(-pi/2 + 2 .* X .* f .* pi) .* exp(-(X.^2)/(2 * sd ...
						  *sd));
A = amp .* cos(-9*pi/20 + 2 .* (X-off) .* f .* pi) .* exp(-((X-off).^2)/(2 * sd ...
						  *sd));
COS = amp .* cos(-pi/2 + 2 .* X .* f .* pi);
E = amp .* cos(2 .* X .* f .* pi) .* exp(-(X.^2)/(2 * sd *sd));
[peak, opeaki] = max(G);

D = abs(G-peak/2);
[half, halfi] = min(D(1:opeaki));
a = halfi;
[half, halfi] = min(D(opeaki:length(X)));
b = halfi + opeaki-1;

[peak, peaki] = max(E);
D = abs(E-peak/2);
[half, halfi] = min(D(1:length(X)));
c = halfi;



[peak, apeaki] = max(A);
D = abs(A-peak/2);
[half, halfi] = min(D(1:apeaki));
e = halfi;
[half, halfi] = min(D(apeaki:length(X)));
g = halfi + apeaki-1;

subplot(1,2,1);
hold off;
plot(X, G,'b');
hold on;
plot([X(a) X(b)],[G(a) G(b)],'b'); 
plot([X(e) X(g)],[A(e) A(g)],'r'); 
plot([0 X(c)],[E(c) E(c)],'k'); 
plot(X, A,'r');
plot(X, E,'k');
%plot([1/(12 * f) 5/(12 * f)],[amp/2 amp/2],'r'); 
drawnow;
%fprintf('Odd HW at HH %.3f %.3f %.3f %.3f',X(opeaki) - X(b),X(a), ...
%	X(opeaki),X(b));
%fprintf(' FW at HH %.3f\n',X(b) - X(a));
fprintf('F %.2f HWHH %.3f %.3f Cos %.3f\n',f,(X(b) - X(a))/2,(X(g) - X(e))/2,1/(6 * f));
%A Gaussian has a halfwidth at half-height of 1.1774 * SD. So any
%given halfwidth, w  is equivalent to A Gaussian with an SD (sd')
%of w/1.1776.
sda = (X(b) - X(a))/(2 * 1.1774);

EW= [EW X(c)];
OW= [OW (X(b) - X(a))/2];
EWR= [EWR sd/(X(c)/1.1774)];
OWR= [OWR sd/((X(b) - X(a))/(2 * 1.1774))];
AWR= [AWR sd/((X(g) - X(e))/(2 * 1.1774))];
F = [F f * sd];
end

Xa = X - X(opeaki);
Ga = amp .* exp(-(Xa.^2)/(2 * sda * sda));
plot(X,Ga,'k');

subplot(1,2,2);

hold off;
%plot sd/sd'  vs f.sd, this is then scale invariant, but the same
%as f vs 1/sigma'
plot(F,OWR);
hold on;
plot(F,EWR,'r');
plot(F,AWR,'k');
plot(F,F * 1.1776 * 6);
save /d/bgc/matlab/gabgauss/gabgauss.mat F OWR EWR AWR;
ofit = polyfit(F,OWR,3)
fit = polyval(ofit, F);
%plot(F, fit, 'g');
efit = polyfit(F,EWR,3)
fit = polyval(efit, F);
%plot(F, fit, 'g');
fsd = 0.2149;
fit= exp(-(F.^2)/(2 * fsd *fsd));
%plot(F,fit,'g');
%plot(F,(1.1776 * 6 .* F),'b');
