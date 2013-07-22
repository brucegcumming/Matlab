%calculate bode plots and step reponses for 
% Gaussian, Halfgaussian, and exponential kernels
%
%

len = 10000;
t = 1:len;
in = [];
out = [];
kl=200;
x = -kl:kl;
kernel = {};
kernel{1} = Gauss(20,-kl:kl);
kernel{2} = Gauss(40,-kl:kl);
kernel{2}(1:kl) = 0;
kernel{3} = exp(-x/40);
kernel{3}(1:kl) = 0;
for j = 1:length(kernel)
    kernel{j} = kernel{j} / sum(kernel{j});
end
hold off;
for j = 1:20
    f(j) = j/1000;
    in(j,:) = sin(t .* 2 .* pi .* f(j));
    tmp = conv(in(j,:),kernel{1});
    out(j,:) = tmp(kl:end-(kl+1));
    tmp = conv(in(j,:),kernel{2});
    gout(j,:) = tmp(kl:end-(kl+1));
    tmp = conv(in(j,:),kernel{3});
    eout(j,:) = tmp(kl:end-(kl+1));
%    plot(t',in(j,:));
    [amp(j), c] = famp(2500:7500,out(j,2500:7500),f(j));
    phase(j) = angle(c) - pi/2;
    [gamp(j), c] = famp(2500:7500,gout(j,2500:7500),f(j));
    gphase(j) = angle(c) - pi/2;
    [eamp(j), c] = famp(2500:7500,eout(j,2500:7500),f(j));
    ephase(j) = angle(c) - pi/2;
    hold on;
end

hold off;
subplot(3,1,1);
hold off;
plot(f,amp);
hold on;
plot(f,gamp,'r');
plot(f,eamp,'g');
ylabel('Gain');
subplot(3,1,2);
hold off;
plot(f,phase);
hold on;
plot(f,gphase,'r');
plot(f,ephase,'g')

legend('Gaussian','HalfGauss','Exp');
ylabel('Phase');
xlabel('Frequency');

subplot(3,1,3)
hold off;
step(1:1000) = 0;
step(1001:2001) = 1;

sout = conv(step,kernel{1});
plot(-kl:kl,sout(1000:1000+(2*kl)));
hold on;
sout = conv(step,kernel{2});
plot(-kl:kl,sout(1000:1000+2*kl),'r');
sout = conv(step,kernel{3});
plot(-kl:kl,sout(1000:1000+2*kl),'g');
title('Step Response');
xlabel('Time');


