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
kernel{2} = kernel{1}.^2;
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
    tmp2 = conv(in(j,:),kernel{2});
    resp = tmp(kl:end-(kl+1));
    sqout(j,:) = tmp2(kl:end-(kl+1));
    out(j,:) = resp;
    rresp = resp;
    rresp(find(resp < 0)) = 0;
    oout(j,:) = (resp-min(resp)).^2;
    rout(j,:) = rresp;
    r2out(j,:) = rresp.^2;
    [amp(j), c] = famp(2500:7500,out(j,2500:7500),f(j));
    [oamp(j), c] = famp(2500:7500,oout(j,2500:7500),f(j));
    [sqamp(j), c] = famp(2500:7500,sqout(j,2500:7500),f(j));
    phase(j) = angle(c) - pi/2;
    [ramp(j), c] = famp(2500:7500,rout(j,2500:7500),f(j));
    gphase(j) = angle(c) - pi/2;
    [r2amp(j), c] = famp(2500:7500,r2out(j,2500:7500),f(j));
    ephase(j) = angle(c) - pi/2;
    hold on;
end

sqamp = amp.^2;
hold off;
subplot(1,1,1);
hold off;
plot(f,amp./max(amp),'k');
hold on;
plot(f,sqamp./max(sqamp),'b');
%plot(f,ramp./max(ramp),'r');
plot(f+0.0001,oamp./max(oamp),'g');
ylabel('Gain');

legend('Gaussian','ksq','HalfRect','Halfrectsq');
xlabel('Frequency');

