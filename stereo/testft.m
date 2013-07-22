epos = 20;
n = 256;

leftim = [zeros(1,(n/2)+epos) ones(1,(n/2)-epos)]; %%single step edge
lFT = fft(leftim);


jf = 1;
nbw = 1.5;
freq = 6/n;
xsd = n/10;
x = [1:length(leftim)]-n/2;
sx = sqrt(log(sqrt(2))) / pi ./ freq * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
Gabor = exp(-x.^2./2./sx^2) .* cos( 2*pi.*x.*freq);


FTGabor = fft(Gabor);
limage = real(ifft(conj(lFT).*FTGabor));
hold off;
plot(limage,'r');
hold on;
plot(leftim,'b');