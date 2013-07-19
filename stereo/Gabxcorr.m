% Random dot patterns filtered through the FT of a simple cell at different
% SFs
clear all
close all


%****************************************************************
% Size of images (nxn)
n = 128;


% Number of sets of images to generate
nimages = 1;
% seeds to use in generating each image
randomseeds = [234 4324 90327 738566 9163759];


% Simple cell frequencies:
freq = 6/n;
nbw=1.5;
orbw=pi/5;


%****************************************************************


freqx = fftshift([-n/2:n/2-1])/n; 
freqy = fftshift([-n/2:n/2-1])/n; 
[freqx2d,freqy2d] = meshgrid(freqx,freqy);
absfreq = sqrt(freqx2d.^2+freqy2d.^2);
x = [-n/2:n/2-1]; y = [-n/2:n/2-1]; [x2d,y2d] = meshgrid(x,y);


% ---- GENERATE BASIC IMAGE IN FOURIER DOMAIN
im = ((rand(n)>0.5)-0.5)*2; % binary random-dot pattern with no DC
FT = fft2(im);


for jf = 1:length(freq)
    sx = sqrt(log(sqrt(2))) / pi ./ freq(jf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
    sy = sqrt(log(4)) / pi ./ freq(jf) / orbw ;
    Gabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( 2*pi.*x2d.*freq(jf));   
    OddGabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* sin( 2*pi.*x2d.*freq(jf));   
    FTGabor = fft2(Gabor);
    FTOddGabor = fft2(OddGabor);
        
    % make the filtered image
    image = fftshift(real(ifft2(FTGabor.*conj(FTGabor))));
    oimage = fftshift(real(ifft2(FTOddGabor.*conj(FTGabor))));
    aimage = fftshift(real(ifft2(FTOddGabor.*conj(FTOddGabor))));
        
    % ------- checking code ---------
    figure;
    imagesc(image); axis equal; 
colormap gray;
    end 

plot(x,image(64,:));
hold on;
plot(x,oimage(64,:));
plot(conv(Gabor(64,:),OddGabor(64,:)));
set(gca,'Xlim',[-64 64]);
figify(gcf,gca);
