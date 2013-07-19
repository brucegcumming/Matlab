% I wanted to check if a single energy-unit with a Gaussian RF
% looked similar to a complex cell with even and odd Gabors.
clear all
close all

% Set up coords:
x = [-20:0.2:20];
y = [-20:0.2:20];
nx = length(x);
ny = length(y);
[x2,y2]=meshgrid(x,y);

% I set up a filter in Fourier space: the FT of a Gabor:
SD = 2;
BW = 1.5;
SF = 1./2/pi./SD.*sqrt(log(2)).*(2^BW+1)/(2^BW-1); % in cycles per deg
FTRF = fft2(exp(-0.5*(x2.^2+y2.^2)./SD^2).*cos(2*pi*SF.*x2));
RF = exp(-0.5*(x2.^2+y2.^2)./SD^2).*cos(2*pi*SF.*x2);

% I make an image
im = zeros(ny,nx);
im(20:22,30:32)=1;
im(40:42,10:12)=-1;
im(60:62,70:72)=1;
im(120:122,130:132)=1;

%FTRF = ones(size(im));

% I convolve the image with the filter in Fourier space:
filtim = ifft2(fft2(im) .* FTRF, 'symmetric');
%filtim = real(ifft2(fft2(im)));

figure
subplot(2,2,1)
imagesc(x,y,im)
title('original image')
subplot(2,2,2)
imagesc(x,y,fftshift(filtim))
title('fftshifted filtered im')
subplot(2,2,3)
imagesc(x,y,RF)
title('Filter')
subplot(2,2,4)
imagesc(x,y,filtim)
title('"filtered" image')

for j=1:4;subplot(2,2,j); axis equal tight;end