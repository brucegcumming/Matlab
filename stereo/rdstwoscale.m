% Random dot patterns filtered through the FT of a simple cell at different
% SFs
clear all
close all


%****************************************************************
% Size of images (nxn)
n = 128;
row = 64;
col = 64;


colors = mycolors;
% Number of sets of images to generate
nimages = 1;
% seeds to use in generating each image
randomseeds = [234 4324 90327 738566 9163759];
rand('seed',randomseeds(5));

% Simple cell frequencies:
freq = [2 4]/n;
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

jf = 1;
kf = 2;
    sx = sqrt(log(sqrt(2))) / pi ./ freq(jf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
    sy = sqrt(log(4)) / pi ./ freq(jf) / orbw ;
    Gabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( 2*pi.*x2d.*freq(jf));   
    OddGabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* sin( 2*pi.*x2d.*freq(jf));   
    FTGabor = fft2(Gabor);
    FTOddGabor = fft2(OddGabor);
    
    sx = sqrt(log(sqrt(2))) / pi ./ freq(kf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
    sy = sqrt(log(4)) / pi ./ freq(kf) / orbw ;
    Gabora = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( 2*pi.*x2d.*freq(kf));   
    OddGabora = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* sin( 2*pi.*x2d.*freq(jf));   
    FTGabora = fft2(Gabora);
    FTOddGabora = fft2(OddGabora);
    
    
        
    figure;
for j = 1
    im = ((rand(n)>0.5)-0.5)*2; % binary random-dot pattern with no DC
    FT = fft2(im);
    % make the filtered image
    image = real(ifft2(FT.*FTGabor));
    oimage = real(ifft2(FT.*FTOddGabor));

    imagea = real(ifft2(FT.*FTGabora));
    oimagea = real(ifft2(FT.*FTOddGabora));
        
    % ------- checking code ---------

%    cross = image .* image + oimage .* oimage;
    colormap gray
    a = (image(row,:) .* image(row,col));
    b = oimage(row,:) .* oimage(row,col);
    norm = (image(row,:).^2 + image(row,col).^2) + (oimage(row,:).^2 + oimage(row,col).^2);

    corrmap(j,:) = 2 * (a+b)./norm;
    aa = (imagea(row,:) .* imagea(row,col));
    ba = oimagea(row,:) .* oimagea(row,col);
    norma = (imagea(row,:).^2 + imagea(row,col).^2) + (oimagea(row,:).^2 + oimagea(row,col).^2);

    plot((1:length(a)) - col,2 *(a+b)./norm,'color',colors{j});
   hold on;
    plot((1:length(aa)) - col,2 *(aa+ba)./norma,'color',colors{j+1});
    
end 

xlabel('Disparity (pixels)')
ylabel('correlation');


%plot((1:length(a)) - col,mean(corrmap,1),'color','k');

figify(gcf,gca);
set(gca,'Xlim',[-64 64]);
