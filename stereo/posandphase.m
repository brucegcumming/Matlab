function posandphase(seedid,varargin)
% Random dot patterns filtered through the FT of a simple cell at different
% SFs

%****************************************************************
% Size of images (nxn)
n = 128;
row = 64;
col = 64;
deg2pix =20;
sf = 1/deg2pix;
sd = 1* deg2pix;
STIMGABOR = 1;
figa = 'Pos/Phase figure';
stim = 0;
showim = 0;
slope = 0.0; % set this to plot a line on top of responses to look at the orientation of the ridge.
dp = 0;
    qianmode = 0;    
colors = mycolors;
% Number of sets of images to generate
nimages = 3;
normalize = 0;
posonly = 0;
% seeds to use in generating each image
randomseeds = [234 4324 90327 738566 9163759 235 ];
rand('seed',4234);
sffilter = 0;
showfilter = 0;
jf = 4; %default harmonic
jsf = 4;
freq = [4 6]/n;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'image',4)
        showim = 1;
    elseif strncmpi(varargin{j},'filter',4)
        sffilter = 1;
    elseif strncmpi(varargin{j},'freq',4)
        j = j+1;
        jf = varargin{j};
    elseif strncmpi(varargin{j},'dp',2)
        j = j+1;
        dp = varargin{j};
    elseif strncmpi(varargin{j},'gabor',4)
        stim = STIMGABOR;
        if length(varargin) > j & isnumeric(varargin{j})
        end
    elseif strncmpi(varargin{j},'norm',4)
        normalize = 1;
    elseif strncmpi(varargin{j},'poso',4)
        posonly = 1;
    elseif strncmpi(varargin{j},'qian',4)
        qianmode = 1;
    elseif strncmpi(varargin{j},'sf',2) %stimulus SF for Gabor, as harmonic
        j = j+1;
        jsf = varargin{j};
    elseif strncmpi(varargin{j},'slope',4)
        j = j+1;
        slope = varargin{j};
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        figa = varargin{j};
    end
    j = j+1;
end
% Simple cell frequencies:
nbw=1.5;
orbw=pi/5;
sf = freq(1);


%****************************************************************

phase = 0;
freqx = fftshift([-n/2:n/2-1])/n; 
freqy = fftshift([-n/2:n/2-1])/n; 
[freqx2d,freqy2d] = meshgrid(freqx,freqy);
absfreq = sqrt(freqx2d.^2+freqy2d.^2);
x = [-n/2:n/2-1]; y = [-n/2:n/2-1]; [x2d,y2d] = meshgrid(x,y);

freq = absfreq;
sf = freq(jsf);
% ---- GENERATE BASIC IMAGE IN FOURIER DOMAIN
im = ((rand(n)>0.5)-0.5)*2; % binary random-dot pattern with no DC
FT = fft2(im);
GetFigure(figa);


if sffilter
    jf = 2;
r = sqrt(x2d.^2 + y2d.^2);
sr = 2;
rp = 5;
sFT = fftshift(exp(-(r-rp).^2/(2 * sr^2)));
subplot(1,1,1);
end
    sx = sqrt(log(sqrt(2))) / pi ./ freq(jf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
    sy = sqrt(log(4)) / pi ./ freq(jf) / orbw ;
    Gabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( 2*pi.*x2d.*freq(jf));   
    OddGabor = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* sin( 2*pi.*x2d.*freq(jf));   
    FTGabor = fft2(Gabor);
    FTOddGabor = fft2(OddGabor);

disp = 0;
    hold off;
    if showfilter
subplot(2,1,1);
imagesc(Gabor);
subplot(2,1,2);
    else
subplot(1,1,1);
        
    end
for k = seedid:seedid
    if stim == STIMGABOR

        [x,y] = meshgrid(1:n,1:n);
        im = cos(2 * pi * sf.*x) .* exp(-((x-n/2).^2/(2 * sd^2) + (y-n/2).^2/(2 * sd^2)));
        imr = cos(2 * pi * sf.*x + dp) .* exp(-((x-n/2).^2/(2 * sd^2) + (y-n/2).^2/(2 * sd^2)));

    else
        if k > length(randomseeds)
            rand('seed',k);
        else
            rand('seed',randomseeds(k));
        end
        im = ((rand(n)>0.5)-0.5)*2; % binary random-dot pattern with no DC
        imr = im;
    end
    FT = fft2(im);
    FTr = fft2(imr);
    if sffilter
        FT = FT .* sFT;
        FTr = FTr .* sFT;
    end
    % make the filtered image
    check = 0;
    if check
        fim = real(ifft2(FT.*FT));
        subplot(2,2,1);
        imagesc(im);
        subplot(2,2,2);
        imagesc(fim);
    end
    image = real(ifft2(FT.*FTGabor));
    oimage = real(ifft2(FT.*FTOddGabor));
    w = size(image,2);
j = 1;
if posonly
    phases = 0;
    midrow = 1;
else
    phases = -180:5:180;
    midrow = 36;
end
for phase = phases .* pi/180;
    Gabora = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( phase+2*pi.*x2d.*freq(jf));   
    OddGabora = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* sin( phase + 2*pi.*x2d.*freq(jf));   
    FTGabora = fft2(Gabora);
    FTOddGabora = fft2(OddGabora);
    aimage = real(ifft2(FTr.*FTGabora));
    oaimage = real(ifft2(FTr.*FTOddGabora));
% this way does not keep constant cyclopean location
if qianmode
    a = (image(row,col) .* aimage(row,:));
    b = oimage(row,col) .* oaimage(row,:);
% this way does.
else
    a = (image(row,:) .* fliplr(aimage(row,:)));
    b = oimage(row,:) .* fliplr(oaimage(row,:));
end    
if normalize
%im(rwo,col) is one eyes pixel;
%im(row,:) is the row of pixels in the other eye.so this gives a vector of
%normalizing constants for each of possible position disparities.
    norm = (image(row,col).^2 + image(row,:).^2 + oimage(row,:).^2 + oimage(row,col).^2);
    %if true disparity is 0, then the symmetry means it doesn't matter
    %which eye is flipped.
    norm = (image(row,:).^2 +  fliplr(image(row,:).^2) + oimage(row,:).^2 + fliplr(oimage(row,:).^2));
else
    norm = 1;
end
	    corrmap(j,:) = (a+b)./norm;
        amap(j,:) = a./norm;
        bmap(j,:) = b./norm;
	    j = j+1;
end 

if showim
    hold off;
    imagesc(corrmap);
    colormap('hot');
    if slope
        pphases = -3600:10:3600;
        ppix = 64 + (pphases * slope);
        hold on;
        plot(ppix, (mod(pphases+180,360))/mean(diff(phases)));
    end
else
    plot(corrmap(midrow,:),'color',colors{k});
    hold on;
    plot(amap(midrow,:),'color',colors{2});
    plot(bmap(midrow,:),'color',colors{3});
end
end
xlabel('Disparity (pixels)')
ylabel('correlation');



%figify(gcf,gca);
%set(gca,'Xlim',[-90 90]);
