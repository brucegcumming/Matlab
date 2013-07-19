function corrmap = oneedge(eposn,varargin)
% Random dot patterns filtered through the FT of a simple cell at different
% SFs

%****************************************************************
% Size of images (nxn)
n = 128;
row = 64;
col = 64;
figid = [];

%setting onoff = 1 manually looks at "on" and "off" even symmtric filters
%separately. 
onoff = 0;
showim = 1;
colors = mycolors;
% Number of sets of images to generate
nimages = 3;
normalize = 0;
% seeds to use in generating each image
randomseeds = [234 4324 90327 738566 9163759];
rand('seed',4234);

j = 1;
while j < nargin
    if strncmpi(varargin{j},'image',4)
        showim = 1;
    elseif strncmpi(varargin{j},'norm',4)
        normalize = 1;
    elseif strncmpi(varargin{j},'onoff',4)
        onoff = 1;
        showim = 0;
    elseif strncmpi(varargin{j},'figid',4)
        j = j+1;
        figid = varargin{j};
    elseif strncmpi(varargin{j},'showim',4)
        if isnumeric(varargin{j+1})
            j = j+1;
            showim = varargin{j};
        else
            showim = 1;
        end
    end
    j = j+1;
end
% Simple cell frequencies:
freq = [4 6]/n;
nbw=1.5;
orbw=pi/5;


%****************************************************************

phase = 0;
freqx = fftshift([-n/2:n/2-1])/n; 
x = [-n/2:n/2-1]; 

for ei = 1:length(eposn)
leftim = [zeros(1,(n/2)+eposn(ei)) ones(1,(n/2)-eposn(ei))]; %%single step edge
rightim = 1-leftim;
%lFT = fft(leftim);
%rFT = fft(rightim);


if ~isempty(figid)
    figure(figid);
else
    figure;
end

jf = 1;
    sx = sqrt(log(sqrt(2))) / pi ./ freq(jf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
    Gabor = exp(-x.^2./2./sx^2) .* cos( 2*pi.*x.*freq(jf));   
    OddGabor = exp(-x.^2./(2 .*sx^2)) .* sin( 2*pi.*x.*freq(jf));   
 %   FTGabor = fft(Gabor);
 %   FTOddGabor = fft(OddGabor);

disp = 0;
hold off;
k = 1;
% make the filtered image
limage = real(ifft(lFT.*FTGabor));
loimage = real(ifft(lFT.*FTOddGabor));
rimage = real(ifft(rFT.*FTGabor));
roimage = real(ifft(rFT.*FTOddGabor));
j = 1;
phases = -180:5:180;

if onoff
    lonimage = limage; 
    lonimage(find(limage < 0)) = 0;
    ronimage = rimage; 
    ronimage(find(rimage < 0)) = 0;
    loffimage = limage; 
    loffimage(find(limage > 0)) = 0;
    roffimage = rimage; 
    roffimage(find(rimage > 0)) = 0;
    a = lonimage(row) .* ronimage(:);
    b = loffimage(row) .* roffimage(:);
    plot(a);
    hold on;
    plot(b);
    plot(limage(row) .* rimage(:),'r');
    plot(limage,'g');
    return;
end
%
%Gabor, FTGabor are for the even symmetric/odd syymetric pair
% Gabora, OddGabora are a quadarture pair around different mean phases.
for phase = phases .* pi/180;
    Gabora = exp(-x.^2./2./sx^2).* cos( phase+2*pi.*x.*freq(jf));   
    OddGabora = exp(-x.^2./2./sx^2) .* sin( phase + 2*pi.*x.*freq(jf));   
    FTGabora = fft(Gabora);
    FTOddGabora = fft(OddGabora);
    limagea = real(ifft(lFT.*FTGabora));
    loimagea = real(ifft(lFT.*FTOddGabora));
    rimagea = real(ifft(rFT.*FTGabora));
    roimagea = real(ifft(rFT.*FTOddGabora));
    a = (limage(row) .* rimagea(:));
    b = loimage(row) .* roimagea(:);
        
if normalize
%im(rwo,col) is one eyes pixel;
%im(row,:) is the row of pixels in the other eye.so this gives a vector of
%normalizing constants for each of possible position disparities.
    norm = (limage(row).^2 + rimage(:).^2 + oimage(row).^2 + roimage(:).^2);
else
    norm = 1;
end
	    corrmap(j,:) = (a+b)./norm;
        limages(j,:) = limage;
        rimages(j,:) = rimagea;
        amap(j,:) = a;
        bmap(j,:) = b;
	    j = j+1;
end 

if showim == 1 
    hold off;
    imagesc(corrmap);
    colormap('hot');
    caxis([-50 50]);
elseif showim == 2
    hold off;
    imagesc(limage);
    colormap('hot');    
else
    plot(corrmap(36,:),'color',colors{k});
    hold on;
    plot(amap(36,:),'color',colors{2});
    plot(bmap(36,:),'color',colors{3});
end
xlabel('Disparity (pixels)')
ylabel('correlation');
drawnow;
end


%figify(gcf,gca);
%set(gca,'Xlim',[-90 90]);
