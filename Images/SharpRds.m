function [im, left, right] = SharpRds(w, h, dx, varargin )


SIDEBYSIDE = 1;
SHARP = 2;
DOTS = 3;
showim = SHARP;
monoc = 0;

dw = round(w/2);
dh= round(h/2);
j = 1;
while j <= nargin - 3    
    if strncmpi(varargin{j},'dots',3)
        cscale = 1.0;
        [im, left, right] = DotIm(w, h, dx, 1, 400);
        tiffl = im2tiff(left, 1);
        tiffr = im2tiff(right, 1);
        im = tiffs2sharp(tiffl, tiffr);
        subplot(1,1,1);
        imagesc(im);
        return;
    end
    j = j+1;
    
end

left = round(rand(h,w));
right = left;
fg = round(rand(dh,dw));

sw = round(dw/2);
sh = round(dh/2);
for j = 1:dw
    right(sh:sh+dh-1,dx+j+sw) = fg(:,j); 
    left(sh:sh+dh-1,j+sw) = fg(:,j); 
end

if showim == SIDEBYSIDE
    subplot(1,2,1);
    imagesc(left);
    subplot(1,2,2);
    imagesc(right);
end
if monoc ==1
    right = zeros(size(right));
end
tiffl = im2tiff(left, 1);
tiffr = im2tiff(right, 1);
im = tiffs2sharp(tiffl, tiffr);

if showim == SHARP
    subplot(1,1,1);
    image(im);
end

function tiff = im2tiff(left, cscale)
tiff(:,:,1) = left .* cscale;
tiff(:,:,2) = left .* cscale;
tiff(:,:,3) = left .* cscale;

function im = tiffs2sharp(lim, rim)
k = 1;
im = zeros(size(lim,1),1+size(lim,2)*2,3);
w = size(lim,2);
for j = 1:w
    im(:,k,1) = lim(:,j,1);  %%Red, L
    im(:,k+1,1) = rim(:,j,1); %% Red, R
    im(:,k+1,2) = lim(:,j,1); %% Green, L
    im(:,k+2,2) = rim(:,j,1);
    im(:,k,3) = lim(:,j,1);
    im(:,k+1,3) = rim(:,j,1);
    k = k+2;
end

function [bim, left, right] = DotIm(w, h, dx, dy,ndots)

left = zeros(w,h);
right = zeros(w,h);

rnd = 10 + ceil(rand(ndots,2) .* (w-26));
for j = 1:ndots
    
    x = rnd(j,1);
    y = rnd(j,2);
    left(y:y+6,x:x+6) = 1;
    if(abs(x- w/2) < w/4 & abs(y - h/2) < h/4)
        right(y+dy:y+6+dy,x+dx:x+6+dx) = 1;
    else
        right(y+dy:y+6+dy,x:x+6) = 1;
    end
end

bim = [left right];