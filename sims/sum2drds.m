function sum2drds(varargin)
%
%make an RDS by summing two 1d noise pattens



seeds = 0;
npix = 256;
seed = 0;

j = 1;
while j <= length(varargin)
    j= j+1;
end


for j = 1:length(seeds)
    seed = seeds(j);
v = repmat(round(rand(1,npix)),npix,1)-0.5;
h = repmat(round(rand(npix,1)),1,npix)-0.5;
hb = repmat(round(rand(npix,1)),1,npix)-0.5;
vb= repmat(round(rand(1,npix)),npix,1)-0.5;
im = v+h;
imb = v-h;
imc = v+hb;
xc = corrcoef(im(:),imb(:))
xc = corrcoef(im(:),imc(:))
scale = 0.5/max(max(abs(im)));
im = 0.5 + im .* scale;
imwrite(im,sprintf('TestIm-%dL.pgm',seed),'PGM');
scale = 0.5/max(max(abs(imb)));
im = 0.5 + imb .* scale;
imwrite(im,sprintf('TestIm-%dR.pgm',seed),'PGM');
end
