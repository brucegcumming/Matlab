function dresp = 1dnoise(varargin)

nreps = 100;
pix2deg = 0.01;
npix = 256;

rnd = rand(nreps,npix);

Ra = Gabor([1 0.5],'pix2deg',0.15);
La = Gabor([1 0.5],'pix2deg',0.15);
Rb = Gabor([1 0.5 pi/4],'pix2deg',0.15);
Lb = Gabor([1 0.5 pi/4],'pix2deg',0.15);


rnd = rand(nreps,npix);
for disp = -20:10:20;
for j = 1:nreps
  iml = rnd(:,j);
  imr(1:end-disp) = iml(disp+1:end);
  resp(1,j) = Ra .* rnd(:,j) + La .* rnd(:,j);
  resp(2,j) = Rb .* rnd(:,j) + Lb .* rnd(:,j);
end
sqresp = resp .^2;
dresp(nd) = sum(sum(sqresp));
nd = nd+1;
end
  
