function [sqresp, im, ids] = stc1d(nreps, varargin)

showplot = 0;

pix2deg = 0.01;
pix2deg = 0.1;
contrast = [1 1];
normalize = [0 0];
ngain = [0 0 0];
dx = 0;
sd  = 0.3;
rnd = [];
disps = [];
j = 1;
npix = 16;
prc = 80' %percentile of resps to exceed threshold

while j < nargin
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'disps',5)
        j = j+1;
        disps = varargin{j};
    elseif strncmpi(varargin{j},'dx',2)
        j = j+1;
        dx = varargin{j}/2;
    elseif strncmpi(varargin{j},'20',2)
        j = j+1;
        contrast(2) = varargin{j};
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        contrast(1) = varargin{j};
    elseif strncmpi(varargin{j},'norm',2)
        j = j+1;
        ngain = varargin{j};
        
    elseif strncmpi(varargin{j},'rnd',2)
        j = j+1;
        rnd = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end

Ra = Gabor([1 sd 0 1 -dx],'pix2deg',pix2deg,'npts',npix);
La = Gabor([1 sd 0 1 dx],'pix2deg',pix2deg,'npts',npix);
Rb = Gabor([1 sd pi/4 1 -dx],'pix2deg',pix2deg,'npts',npix);
Lb = Gabor([1 sd pi/4 1 dx],'pix2deg',pix2deg,'npts',npix);

if isempty(rnd)
    rnd = rand(npix*2,nreps) - 0.5;
end

binary = 1;
if binary
    rnd(find(rnd < 0)) = -1;
    rnd(find(rnd > 0)) = 1;
end

for j = 1:nreps
  iml = rnd(1:npix,j)';
  imr = rnd((npix+1):end,j)';
  im.lpixels(j,:) = iml;
  im.rpixels(j,:) = imr;
  rresp(1,j) = sum(Ra .* imr);
  lresp(1,j) = sum(La .* iml);
  rresp(2,j) = sum(Rb .* imr);
  lresp(2,j) = sum(Lb .* iml);
end

resp = (lresp + rresp);
sqresp = mean(resp .^2,1);

th = prctile(sqresp, prc);
ids = find(sqresp > th);
