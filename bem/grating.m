function [pwr, disps, dresp, details] = grating(contrast, ngain, varargin)
%[pwr, disps, dresp, details] = grating(contrast, ngain, varargin)
%Calculate effects of binoc/monoc normalization on disp tuning
%Contrast contains contrasts for [L R]
%
%ngain contains divisor for [L R B]
%
showplot = 0;
nreps = 100;
pix2deg = 0.01;
normalize = [0 0];
dphases = 0:0.1:pi*2;

j = 1;
dp = 0;
while j < nargin-1
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'lc',2)
        j = j+1;
        contrast(2) = varargin{j};
    elseif strncmpi(varargin{j},'phases',3)
        j = j+1;
        dphases = varargin{j};
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        contrast(1) = varargin{j};
    elseif strncmpi(varargin{j},'dp',2)
        j = j+1;
        dp = varargin{j};
    end
    j = j+1;
end

rfsf = 1;
p2x = 0.015;

Ra = Gabor([rfsf 0.5 0],'pix2deg',p2x);
La = Gabor([rfsf 0.5 0+dp],'pix2deg',p2x);
Rb = Gabor([rfsf 0.5 pi/4],'pix2deg',p2x);
Lb = Gabor([rfsf 0.5 dp + pi/4],'pix2deg',p2x);

x = [1:length(Ra)] .* p2x * rfsf * 2 * pi;
npix = 512;

normscale(1) = sqrt(1+ ngain(1) * contrast(1).^2);
normscale(2) = sqrt(1+ ngain(2) * contrast(2).^2);
normscale(3) = sqrt(1+ngain(3) * mean(contrast.^2));
normscale = ngain;
    nd = 1;
maxdisp = 21;
for k = 1:length(dphases)  % Can't worry about row/column
    disp = dphases(k);
    j = 1;
    for phase = 0:pi/8:2*pi
        iml = contrast(2) .* sin(x+phase+disp/2)./normscale(2);
        imr = contrast(1) .* sin(x+phase-disp/2)./normscale(1);
        rresp(1,j) = mean(Ra .* imr);
        lresp(1,j) = mean(La .* iml);
        rresp(2,j) = mean(Rb .* imr);
        lresp(2,j) = mean(Lb .* iml);
        j = j+1;
    end

resp = (lresp + rresp)./normscale(3);
sqresp = resp .^2;
dresp(nd) = sum(sum(sqresp));
disps(nd) = disp;
nd = nd+1;
end
  
if showplot
    h = plot(disps,dresp);
    if length(dphases) < 20
        set(h,'symbol','o');
    end
end

pwr(1) = famp(disps,dresp,1/(2*pi));
pwr(2) = mean(dresp);
details.gains = normscale;
details.resp = resp;
details.contrast = contrast;