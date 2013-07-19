function [disps, dresp, allresps, antiresps, rnd] = rls(nreps, varargin)

showplot = 0;

pix2deg = 0.01;
pix2deg = 0.1;
contrast = [1 1];
normalize = [0 0];
ngain = [0 0 0];
ralfplot = 1;
dx = 0;
showrf = 0;
sd = 2;
rnd = [];
disps = [];
j = 1;
while j < nargin
    if strncmpi(varargin{j},'plot',4)
        showplot = 1;
    elseif strncmpi(varargin{j},'disps',5)
        j = j+1;
        disps = varargin{j};
    elseif strncmpi(varargin{j},'dx',2)
        j = j+1;
        dx = varargin{j}/2;
    elseif strncmpi(varargin{j},'lc',2)
        j = j+1;
        contrast(2) = varargin{j};
    elseif strncmpi(varargin{j},'pix',3)
        j = j+1;
        pix2deg = varargin{j};
    elseif strncmpi(varargin{j},'rc',2)
        j = j+1;
        contrast(1) = varargin{j};
    elseif strncmpi(varargin{j},'rfplot',4)
        showrf = 1;
        if strncmpi(varargin{j},'rfplotonly',8)
        showrf = 2;
        end
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

Ra = Gabor([1 sd 0 1 -dx],'pix2deg',pix2deg);
La = Gabor([1 sd 0 1 dx],'pix2deg',pix2deg);
Rb = Gabor([1 sd pi/2 1 -dx],'pix2deg',pix2deg);
Lb = Gabor([1 sd pi/2 1 dx],'pix2deg',pix2deg);

if showrf
    x = ([1:length(Ra)]-length(Ra)/2).*pix2deg;
    GetFigure('BEM RFs');
    hold off;
    plot(x,Ra);
    hold on;
    plot(x,Rb,'g');
    plot(x,La,'b:');
    plot(x,Lb,'r');
    if showrf == 2
        return;
    end
end
npix = 512;
if isempty(rnd)
    rnd = rand(npix,nreps) - 0.5;
end

binary = 1;
if binary
    rnd(find(rnd < 0)) = -1;
    rnd(find(rnd > 0)) = 1;
end

normscale(1) = sqrt(1+ ngain(1) * contrast(1).^2);
normscale(2) = sqrt(1+ ngain(2) * contrast(2).^2);
normscale(3) = sqrt(1+ngain(3) * mean(contrast.^2));

nd = 1;
if ~isempty(disps)
    disps = round(disps./pix2deg);
    maxdisp = ceil(length(disps)/2);
    maxdisp = max([1-min(disps) 1]); 
else
    maxdisp = 20;
    disps = -(maxdisp-1):1:(maxdisp-1);
end

if ralfplot
    GetFigure('PosSums');
    hold off;
    GetFigure('PhaseSums');
    hold off;
    colors = mycolors;
end
for disp = disps;
for j = 1:nreps
  iml = contrast(2) .* rnd(maxdisp:maxdisp+length(Ra)-1,j)'./normscale(2);
  imr = contrast(1) .* rnd(maxdisp+disp:maxdisp+disp+length(Ra)-1,j)'./normscale(1);

  rresp(1,j) = sum(Ra .* imr);
  lresp(1,j) = sum(La .* iml);
  rresp(2,j) = sum(Rb .* imr);
  lresp(2,j) = sum(Lb .* iml);
end

aresp = (lresp - rresp)./normscale(3);
resp = (lresp + rresp)./normscale(3);
sqresp = resp .^2;
dresp(nd) = sum(sum(sqresp));
disps(nd) = disp;
allresps(nd,:) = sum(sqresp);
antiresps(nd,:) = sum(aresp.^2);
if ralfplot
    GetFigure('PosSums')
    c = 1+mod(nd-1,size(colors,2))
    plot(rresp(1,:)+lresp(1,:),rresp(2,:)+lresp(2,:),'.','color',colors{c});
    hold on;
    GetFigure('PhaseSums')
    c = 1+mod(nd-1,size(colors,2))
    plot(rresp(1,:)+lresp(2,:),-rresp(2,:)+lresp(1,:),'.','color',colors{c});
    hold on;
end

nd = nd+1;
end

disps = disps .* pix2deg;

if showplot
    plot(disps,dresp);
end