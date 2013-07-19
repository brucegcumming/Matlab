function [resp, sd, r,p] = posandphasespace(nreps, dxs, varargin)

j = 1;
while j <= length(varargin)
    j = j+1;
end

%comparison is not trivial
%cant use sd, becuase the space is circular. So, for example if disparity
%is close to half a period, phase disparity is +- pi, generates large SD. 
%but can't use a cicular variacne readily either. RF position plots tend to
%have a number where Max = biggest disp in population (becuase is it part
%of a rising limb that is missed out).  Where this bin is in the circular
%mapping relative to the true disp changes the circular variance. And SD. 
for j = 1:length(dxs)
[resp(:,:,j), sd(j), r(j), p(:,j)] = OneDRds(nreps, varargin{:},'stimdisp',dxs(j));
end

function [resp, sd, r, p] = OneDRds(nreps,varargin)
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
        
colors = mycolors;
% Number of sets of images to generate
nimages = 3;
normalize = 0;
posonly = 0;
% seeds to use in generating each image
sffilter = 0;
showfilter = 0;
jf = 4; %default harmonic
jsf = 4;
dx = 0;
j = 1;
usephasedisp = 0;
symmetric = 1;

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
    elseif strncmpi(varargin{j},'rfphase',4)
        usephasedisp = 1;
    elseif strncmpi(varargin{j},'asymmetric',4)
        symmetric = 0;
    elseif strncmpi(varargin{j},'poso',4)
        posonly = 1;
    elseif strncmpi(varargin{j},'sf',2) %stimulus SF for Gabor, as harmonic
        j = j+1;
        jsf = varargin{j};
    elseif strncmpi(varargin{j},'stimdisp',7)
        j = j+1;
        dx = varargin{j};
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        figa = varargin{j};
    end
    j = j+1;
end
% Simple cell frequencies:
jf = 1;
freq = [4 6]/n;
nbw=4;
orbw=pi/5;
sf = freq(1);


%****************************************************************

phase = 0;
dphases = -pi:pi/16:pi;
disps = [-16:16]; %+- half period, for quivalanece with dp.
freqx = fftshift([-n/2:n/2-1])/n; 
x = [-n/2:n/2-1]; 
% ---- GENERATE BASIC IMAGE IN FOURIER DOMAIN
im = ((rand(n)>0.5)-0.5)*2; % binary random-dot pattern with no DC
FT = fft2(im);
GetFigure(figa);

    sx = sqrt(log(sqrt(2))) / pi ./ freq(jf) * ( 2.^nbw + 1) ./ ( 2.^nbw - 1) ;
            stimph = dx .* -2 .* pi .* freq(jf);
        if usephasedisp
            d = 0;
            for j = 1:length(dphases)
                GaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* cos( 2*pi.*(x-d).*freq(jf));
                OddGaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* sin( 2*pi.*(x-d).*freq(jf));
                dp = dphases(j);
                GaborR(j,:) = exp(-(x+d).^2./(2.*sx)^2) .* cos( 2*pi.*(x+d).*freq(jf)+dp);
                OddGaborR(j,:) = exp(-(x+d).^2./(2.*sx)^2) .* sin( 2*pi.*(x+d).*freq(jf)+dp);
            end
        elseif symmetric
            dphases = 2 .* pi .* disps .* freq(jf);
            for j = 1:length(disps)
%small offset breaks symmetry, so that mathces fall equally either side. 
                d = disps(j)/2 -0.01;
                GaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* cos( 2*pi.*(x-d).*freq(jf));
                OddGaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* sin( 2*pi.*(x-d).*freq(jf));
                GaborR(j,:) = exp(-(x+d).^2./(2.*sx)^2) .* cos( 2*pi.*(x+d).*freq(jf));
                OddGaborR(j,:) = exp(-(x+d).^2./(2.*sx)^2) .* sin( 2*pi.*(x+d).*freq(jf));
            end
        else
            dphases = 2 .* pi .* disps .* freq(jf);
            for j = 1:length(disps)
                d = disps(j);
            GaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* cos( 2*pi.*(x-d).*freq(jf));
            OddGaborL(j,:) = exp(-(x-d).^2./(2.*sx)^2) .* sin( 2*pi.*(x-d).*freq(jf));
            GaborR(j,:) = exp(-(x).^2./(2.*sx)^2) .* cos( 2*pi.*(x).*freq(jf));
            OddGaborR(j,:) = exp(-(x).^2./(2.*sx)^2) .* sin( 2*pi.*(x).*freq(jf));
        end
    end
disp = 0;
    hold off;
    if showfilter
subplot(2,1,1);
imagesc(Gabor);
subplot(2,1,2);
    else
subplot(1,1,1);
    end

    im = ((rand(n,nreps)>0.5)-0.5)*2; %nreps 1 d images
    if dx > 0
        imr(1+dx:size(im,1),:) = im(1:end-dx,:);
        imr(1:dx,:) = ((rand(dx,nreps)>0.5)-0.5)*2;
    else
        imr = im;
    end
    eimage = (GaborL*im + GaborR*imr).^2;
    oimage = (OddGaborL*im + OddGaborR*imr).^2;
    resp = eimage+oimage;

    [a,b] = max(resp);
    subplot(2,1,1);
    hist(b,nreps/20);
    subplot(2,1,2);
    imagesc(eimage+oimage);
    sd = std(b);
    r = abs(sum(cos(dphases(b))) + i*sum(sin(dphases(b))))./length(b);
    dps = sort(abs(dphases));
    perr = cos(stimph - dphases(b));
    for j = 1:length(dps)
        p(j) = sum(perr >= cos(dps(j)));
    end
xlabel('Disparity (pixels)')
ylabel('correlation');



%figify(gcf,gca);
%set(gca,'Xlim',[-90 90]);
