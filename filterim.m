function [nbuilt, details] = filterim( sf, or, wi, varargin)

%filterim( sf, or, rsd, varargin)
% make SF/ori bandpass filtered image
% sf(1) is peak, sf(2) is SD of Gaussian in frequency domain (cpd)
% or(1) is peak orientaiton, or(2) is Gaussian sd.
% rsd is a Gaussian envelope applied to the result.
%
%filterim( sf, or, rsd, 'save', dir) saved PGM images in dir
%
%filterim( sf, or, rsd, 'nseed', n), builds n images with different seeds.
%
%filterim( sf, or, rsd, 'imzise', s) maxes the images s x s pixels. [256]
%
%filterim( sf, or, rsd, 'deg2pix',x ) sets up the pixel size. 


savedir = [];
gamma = 0;
nseed = 1;
imsize = 256;
pix2deg = 0.0292;
seedoffset = 0;
replace = 1;
checkreplace = 0;
version = '$Revision: 1.4 $';
version = version(11:end-1); % just the numbber
[err, host] = unix('hostname');
quietmode = 0;
nuilt = 0;
envelopetype = 'Gauss';
envelopeid = 0;
rfenvelope = [];
envsmw = 1; %widht
parallel = 0;

savedir = [];
getft = 0;
state.plotim =1;
j = 1;
while j < nargin -2
    if(strncmpi(varargin{j},'save',4))
        j = j+1;
        savedir = varargin{j};
	mkpath([savedir '/']);
    elseif(strncmpi(varargin{j},'checkreplace',4))
        checkreplace = 1;
    elseif(strncmpi(varargin{j},'rfenvelope',4))
        j = j+1;
        rfenvelope = varargin{j};
    elseif(strncmpi(varargin{j},'envelope',4))
        j = j+1;
        envelopetype = varargin{j};
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            envsmw = varargin{j};
        end

    elseif(strncmpi(varargin{j},'getft',4))
        getft = 1;
    elseif(strncmpi(varargin{j},'nseed',4))
        j = j+1;
        nseed = varargin{j};
    elseif(strncmpi(varargin{j},'noplot',4))
        state.plotim = 0;
    elseif(strncmpi(varargin{j},'imsize',4))
        j = j+1;
        imsize = varargin{j};
    elseif(strncmpi(varargin{j},'deg2pix',4))
        j = j+1;
        pix2deg = varargin{j};
    elseif(strncmpi(varargin{j},'pix2deg',4))
        j = j+1;
        pix2deg = varargin{j};
    elseif(strncmpi(varargin{j},'parallel',4))
        parallel = 1;
    elseif(strncmpi(varargin{j},'quieter',6))
        quietmode = 2;
    elseif(strncmpi(varargin{j},'quiet',4))
        quietmode = 1;
    elseif(strncmpi(varargin{j},'gamma',4))
        j = j+1;
        gamma = varargin{j};
    elseif(strncmpi(varargin{j},'noreplace',4))
        replace = 0;
    elseif(strncmpi(varargin{j},'seed',4))
        j = j+1;
        seedoffset = varargin{j};
    end
    j = j+1;
end



if checkreplace
    propfile = [savedir '/paramlist'];
    binocfile = [savedir '/lastused'];
    a = dir(propfile);
    b = dir(binocfile);
    if ~isempty(a) & ~isempty(b) & b.datenum > a.datenum
        replace = 1;
        fprintf('Replacing %s\n',savedir);
    else
        replace = 0;
    end
end


%Equivalent width of Gaussian is 2.5066 sigma
rsd = wi/(2.5066 * pix2deg);
center = 1+ (imsize/2);
[x,y] = meshgrid(1:imsize,1:imsize);
r = sqrt((x-center).^2 + (y-center).^2);
envelopeid = find(strcmp(envelopetype,{'Gauss' 'cosine'}));
if strcmp(envelopetype,'cosine')
    edger = ((wi/pix2deg)-r) * pi*pix2deg/envsmw;
    envelope = (1+sin(edger))/2;
    envelope(edger < -pi/2) = 0;
    envelope(edger > pi/2) = 1;
else
    envelope = exp(-(r.^2)/(2*rsd^2));
end
if ~isempty(rfenvelope)
    envelope = envelope .* rfenvelope;
end
angle = atan2(y-center,x-center);

orpeak = (90 - or(1)) * pi/180;
dangle = atan(tan(angle-orpeak));
dangleb = dangle + pi;
danglea = dangle -pi;
danglec = dangle +pi*2;
dangled = dangle -pi*2;
orsd= or(2) * pi/180;
sfsd = sf(2) * imsize * pix2deg;
sfpeak = sf(1) * imsize * pix2deg;
if(orsd > 0.7* pi)  %uniform ori is sd > 126 degrees (monkey can do 90!!)
  filter = exp(-((r-sfpeak).^2)./(2 * sfsd^2));
elseif sf(2) >= 100 | isinf(sfpeak)  %All sfs  is sfsd > 100, before convert to pixels
  orf = exp(-dangle.^2/(2 * orsd^2));
  orfa = exp(-danglea.^2/(2 * orsd^2));
  orfb = exp(-dangleb.^2/(2 * orsd^2));
  orfc = exp(-danglec.^2/(2 * orsd^2));
  orfd = exp(-dangled.^2/(2 * orsd^2));
  orfall = (orf+orfb+orfa+orfc+orfd)./3;
  orfall = orfall ./ max(orfall(:));
  id = (find(r < 0.3));
  filter = orfall ./(abs(r));
  filter(id) = 0;
  filter = orfall;
else
  orf = exp(-dangle.^2/(2 * orsd^2));
  orfa = exp(-danglea.^2/(2 * orsd^2));
  orfb = exp(-dangleb.^2/(2 * orsd^2));
  orfc = exp(-danglec.^2/(2 * orsd^2));
  orfd = exp(-dangled.^2/(2 * orsd^2));
  orfall = (orf+orfb+orfa+orfc+orfd)./3;
  orfall = orfall ./ max(orfall(:));
  filter = exp(-((r-sfpeak).^2)./(2 * sfsd^2)) .* orfall;
end

% for some reason FTa .*FTb does not quite convolve, need FTa .*
% fftshift(FTb), so do this by just NOT shifting the fitler....
details.filter = filter;
filter = fftshift(filter);
    tic;
    nbuilt = 0;
    imsum = zeros(256);
    fts = zeros([nseed 256 256]);

    
    
    
    
    
    state.seedoffset = seedoffset;
    state.quietmode = quietmode;
    state.gamma = gamma;
    state.buildim = 1;
    state.replace = replace;
    if parallel
        parfor seed = (1:nseed)
            im = BuildIm(savedir,seed, state, filter, envelope);
            if ~isempty(im)
            nbuilt = nbuilt+1;
            end
        end
    else
        for seed = 1:nseed
            [im, ft] = BuildIm(savedir,seed, state, filter, envelope);
            if ~isempty(im)
            nbuilt = nbuilt+1;
            if getft
                fts(nbuilt,:,:) = fft2(im-0.5);
            end
            end
        end
    end
    
details.image = im;
details.ft = ft;

if ~isempty(savedir) & nbuilt
    fid = fopen(sprintf('%s/paramlist',savedir),'w');
    if fid > 0
        fprintf(fid,'Build at %s on %s\n',datestr(now),host);
        fprintf(fid,'Version %s\n',version);
        fprintf(fid,'SeedOffset %d\n',seedoffset);
        fprintf(fid,'pix2deg %.4f\n',pix2deg);
        fprintf(fid,'sf %.4f sd %.4f\n',sf(1),sf(2));
        fprintf(fid,'or %.2f sd %.2f\n',or(1),or(2));
        fprintf(fid,'gamma %.2f\n',gamma);
        fprintf(fid,'imsz %.0f sd %.2f\n',imsize,wi/2.5066);
        fprintf(fid,'envelope %d %.1f\n',envelopeid,envsmw);
        fclose(fid);
    end
    fid = fopen(sprintf('%s/paramhistory',savedir),'a');
    if fid > 0
        fprintf(fid,'%s V%s se%d px %.4f',datestr(now),version,seedoffset,pix2deg);
        fprintf(fid,'sf %.4f %.4f or%.1f env%d',sf(1),sf(2),or(1),or(2), envelopeid);
        fclose(fid);
    end
end
if nbuilt >1 | quietmode < 2
fprintf('Took %.3f for seed %d %s (%d new images)\n',toc,seedoffset,savedir,nbuilt);
end

if getft
    nbuilt = fts;
end


function [im, ft] = BuildIm(savedir, seed, state,filter, envelope)
    if ~isempty(savedir)
        imname = sprintf('%s/se%d.pgm',savedir,seed);
        if exist(imname,'file') && ~state.replace
            if state.quietmode < 1
                fprintf('%s exists\n',imname);
            end
            state.buildim = 0;
        end
    end
    if state.buildim
        rand('state',seed+state.seedoffset);
        im = rand(256) - 0.5;
        ft = fft2(im);
        
        ft = fftshift(ft).* filter;
        im = real(ifft2(ft)) .* envelope;
 %       im = real(ifft2(ft));
        scale = 0.5/max(max(abs(im)));
        im = 0.5 + im .* scale;
        if state.gamma
            im = im .^(1/state.gamma);
        end
        if ~isempty(savedir)
            imwrite(im,sprintf('%s/se%d.pgm',savedir,seed),'PGM');
        elseif state.plotim
            GetFigure('Filtered Images');
            subplot(1,2,1);
            imagesc(abs(fftshift(ft)));
            subplot(1,2,2);
            imagesc(real(im));
            colormap('gray');
        end
    else
        im = [];
        ft = [];
    end
