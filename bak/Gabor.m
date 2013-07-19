function [g, x, y] = Gabor(params, varargin)

%
% Gabor(params, ...) generates a Gabor function, in one or 2 dimensions
%   params is a vector containing
%     [freq sd phase amplitude xmean baseline ori sdy]
%   If the length of params is shorter than 6, default parameters
%   are filled in.
%
%   A 2-dimensional Gabor is generated if the length of params is
%   longer than 6, or any parameters are set (orientation, SDy)
%   which imply a d-D function
%
% units are in degrees, each pixel = 0.0292 degreees, default width 256
% pixels


dims = 1;
pix2deg = 0.1;
imsize = 256;
ori = 0;
f = params(1);
npix = 256;
ym = 0;
x = [];

j = 1;
while j < nargin
    if isnumeric(varargin{j})
        if (length(varargin{j}) == 1)
            imsize = varargin{j};
        end
    elseif strncmpi(varargin{j},'pos',3) | strncmpi(varargin{j},'xpos',3)
        j = j+1;
        xm = varargin{j};
    elseif strncmpi(varargin{j},'npts',3)
        j = j+1;
        imsize = varargin{j};
    elseif strncmpi(varargin{j},'pix',3)
        j = j+1;
        pix2deg = varargin{j};
    elseif strncmpi(varargin{j},'xv',2)
        j = j+1;
        x = varargin{j};
    end
    j = j+1;
end

if length(params) < 2
  sx = 0.5/f;
else
  sx = params(2);
end

if length(params) < 3
  phi = 0;
else
  phi = params(3);
end
if length(params) < 4
  A = 1;
else
  A = params(4);
end

if exist('xm');
elseif length(params) < 5 
  xm = 0;
else
  xm = params(5);
end

if length(params) < 6
  base = 0;
else
  base = params(6);
end


if length(params)> 6
    dims = 2;
  ori = params(7);
  if length(params) > 7
    sy = params(8);
  else
    sy = sx; 
  end
end
npts = imsize;

if isempty(x)
x = pix2deg .* (ceil([1:npts] - npts/2) -1);
end
y = pix2deg .* (ceil([1:npts] - npts/2) -1);

if dims == 1
g = base + A .* exp(-(x-xm).^2/(2.*sx.^2)) .* cos(2 * pi * f .* (x-xm) + phi);
else
    [x2d,y2d] = meshgrid(x,y);
    x2 = x2d .* cos(ori) + y2d .* sin(ori);
    y2 = y2d.*cos(ori) - x2d .* sin(ori);
    g = base + A .* exp(-(x2-xm).^2./(2.*sx.^2)) .* cos(2 * pi * f * (x2-xm) + phi) ...
        .* exp(-(y2-ym).^2./(2.*sy.^2));
end
