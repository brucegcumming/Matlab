function edgetuning(varargin)

%egdetuning(...)
%simple simulation of difference two 2-D Gaussians for simulationg
%sqcorrug tuning. 
%egdetuning('type',x) sets the type of plo:
%          0 is pseudocolor
%          1 is an orientation curve for each phase
%          2 is a phase resposne curve for each ori
%egdetuning('separation',x) gives a vector of separation values to evaluate
%egdetuning('offset',x) gives a vector of offsets (mis-centering)
%egdetuning('mode',x) uses 2 Gaussians if x == 0, two DOGs if x == 1
%
%

plottype = 0;
separations = 10;
offsets = [0 10 20];
gains = [1 1];

GAUSSIAN = 0;
DOG = 1;
mode = GAUSSIAN;

SUM = 0;
PRODUCT = 1;
SUMSQUARE = 2;
imode = SUM;

j = 1;
while j < nargin+1
    if strncmpi(varargin{j},'type',3)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'separations',3)
        j = j+1;
        separations = varargin{j};
    elseif strncmpi(varargin{j},'offsets',3)
        j = j+1;
        offsets = varargin{j};
    elseif strncmpi(varargin{j},'gains',3)
        j = j+1;
        gains = varargin{j};
    elseif strncmpi(varargin{j},'mode',3)
        j = j+1;
        modes = varargin{j};
        mode = modes(1);
        if length(modes) > 1
            imode = modes(2);
        end
    end
    j = j+1;
end


sd = [10 10];
fsize = [256 256];
nplot = 1;

lims = (fsize-1)/2;
[x, y] = meshgrid(-lims(1):lims(1),-lims(2):lims(2));


for separation = separations;
    for offset = offsets;

c = [separation/2 0 -separation/2 0] + [offset 0 offset 0]
g{1} = exp(-(((x-c(1))./sd(1)).^2 + ((y-c(2))./sd(2)).^2)/2);
g{2} = exp(-(((x-c(3))./sd(1)).^2 + ((y-c(4))./sd(2)).^2)/2);

if mode == GAUSSIAN
    ga = g{1} .* gains(1);
    gb = g{2} .* gains(2);
elseif mode == DOG
    ga = gains(1) .* (g{1} - g{2});
    gb = gains(2) .* (g{2} - g{1});
end
figtag = sprintf('EDGEFig%d',nplot);
GetFigure(figtag);
subplot(1,2,1);
imagesc(ga - gb);
axis('image');

ors = -pi:pi/10:pi;
steps = -30:5:30;
Z = [];
Za = [];
Zb = [];
for j = 1:length(ors);
  for k = 1:length(steps)
      im = edge(ors(j),'phase',steps(k));
      Za(k,j) = sum(sum(im .* ga));
      Zb(k,j) = sum(sum((~im) .* gb));
  end
end

if imode == SUM
    Z = Za+Zb;
elseif imode == SUMSQUARE
    Z = (Za + Zb).^2;
else
    Z = Za .* Zb;
end

[X,Y] = meshgrid(ors, steps);
subplot(1,2,2);
if plottype == 0
    pcolor(X,Y,Z);
elseif plottype == 1
    plot(steps,Z');
else
    plot(ors,Z);
end
title(sprintf('separation %d, offset %d',separation,offset));

nplot = nplot+1;
end
end