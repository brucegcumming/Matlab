function im = edge(angle, varargin)

isize = [256 256];
j = 1;
plotedge = 0;
ostep = 0;
while j < nargin
    if(strncmpi(varargin{j},'size',4))
        isize(1) = varargin(j+1);
        isize(2) = varargin(j+2);
        j = j+2;
    elseif(strncmpi(varargin{j},'phase',4))
        j = j+1;
        ostep = varargin{j};
    elseif(strncmpi(varargin{j},'plot',4))
        plotedge = 1;
    end
    j = j+1;
end

im = zeros(isize);
theta = angle;
off(1) = ostep * cos(theta);
off(2) = ostep * sin(theta);

lims = (isize-1)/2;
[x, y] = meshgrid(-lims(1):lims(1),-lims(2):lims(2));
z = (x+off(1)) .* cos(theta) + (y + off(2)) .* sin(theta);

idx = find(z >= 0);
im(idx) = 1;

if plotedge
    imagesc(im);
    axis('image');
end
