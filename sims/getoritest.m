function [a,b] = getoritest(varargin)

j = 1;
nsd = 0;
f = 0.2;
npts = 100;
while j <= length(varargin)
    if strncmpi(varargin{j},'noise',3)
        j = j+1;
        nsd = varargin{j};
    elseif strncmpi(varargin{j},'npts',3)
        j = j+1;
        npts = varargin{j};
    elseif strncmpi(varargin{j},'f',1)
        j = j+1;
        f = varargin{j};
    end
    j = j+1;
end
x = (rand(npts,1)-0.5).*2;
y = (rand(size(x))-0.5).*2;
X = x;
Y = y;
noise = randn(size(x)) .* nsd;
z = cos(2 * pi * f * y) + noise;
[a,b] = getori(x,y,z-mean(z(:)));
subplot(3,1,1);
imagesc(a);
subplot(3,1,2);
imagesc(b);
[xi,yi] = meshgrid([-1:0.1:1],[-1:0.1:1]);
Z = interpf(x,y,z,xi,yi,1,0.2);
subplot(3,1,3);
imagesc(Z);

%plot3(x(:),y(:),z(:),z(:),'o');