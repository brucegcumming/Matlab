function [x,y,z] = GaborNoise(sd,varargin)

gb = [1 0.1 (pi+0.8) 0.8 0 0 0 0.2];
j = 1;
while j < nargin
    if strncmpi(varargin{j},'gabor',3)
        j = j+1;
        gb = varargin{j};
    end
    j = j+1;
end
x = Gabor(gb,'npts',16);
z = abs(fftshift(fft2(x)));
x = x + randn(size(x)) .* sd;
y = abs(fftshift(fft2(x)));

subplot(2,1,1);
cmax = max(abs(x(:)));
imagesc(x); caxis([-cmax cmax]);
colormap('hot');
axis('image');
colorbar;
subplot(2,1,2);
hold off;

plot(y(10:end,9),'o-');
hold on;
plot(z(10:end,9),'ro-');
%plot(y(8,:),'r');

