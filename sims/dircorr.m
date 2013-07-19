function dircor(varargin)



[x,y] = meshgrid([0:0.01:2*pi],[0:0.01:2*pi]);
c = cos(x-y+pi/2);
imagesc(c>0);

