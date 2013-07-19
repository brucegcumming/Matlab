function PlotPGM(name,varargin);


im = imread(name, 'PGM');
imagesc(im);
colormap('gray');