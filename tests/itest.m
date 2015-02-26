
[x, y ] = meshgrid([1 2 3 4 5 6 7 8 9 10], [1 2 3 4 5 6 7 8 9 10]);
z = x + y;
subplot(2,1,1);
pcolor(x,y,z);
[xi, yi] = meshgrid(linspace(1,5,40), linspace(1,5,40));
zi = Interpf(x,y,z,xi,yi,1,1);
subplot(2,1,2);
pcolor(xi,yi,zi);
  shading('interp');

