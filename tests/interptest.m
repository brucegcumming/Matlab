X = 0:20;
X = repmat(X,21,1);
Y = X';
Z = ones(size(Y));
pcolor(X,Y,Z);
[xi, yi] = meshgrid(linspace(min(min(X)),max(max(X))),linspace(min(min(Y)),max(max(Y))));
zi = Interpf(X,Y,Z,xi,yi,1,1);
pcolor(xi,yi,zi);
colorbar;


X = X .*10;
Y = Y .* 10;

[xi, yi] = meshgrid(linspace(min(min(X)),max(max(X))),linspace(min(min(Y)),max(max(Y))));
zi = Interpf(X,Y,Z,xi,yi,1,1);
mean(mean(zi))
zi = Interpf(X,Y,Z,xi,yi,1,10);
mean(mean(zi))
pcolor(xi,yi,zi);
colorbar;
