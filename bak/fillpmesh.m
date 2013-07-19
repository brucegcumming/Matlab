function [X, Y, Z] = fillpmesh(x, y, z, varargin)

%[X , Y, Z] = fillpmesh(x,y,z)
%adds an extra row and column of values so that pcolor(X,Y,Z), shows all
%the data(!!) even with no shading.

plotting = 0;

j = 1;
while j < nargin - 2
    if strncmpi(varargin{j},'plot',4)
        plotting = 1;
    elseif strncmpi(varargin{j},'smplot',4)
        plotting = 2;
        j = j+1;
        sigma = varargin{j};
    end
    j = j+1;
end
if(min(size(x)) == 1)
  [X, Y] = meshgrid(x,y);
else
  X = x;
  Y = y;
end
if size(X,1) > 1
    dx = X(size(X,1),1) -  X(size(X,1)-1,1);
else
    dx = 1;
end
X = [X; X(size(X,1),:)+dx];

if size(X,2) > 1
    dx = X(1,size(X,2)) -  X(1,size(X,2)-1);
else
    dx = 1;
end
X = [X X(:,size(X,2))+dx];

if size(Y,1) > 1
    dy = Y(size(Y,1),1) -  Y(size(Y,1)-1,1);
else
    dy = 1;
end
Y = [Y; Y(size(Y,1),:) + dy];


if size(Y,2) > 1
    dy = Y(1,size(Y,2)) -  Y(1,size(Y,2)-1);
else
    dy = 1;
end
Y = [Y Y(:,size(Y,2))+dy];

Z = [z; z(size(z,1),:)];
Z = [Z Z(:,size(Z,2))];

if plotting == 2
    xi = linspace(min(X(:)),max(X(:)));
    yi = linspace(min(Y(:)),max(Y(:)));
    [lx,ly] = meshgrid(xi,yi);
    lz = interpf(X,Y,Z,lx,ly,1,sigma);
    pcolor(lx,ly,lz);
    shading('flat');
    
elseif plotting
    pcolor(X,Y,Z);
    shading('flat');
end
return;
