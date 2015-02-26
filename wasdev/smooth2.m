function Z = smooth2i(x,y, z, X, Y, sd, varargin);


if min(size(x)) == 1
    [x,y] = meshgrid(x,y);
end
if min(size(X)) == 1
    [X,Y] = meshgrid(X,Y);
end
Z = parzenSurf2d(x,y,z,X,Y,sd);
        