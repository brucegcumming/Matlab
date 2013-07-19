function Z = interpf(x, y, z, X, Y, interptype, smoothParzen)
% Z = interpf(x, y, z, X, Y, type, smooth)
% interpolate z = f(x,y) on X, Y grids

if (interptype == 0)
	interpmode = 'Linear';
	Z = griddata(x, y, z, X, Y, interpmode);
elseif (interptype == 1)
	interpmode = 'Gaussian';
	[rows, cols] =  size(x);
	if(cols > 1)
        Z = parzenSurf2d(x, y, z, X, Y, smoothParzen);
    else
        Z = parzenSurf(x, y, z, X, Y, smoothParzen);
    end
end

