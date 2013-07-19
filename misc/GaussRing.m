function Z = GaussRing(npix,r, sd)

[x,y] = meshgrid(linspace(-1,1,npix),linspace(-1,1,npix));

dr = abs((x.^2 + y.^2)-(r.^2)); 
Z = exp(-dr/(2 * sd *sd));
