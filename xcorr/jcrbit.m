[x2d,y2d]=meshgrid(x,y);
for j=1:ntheta
xprime = x2d.*  cos(theta) - y2d.*sin(theta);
grating = sin(2*    pi*freq.*xprime);
    
FT(jth) = sum(sum(grating.*RF));
end