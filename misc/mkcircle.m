az = rand(100,1) .* pi * 2;
el = zeros(size(az));
r = ones(size(az));

[x,y,z] = sph2cart(az,el,r);
x2d = [x(:),y(:),z(:), ones(size(x))];
A = viewmtx(0,45);
xr = A*x2d';
A = viewmtx(45,0);
xa = A*xr;
plot3(xa(1,:),xa(2,:),xa(3,:),'o');