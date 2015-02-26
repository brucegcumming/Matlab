function dprime = ClusterDprime(x, DATA, nspks)
global dprimes;
global ncall;
global ycs;

eid = DATA.currentexpt(1);
C = DATA.cluster{DATA.currentcluster,DATA.probe};
if length(x) == 2
C.x(2) = x(1);
C.y(2) = x(2);

else
C.x(1) = x(1);
C.x(2) = x(2);
C.y(1) = x(3);
C.y(2) = x(4);
end
ispk = DATA.ispk;

DATA.cluster{1,DATA.probe} = C;
[DATA, dprime, details] = SetSpkCodes(DATA,ispk,DATA.probe, 0);

ns = details.nc(1);
ncall = ncall+1;
dprimes(ncall) = dprime;
ycs(ncall,:) = x;
if ns < nspks * 0.5 | ns > nspks * 1.5 | C.x(1) > details.maxx | C.y(1) > details.maxy 
dprime = 1000;
else
dprime = -dprime;
end




