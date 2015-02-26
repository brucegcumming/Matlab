function C = OptimizeDprime(DATA)
global ncall;
global ycs;
ncall = 0;
ycs = [];
eid = DATA.currentexpt(1);
cl = DATA.currentcluster;
C = DATA.cluster{cl,DATA.probe};

byforce = 1;
%    options = optimset();
if byforce
x = [C.x(1) C.x(2) C.y(1) C.y(2)];
for k= 1:3
scales = [0.9:0.2:1.1];
n = 0;
for a = scales
for b = scales
for c = scales
for d = scales;
n = n+1;
xs(n,:) = x .* [a b c d];
end
end
end
end
for j = 1:n
DATA.cluster{1,DATA.probe}.x = [xs(j,1) xs(j,2) C.x(3)];
DATA.cluster{1,DATA.probe}.y = [xs(j,3) xs(j,4) C.x(3)];
if isfield(DATA,'ispk')
[DATA, dprimes(j), details] = SetSpkCodes(DATA,DATA.ispk,DATA.probe, 0);
else
[DATA, dprimes(j), details] = SetSpkCodes(DATA,DATA.spklist,DATA.probe, 0);
end
[dprime, best] = max(dprimes);
x = xs(best,:);
end
end
else 

%for some reason this seems to get suck in very local minima. I
%gueess because dprime is not a smoth function of the variables on
%a fine scale. see lemM0600 probe 4 expt 39
options = optimset('TolFun',0.1,'TolX',0.01,'Largescale','off','simplex','on');
[x, dp, flag, output] = fminsearch(@cmb.ClusterDprime, [C.x(1) C.x(2) C.y(1) C.y(2)],options, DATA, C.nspk);
%    [x, dp, flag] = fminsearch(@cmb.ClusterDprime, x,options, DATA, C.nspk);
dprime = -dp;
end
C.x = [x(1) x(2) 1];
C.y = [x(3) x(4) 1];
C.dprime = dprime;



