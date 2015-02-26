function [counts, bcounts, details]=twopoolcounts(r,rb,ntrials,poolsize,varargin)
% [acounts, bcounts, details]=twopoolcounts(r,rb,ntrials,poolsize,varargin)
% based on Shadlen et al. 1996: returns correlated (mean correlation 
% in the pool = r) spike counts (m-by-n-matrix) of m (poolsize) neurons 
% in a pool for n trials (ntrials).
% default: mean firing of all the neurons is 50, fano factor: 1.5
% options: 'means',vector of size(1,poolsize), counts are calculated around
% these input means for each neuron in the pool.
% 'fano' fano-factor
% hn 08/22/06

check = 0;
j = 1;
fano=1.5;
means = ones(1,poolsize)*50;

CMatrix = [];
mx = round((r - rb) * poolsize/r); %%used to be m. Not sure what its for
details.corr = [r r * (poolsize-mx)/poolsize];
m = 0;
while j  <= nargin-4
    str = varargin{j};
    if strcmp('means',str)
        means = varargin{j+1};
        j=j+1; 
    elseif strncmpi('off',str,3)
        m = varargin{j+1};
        j=j+1;
    elseif strncmpi('check',str,3)
        check = 1;
    elseif strncmpi('cmatrix',str,6)
        j = j+1;
        CMatrix = varargin{j};
    elseif strcmp('fano',str)
        fano = varargin{j+1};
        j=j+1;
    end   
  j= j+1;
end


N = poolsize;
counts=[];
oldway = 0;
if length(CMatrix)
Q = sqrtm(CMatrix);
vv(1) = Q(1,1);
vv(2) = Q(1,2);
vv(3) = Q(1,end);
elseif oldway
vv = Finduvw(N,r,rb);
Q = ones(2* N, 2* N);
Q = Q * vv(2) + diag(diag(Q)) .* (vv(1)-vv(2));
Q(N+1:end,1:N) = vv(3);
Q(1:N,N+1:end) = vv(3);
else
CM = ones(2.*N) .* r;
CM(N+1:end,1:N) = rb;
CM(1:N,N+1:end) = rb;
for j = 1:2*N
    CM(j,j) = 1;
end
Q = sqrtm(CM);
vv(1) = Q(1,1);
vv(2) = Q(1,2);
vv(3) = Q(1,end);
end



if check == 2
    for n = ntrials:-1:1
        z = randn(N,1);
        counts(:,n)=Q*z;
        bcounts(:,n)= sum(v.*z) + (u-v) .* z;;
    end
    details.diffs = counts - bcounts;
    return
end

for n = ntrials:-1:1
     z = randn(2 * N,1);
%    rnd = sum(v*z) + (u-v) .* z;
%    arnd = rnd(1:poolsize) + sum(w*z(1:poolsize));
%    brnd = rnd(poolsize+1:end) + sum(w*z(poolsize+1:end));
    rnd(:,n) = Q*z;
    counts(:,n) = rnd(1:N,n).*sqrt(fano*means')+means';
    bcounts(:,n)= rnd(1+N:end,n).*sqrt(fano*means')+means';
end

if m >=1
    details.offset = m;
    else
    details.offset = 0;
    m = 0;
end
details.uvw = [vv(1) vv(2) vv(3)];
if check
 details.meancorr = checkcorrs(counts,bcounts,details.offset);
end

%
%  Q = a b b b  * z = [w x y z]  = wa+xb+yb+zb
%      b a b b                     wb+xa+yb+zb
%      b b a b                     wb+xb+ya+zb
%      b b b a                     wb+xb+yb+za
%
% so every pair of numbers has (n-2) random vector* b in common.
%
%
%  Q = a b c c  * z = [w x y z]  = wa+xb+yb+zb
%      b a c c                     wb+xa+yb+zb
%      b b d c                     wc+xc+yd+zc
%      b b c d                     wc+xc+yc+zd
%
% to do pools, make z n+m elements long. make pool 1 z(1:n) * Q, pool2
% z(m:end) * Q;
% but then element 1 of pool2 is highly correlated with element m of pool1
%  Q = 0 b b b  * z = [w x y z]  = xb+yb+zb
%      b 0 b b                     wb+yb+zb 
%      b b 0 b                     wb+xb+zb
%      b b b 0                     wb+xb+yb
%
%   b b b b  * z = [w x y z] + (a-b) * [w x y z] is the same.
%   b b b b                      
%   b b b b                     
%   b b b b   
% but key may be that b *z is weakly correlated with w,x,y,z, not indep. 
% one possibility is to calulate a pool of 2N at the lower correlation,
% then add to each individual pool a signal derived from upper/lower
% corners. 
%  a b c c  * z = [w x y z]  = wa+xb+yb+zb
%  b a c c                     wb+xa+yb+zb
%  b b d c                     wc+xc+yd+zc
%  b b c d                     wc+xc+yc+zd
% Then add to [w x]
%  d d  * [w x]  
%  d d
% And add ot [y z]
%  d d  * [y y]  
%  d d

function cval = checkcorrs(a,b,offset)

%for each element in pool, calculate corr with next one
for j = 2:size(a,1)-offset
    x = corrcoef(a(j-1,:),a(j,:));
    corrs(j-1,1) = x(1,2);
    x = corrcoef(b(j-1,:),b(j,:));
    corrs(j-1,2) = x(1,2);
    x = corrcoef(b(j,:),a(j+offset,:));
    corrs(j-1,3) = x(1,2);
end
cval = mean(corrs);
%fprintf('Corrs %.3f\n',cval);


function [x, check] = Finduvw(N,r,s)

options = optimset('TolFun',1e-12);
x = [1 0 0];
x = fminsearch(@evaluvw, x, options, N,r,s);

dd = (2 * x(1) + 2 * (N-1) * x(2));
x(3) = s/dd;


function sumsq = evaluvw(x,N,r,s)


dd = (2 * x(1) + 2 * (N-1) * x(2));
Ns =  N * s ^2/dd^2;  %%N * c^2
err(1) = x(1).^2 + (N-1) * x(2)^2 + Ns -1;
err(2) = 2 * x(1) * x(2) + (N-2) * x(2)^2+ Ns - r;

sumsq = sum(err.^2);

