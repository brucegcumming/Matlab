function [counts,rnd]=corr_counts(r,ntrials,poolsize,varargin)
% [counts, rnd]=corr_counts(r,ntrials,poolsize,varargin)
% based on Shadlen et al. 1996: returns correlated (mean correlation 
% in the pool = r) spike counts (m-by-n-matrix) of m (poolsize) neurons 
% in a pool for n trials (ntrials).
% default: mean firing of all the neurons is 50, fano factor: 1.5
% options: 'means',vector of size(1,poolsize), counts are calculated around
% these input means for each neuron in the pool.
% 'fano' fano-factor
% hn 08/22/06

j = 1;
fano=1.5;
CMatrix = [];

means = ones(1,poolsize)*50;
while j  <= nargin-3
    str = varargin{j};
    if strcmp('means',str)
        means = varargin{j+1};
        j=j+1; 
    elseif strcmp('count',str)
        means = ones(1,poolsize).*varargin{j+1};
        j=j+1;         
    elseif strcmp('fano',str)
        fano = varargin{j+1};
        j=j+1;
    end   
  j= j+1;
end

counts=[];
N = poolsize;
if length(CMatrix)
    Q = sqrtm(CMatrix);
else
if r == 0
    u = 1;
    v = 0;
else
    u = 1/(r*sqrt(N))*sqrt(2/N+ r-2*r/N - 2/N*sqrt((1-r)*(1-r+r*N)))* (1+sqrt((1-r)*(1-r+r*N)));
    v=1/sqrt(N)*sqrt(2/N+r-2*r/N-2/N*sqrt((1-r)*(1-r+r*N)));
end
Q = ones(N,N)*v;
for n = 1:N;
    Q(n,n) = u;
end
end
if fano > 0
    for n = 1:ntrials
        z = randn(N,1);
        if length(CMatrix)
            rnd(:,n) = Q * z;
        else
            rnd(:,n) = sum(v.*z) + (u-v) .* z;
        end
        counts(:,n)=rnd(:,n).*sqrt(fano*means')+means';
    end
else %just calculate the normal rands
%    z = randn(N,1, ntrials);
    for n = 1:ntrials:1
        z = randn(N,1);
        rnd(:,n) = sum(v.*z) + (u-v) .* z;
%        rnd(:,n) = sum(v.*z(:,:,n)) + (u-v) .* z(:,:,n);
        %rnd(:,n) = Q * z;
    end
    counts = rnd;
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