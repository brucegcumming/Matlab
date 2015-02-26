function [counts,rnd]= MatrixCorrCounts(r,ntrials,poolsize,varargin)
% [counts, rnd]=corr_counts(r,ntrials,poolsize,varargin)
% returns correlated spike counts, correlation matrix given by r, 
% for n trials (ntrials).
% default: mean firing of all the neurons is 50, fano factor: 1.5
% options: 'means',vector of size(1,poolsize), counts are calculated around
% these input means for each neuron in the pool.
% 'fano' fano-factor
% hn 08/22/06

j = 1;
fano=1.5;
means = ones(1,poolsize)*50;
while j  <= nargin-3
    str = varargin{j};
    if strcmp('means',str)
        means = varargin{j+1};
        j=j+1; 
    elseif strcmp('count',str)
        means = ones(1,poolsize).*varargin{j+1};
        j=j+1;         
    elseif strcmp('fourpools',str)
        r = ones(poolsize) .* 0.1;
    elseif strcmp('fano',str)
        fano = varargin{j+1};
        j=j+1;
    end   
  j= j+1;
end

counts=[];
N = poolsize;

if size(r,1) == 1
    Q = ones(N,N) .* r;
else
    Q = r;
end
for n = 1:N;
    Q(n,n) = 1;
end
Q = Q^0.5;


if fano > 0
    for n = 1:ntrials
        z = randn(N,1);
        rnd(:,n) = Q * z;
        %rnd(:,n) = Q * z;
        counts(:,n)=rnd(:,n).*sqrt(fano*means')+means';
    end
else %just calculate the normal rands
%    z = randn(N,1, ntrials);
    for n = 1:ntrials:1
        z = randn(N,1);
        rnd(:,n) = Q * z;
    end
    counts = rnd;
end
checkcorrs(rnd,1);
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

function cval = checkcorrs(a,offset)

%for each element in pool, calculate corr with next one
for j = 2:size(a,1)-offset
    x = corrcoef(a(j-1,:),a(j,:));
    corrs(j-1,1) = x(1,2);
end
cval = mean(corrs);
fprintf('Corrs %.3f\n',cval);
