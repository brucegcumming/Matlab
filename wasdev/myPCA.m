function [E,V, pcs] = myPCA(X, varargin)
%Do quick PCA on matrix, plot clustersing;
nc = 4;
noplot = 0;
normalize = 0;

j=1;
while j <= length(varargin)
    if strncmpi(varargin{j},'normalize',3)
        normalize =1;
    end
    j = j+1;
end

if normalize
    X = myNormalize(X,'std');
end

[E,V] = eig(cov(X));
V = diag(V);
if nc > length(E)
    nc = length(E);
end
for j = 1:nc
pcs(:,j) = X * E(:,end+1-j);
end
if noplot == 0
    PlotND(pcs,[],varargin{:});
end