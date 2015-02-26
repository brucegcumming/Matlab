function res = FitGaussMeans(X,N, varargin)
verbose = 0;
dprime = 0;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'verbose',4)
verbose = 1;
end
j = j+1;
end

G = gmdistribution.fit(X,N,'Options',statset('MaxIter',1000));
res.obj = G;
distance = mahal(G,G.mu);
distance = sqrt(2./((1./distance(1,2))+(1./distance(2,1))));
sigmas = diag(mean(G.Sigma,3)); % sum variances across componets, in each dimenions

if N == 2
nsd = diff(G.mu)./sqrt(sigmas)';
dprime = sqrt(sum(nsd.^2));
res.dprime = dprime;
end
res.mahal = distance;
if verbose
fprintf('Distance %.2f (%.2f)\n',distance,dprime);
end



