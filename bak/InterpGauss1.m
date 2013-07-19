function Y = InterpGauss1(x,y,X,sd)

step = range(X)./length(X); %mean step size;
v = sd.^2;
for j = 1:length(X)
    d = (X(j) - x).^2;
    G = exp(-2.*d./v);
    Y(j) = sum(G .* y)./sd;
end
