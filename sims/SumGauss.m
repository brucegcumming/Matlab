function SumGauss(A,u,sd)

x = -100:100;

for j = 1:length(A)
    G(j,:) = Gauss([u(j) sd(j) A(j)],x);
end
plot(sum(G,1));