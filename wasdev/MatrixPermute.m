function X = MatrixPermute(M, order)
%for square matrix M, re-oder rows and columnts according to order

n = size(M,1);
for j = 1:n
    for k = 1:n;
        X(j,k) = M(order(j),order(k));
    end
end
