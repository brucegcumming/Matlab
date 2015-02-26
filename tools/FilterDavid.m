function X = FilterDavid(X,smoothw)

ts = now;
for j = 1:size(X,1)
    sm = smooth(double(X(j,:)),smoothw);
    V = double(X(j,:)) - sm;
    m = max(abs(V));
    X(j,:) = int16(round(V .* 32000./m));
end
fprintf('Took %.2f\n',mytoc(ts));

