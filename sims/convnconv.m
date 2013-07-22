function convnconv(n)

k = zeros(1,20);
k(5:15) = 1;

for j = 1:n
    newk = conv(k,k);
    newk = newk./max(newk);
    plot(newk);
    hold on;
    k = newk;
end
