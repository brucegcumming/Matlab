function r = Consist(nsd)




nt = 10000;
for j = 1:length(nsd)
signal = randn(nt,1);
noise = randn(nt,2) .* nsd(j);

a = corrcoef(signal, signal+noise(:,1));
r(1,j) = a(1,2);
pc(1,j) = sum(sign(signal) == sign(signal+noise(:,1)))./nt;
pc(2,j) = sum(sign(signal+noise(:,2)) == sign(signal+noise(:,1)))./nt;

a = corrcoef(signal+noise(:,2), signal+noise(:,1));
r(2,j) = a(1,2);
end

if j > 1
    plot(pc(1,:),pc(2,:));
end