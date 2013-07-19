
if(~exist('len'))
  len = 1000;
end

clear ft;
for loop = 1:100
path = cumsum(binornd(1,0.5,1,len) - 0.5);
ft(loop,:) = abs(fft(path - mean(path)));
end
mft = mean(ft) / trapz(mean(ft));
plot(mft(1:100));
hold on;
