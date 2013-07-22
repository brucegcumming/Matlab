function mpwr = rdsft(nrpt)

w=401;
h=401;

for j = 1:nrpt;
    pwr(:,:,j) = ifft2(fft2(round(rand(w,h))));
end
mpwr = mean(pwr,3);
mpwr(1,1) = 0;
imagesc(mpwr);
    