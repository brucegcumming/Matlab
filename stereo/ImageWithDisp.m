function ImageWithDisp(image,disp)

n = size(image,1);
offset = min(min(image));
scale = max(max(image)) - offset;
rimage = [zeros(n,disp) image];
limage = [image zeros(n,disp)];
cimage(:,:,1) = (rimage - offset)./scale;
cimage(:,:,2) = (limage - offset)./scale;
cimage(:,:,3) = (limage - offset) ./scale;
imagesc(cimage);
