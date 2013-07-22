hold off;
pix2deg = 0.025;
x = 1:1024;
lf = 1;
x = pix2deg *(x-length(x)/2);
y = cos(x .* 2 * pi .* lf);
%plot(x,y);
z = (1+cos(x .* 2 * pi .* lf/8))/2;
plot(x,z.*y);
hold on;

lrfs(1,:) = Gabor(1,'pos',-4+pd,'npts',length(x),'pix',pix2deg);
lrfs(2,:) = Gabor(1,'pos',pd,'npts',length(x),'pix',pix2deg);
lrfs(3,:) = Gabor(1,'pos',4+pd,'npts',length(x),'pix',pix2deg);
rrfs(1,:) = Gabor(1,'pos',-4-pd,'npts',length(x),'pix',pix2deg);
rrfs(2,:) = Gabor(1,'pos',-pd,'npts',length(x),'pix',pix2deg);
rrfs(3,:) = Gabor(1,'pos',4-pd,'npts',length(x),'pix',pix2deg);


plot(x,rfs(1,:),'r');
plot(x,rfs(2,:),'r');


edisps = -2:0.1:2
%edisps = 0;
resps = [];
br = [];
weights = [-0.5 1 -0.5];
for n = 1:length(edisps)
    zr = (1+cos((x-edisps(n)) .* 2 * pi .* lf/8))/2;
    zl = (1+cos((x+edisps(n)) .* 2 * pi .* lf/8))/2;
    for j = 1:size(rrfs,1)
      resps(:,j) = [sum(rrfs(j,:).*zr.*y) sum(lrfs(j,:).*zl.*y)];
    end
    br(n) = sum(sum(resps) .* weights);
end

hold off;
plot(edisps,br);
    