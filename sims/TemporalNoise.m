function TemporalNoise(varargin)

npts = 1000;
nev = 10000;
len = 50;
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'nevents',3)
        j = j+1;
        nev = varargin{j};
    end
    j = j+1;
end

t = -len:len;
y = exp(-abs(t));

V = zeros(1,npts);
r = 50+ceil(rand(1,nev).*(npts-100));
for j = 1:length(r)
    idx = r(j)-50:r(j)+50;
    V(idx) = V(idx)+y;
end
V = V(len:end-len);
plot(V);
title(sprintf('%d events',nev));