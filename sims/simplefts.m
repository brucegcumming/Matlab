function simplefts(varargin)

nf = 2;
w = 256;
bw = 2;
j = 1;
group = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'bw',2)
        j = j+1;
        bw = varargin{j};
    elseif strncmpi(varargin{j},'group',4)
        j = j+1;
        group = varargin{j};
    end
    j = j+1;
end

if group == 1
    nf = 4;
    im = zeros(w,1);
    im(round(w/2)) = 1;
    im = im - mean(im);
    showft(im,nf, 1);

    im = zeros(w,1);
    im(round(w/2)-1:round(w/2)) = 1;
    im = im - mean(im);
    showft(im,nf, 2);

    im = zeros(w,1);
    im(round(w/2)-1:round(w/2)+2) = 1;
    im = im - mean(im);
    showft(im,nf, 3);

    im = zeros(w,1);
    im(1:round(w/2)) = 1;
    im = im - mean(im);
    showft(im,nf, 4);
else    
    
    
im = zeros(w,1);
im(round(w/2)-bw:round(w/2)) = 1;
im = im - mean(im);
showft(im,nf, 1);

for j = 1:100
im = rand(w,1) > 0.95;
im = im - mean(im);
pwr(j,:) = abs(fft(im));
end
subplot(nf,2,3);
plot(im);
subplot(nf,2,4);
plot(fftshift(mean(pwr)));
%showft(im,nf, 2);
end

function showft(im, nf, n)
n = n-1;
subplot(nf,2,n*2+1);
plot(im,'o-');
subplot(nf,2,(n*2)+2);
plot(fftshift(abs(fft(im))));
