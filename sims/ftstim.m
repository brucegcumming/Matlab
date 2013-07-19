function sumpwr = ftstim(varargin)
imw = 256;
barh = 15;
barw = 3;
%16 pixels x 4 pixels  = 4x1 degree, seems to be normal.
sumbars = 0;
nreps = 100;
%aim for 
nbars = 100;

bp = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'bar',3)
        j = j+1;
        barh = varargin{j}(1);
        if length(varargin{j}) > 1
            barw = varargin{j}(2);
        end
    elseif strncmpi(varargin{j},'sum',3)
        sumbars = 1;
    elseif strncmpi(varargin{j},'nbars',3)
        j = j+1;
        nbars = varargin{j}(1);
        if nbars == 1 & length(varargin{j}) == 3
            bp = varargin{j}(2:3);
        end
    elseif strncmpi(varargin{j},'nrep',3)
        nreps = varargin{j};
    end
    j = j+1;
end
         
w = barh;
im = zeros(imw+w,imw+w);
sumpwr = zeros(imw,imw);
pos = ceil(rand(nbars,2) .* imw);
if nbars == 1 && length(bp)
    pos = bp;
end

for n = nreps:-1:1
im = zeros(imw+w,imw+w);
if sumbars
for j = 1:nbars
    im(pos(j,1):pos(j,1)+barw,pos(j,2):pos(j,2)+barh) = im(pos(j,1):pos(j,1)+barw,pos(j,2):pos(j,2)+barh)+1;
end
else
for j = 1:nbars
    im(pos(j,1):pos(j,1)+barw,pos(j,2):pos(j,2)+barh) = 1;
end
end
im = im - mean(im(:));
im = im(1:imw,1:imw);
pwr = abs(fft2(im));
sumpwr = sumpwr + pwr;

end
subplot(2,1,1); 
imagesc(im);
subplot(2,1,2); 
imagesc(fftshift(sumpwr));
%imagesc((sumpwr));
drawnow;
