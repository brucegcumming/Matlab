function [oresp, dirs] = CalcResp(pwr, speed, varargin)

pixperdeg = 4;
usewov = 0;
%given a 2-D power spectrum for a simulus, calculate summed response along
%each orientation, for a given speed
% M+C use Sf 0.05 to 0.2 (4, log spaced) cycles per degree 2 - 8 Hz
% and sd = 1/3 of peak

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sfonly',3)
        usewov = 1;
    elseif strncmpi(varargin{j},'tfonly',3)
        usewov = 2;
    end
    j =j+1;
end
x = [1:size(pwr,1)]- floor(size(pwr,1)/2);
y = [1:size(pwr,2)]- floor(size(pwr,2)/2);
[X,Y] = meshgrid(x,y);

%pwrspec is just powers, so Nyquist limit is 0.5 cycle pixel
angles = atan2(Y,X);
sfs = abs(X+i.*Y); %cycles per image
sfs = sfs .* pixperdeg./length(pwr);
%speed = tf/sf;
tf = speed .* cos(angles).* sfs;

twov = SetWov(tf, [ 2 4 6 8 ], 0.33);
swov = SetWov(sfs, [ 0.05 0.1 .15 0.2], 0.33);

if usewov == 1
resps = pwr .* swov;
else
resps = pwr .* twov .* swov;
end
dirs = -pi:pi/20:pi-0.1;
for n = 1:length(dirs)
    d = dirs(n);
    w = cos(angles-d).^4;
    id = find(abs(angles-d) < pi/40 & sfs < 120);
    oresp(n) = mean(resps(id));
    oresp(n) = sum(resps(:) .*w(:));
end
subplot(2,1,1);
plot(dirs,oresp);
subplot(2,1,2);
imagesc(resps);

%imagesc(tf);

function wov = SetWov(f, c, k)

wov = zeros(size(f));
for j = 1:length(c)
    sd = k .* c(j);
    wov = wov + exp(-(abs(f)-c(j)).^2./(2 .* sd.^2));
end