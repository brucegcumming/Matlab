function [allresps, details] = runbas(varargin)


dornd = 0;
dogauss = 0;
plotstim = 0;
simgauss = 1;
allresps = [];
sfs = [0.05  0.8 1.28 0.2];
tfs = [2 2.8 4 5.7 8]
speeds = [1:5:100]; %need to go up to 100 deg/sec
  %speed is pixels(0.1 deg) frames
  %10ms per frames. 1pixel/frame = 10 degress/sec
hold off;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'speeds',4)
        j = j+1;
        speeds = varargin{j};
    end
    j = j+1;
end

if dornd
pwrspec = ftstim('nrep',1000);
CalcResp(fftshift(pwrspec),1);
end

if dogauss
x = -127:128;
[X,Y,Z] = gauss2d([2 4],x);
for j = 1:length(speeds)
[resp, angles] = CalcResp(Z,speeds(j));
id = find(angles > -3.*pi/4 & angles <pi/4);
[a,b] = max(resp(id));
prefo(j) = angles(id(b));
end
plot(speeds,prefo);
end



tsd = 0.16./tfs;
if simgauss
    sy=0.5;
    sx=0.5;
x = -127:128;
    oris = [0:pi/8:pi/2];
%    build a space-time stimulus
for tfi = 1:length(tfs)
    tf = tfs(tfi);
% time, use 10ms per pixel.
    tk = Gabor([0.1.*tsd(tfi) tsd(tfi) pi/2 1 0 0],'pix',0.01); %Odd symmetric Gabor for T
    tk = tk./sum(abs(tk)); % normalize ;
%+-128 time units =2.56 sec = 0.39 Hz fundamental
% SD of 0.08 is harmonic 5, = 2Hz
    tks(tfi,:) = tk;
    ts = [1:length(tk)] - length(tk)/2;
    tpwr(tfi,:) = abs(fft(tk));
    details.tftfrq = [1:length(tk)] ./(256 * 0.01);
    details.tk(tfi,:) = tk;
end
plot(tpwr);

for tfi = 1:length(tfs)
    tk = tks(tfi,:);
    for sfi = 1:length(sfs)
    sf = sfs(sfi);
    sy=0.5./(2*sf);
    sx=sy;
for j = length(oris):-1:1
   rfs{j} = Gabor([sf sx 0 1 0 0 oris(j) sy],'pix2deg',0.1);
end

    for si = 1:length(speeds)
        speed = speeds(si);

        for t = 1:length(ts)
 %SD of 50,100  pixels = 5, 10 degrees 
    [X,Y,Z] = gauss2d([50 100],x ,'mean',[0 ts(t).*speed]);
    if plotstim
    imagesc(Z); drawnow;
    end
    for j = length(oris):-1:1
        resps(j,t) = sum(rfs{j}(:) .* Z(:));
    end
end
figure(1)
subplot(2,1,1);
imagesc(resps.^2);
subplot(2,1,2);
for j = length(oris):-1:1
    cresps(j,:) = conv(resps(j,:),tk);
end
imagesc(cresps.^2);

figure(2);
allresps(si,sfi,tfi,:) = sum(cresps.^2,2);
plot(squeeze(allresps(si,sfi,tfi,:)));
title(sprintf('Speed %.2f',speed));
    end
end
end
end
details.tfpwr = tpwr;
