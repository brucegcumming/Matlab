function [angles, fdist] = RobAngles(varargin)

angles = [];
nresample = 10000;
mode = 0;

j = 1;
while j <= nargin
    if ischar(varargin{j})
        if strncmpi('gauss',varargin{j},3)
            GaussSimulate([pi/4:pi/4:pi]);
            return;
        elseif strncmpi('fitgauss',varargin{j},3)
            mode =1;
        end
    end
    j = j+1;    
end


a = textread('new.txt');  %unwrapped, smoothed data
fdist = a(:,2);
sd = std(fdist);
if mode == 1
    % robs units are 370 to the circle
    xvals = [1:length(fdist)] * 2 * pi/370;
    fits{1} = FitGauss(xvals,fdist');
%    angles = f2dist(fdist,1000,'xvals',[1:length(fdist)]);
    fits{2} = FitDist(xvals,fdist',1); %%Skewed Gaussians
    fits{3} = FitDist(xvals,fdist',0); %%two Gaussians
    angles = fits;
    return;
end
adist = subsample(fdist,27) * 27;  %estimate unsmoothed data
angles = []; %convert to list of angles
for j = 1:length(adist)
    for k = 1:adist(j)
        angles = [angles j]; 
    end
end
sd = std(angles)

%anlges are in 1/2.7 of a cup, but unwrapped. Convert these to radians
%370 units is one lap, if the rotor is stationary
wangles = angles * 2 * pi/(370/27);
[a,b] = Rayleigh(wangles,'resample',nresample);
fprintf('Stationary rotor r = %.4f p = %.3f\n',a,b);
%for a typical rotor speed, a ball travelling 20 cups of the stator
%completes on lap of the rotor
wangles = angles * 2 * pi/(370/27);
[a,b] = Rayleigh(wangles,'resample',nresample);
fprintf('20 Cup lap r = %.4f p = %.3f\n',a,b);


function GaussSimulate(sd)

x = -4*pi:0.1:4*pi;
for k = 1:length(sd)
    y = Gauss(sd(k),x);
    adist = round (y .* 1000/sum(y));
    ix = 1;
    angles = [];
    for j = 1:length(adist)
        if(adist(j))
            angles(ix:ix+adist(j)) = x(j);
            ix=ix+adist(j);
        end
    end
    r(k) = Rayleigh(angles,0);
end
plot(sd,r);
xlabel('SD of wrapped Gaussian');
ylabel('vector length');