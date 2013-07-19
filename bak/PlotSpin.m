function PlotSpin(Expt, varargin)

latency = 500;
delay = 500;
width = 2;
duration = 0;
POLAR = 1;
plottype = 0;
j = 1;
while j < nargin
    if strncmpi(varargin{j},'width',2)
        j = j+1;
        width = varargin{j};
    elseif strncmpi(varargin{j},'delay',2)
        j = j+1;
        delay = varargin{j};
    elseif strncmpi(varargin{j},'polar',2)
        plottype = POLAR;
    end
    j = j+1;
end

duration = mean([Expt.Trials.End] - [Expt.Trials.Start]);
aspikes = zeros(180,1);
areps = zeros(180,1);
breps = zeros(180,1);
bspikes = zeros(180,1);
ospikes = [];
for j = 1:length(Expt.Trials)
    ispks = find(Expt.Trials(j).Spikes > delay & Expt.Trials(j).Spikes <= duration+delay);
    spikes = Expt.Trials(j).Spikes(ispks)';
    spinrate = Expt.Trials(j).Sr;
    or = Expt.Trials(j).or;
    bins =  1+ mod(round(3600 + or + spinrate * (spikes-delay)/(10000)),180);
    obins =  1+ mod(round(3600 + or + spinrate * (1:duration)/10000),180);
    if Expt.Trials(j).st == 0 %blank
        ospikes = [ospikes spikes];
    elseif Expt.Trials(j).Sr > 0
        aspikes(bins) = aspikes(bins)+1;
        areps(obins) = areps(obins)+1;
    else
        bspikes(bins) = bspikes(bins)+1;
        breps(obins) = breps(obins)+1;
    end
end
iblnk = find([Expt.Trials.st] ==0);
bsdf = trigsdf(Expt.Trials(iblnk),200,delay:100:duration+delay);
brate = 10000 * length(ospikes) ./(length(iblnk) * duration);
spinrate = mean(abs([Expt.Trials.Sr]));
%Each bin = 1 degree, = 1/Sr seconds. 
aspikes = spinrate .* aspikes ./ areps;
bspikes = spinrate .* bspikes ./ areps;
x = -width*3:width*3;
kl = exp(-(x.^2)/(2 * width^2));
kl = kl ./sum(kl);
rates = conv([aspikes; aspikes],kl);
rates = [rates(width+90:width+269); rates(width+90:width+269)];
angles = [-90:269]';
if plottype == POLAR
    [x,y] = pol2cart(angles*pi/180, rates);
    plot(x,y);
else
    plot(angles,rates);
end
hold on;
rates = conv([bspikes; bspikes],kl);
rates = [rates(width+90:width+269); rates(width+90:width+269)]; 
bt = angles(1) + range(angles) .* [1:length(bsdf)]./length(bsdf);

if plottype == POLAR
    [x,y] = pol2cart(angles*pi/180, rates);
    plot(x,y,'r');
    [x,y] = pol2cart(bt*pi/180, bsdf');
    plot(x,y,'g');
    [x,y] = pol2cart(angles*pi/180, ones(size(angles)) .* brate);
    plot(x,y,'g');
else
    plot(angles,rates,'r');
    plot(bt,bsdf,'g');
end
legend('Sr+','Sr-','blank');
title(sprintf('%s Sr %.0f',splitpath(Expt.Header.Name),mean(abs([Expt.Trials.Sr]))));
if plottype == POLAR
    axis image;
end