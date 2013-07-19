function simplesimple(km,kb,varargin)
%simplesimple(km,kb,...)  simplified simulation of a binocular simple cells
%disparity tuning. (Actaully a simple cell pair - squaring not half.
%Contrast gain control is appliied to monocular responses, before summation, according to
%constant km, then to the summed response according to the value of kb,
%then squared (like Truchard et al 2000).
% The relationship between gain (G) and contrast (C) is defined so that
% gain is 1 when contrast is one.
%
% G = (C+k)/(C * (1+k))
%
% for binocular normalization, C refers to the mean of the contrasts in the
% two eyes, k is kb
%
%simplesimple(km,kb,'power') uses a power law to set gain, where
% G = C ^k
% for comparison with Truchard et at.
% simplesimple(-0.45,-0.5,'power', gives a nice replication of typical DPIC
% data, with the (low,low) response at PREF slightly greater than
% (low,high) at pref.
%
%Dashed lines show responses with NO binocular normalization for
%comparison. N.B. these are lower overall since the total gain (km * kb) is
%lower. Should probably adjust s that gl*gb + gr*gb is a constant for the
%comparison.

step = pi/20;
x = -pi:step:pi-step;
sl = sin(x);
colors = mycolors;
hold off;
np = 1;
cl = 1;

powerlaw = 0;
j = 1;
while (j < nargin-1)
    if(strncmpi(varargin{j},'power',3))
        powerlaw = 1;
    end
    j = j+1;
end

for cl = [0.1 0.95]; %left contrast
for cr = [0.1 0.95]; %right contrast
j = 1;
cb = (cl+cr)/2;
gb = (cb+kb) / (cb * (1+kb));
gl = (cl+km) / (cl * (1+km));
gr = (cr+km) / (cr * (1+km));

if powerlaw
    gb = cb^kb;
    gl = cl^km;
    gr = cr^km;
end

%
%response is a sinusoidal function of phase in each eye. Evaluate (l+r)^2
%over a range of phases and sum. Do this sum for each interocular phase
%difference.
%

for dp = -pi:step:pi
  sl = sin(x) .* cl .* gl;
  sr = (sin(x+dp) .* cr) .* gr;
  b(j) = mean((gb .*(sr + sl)).^2);
  sr = (sin(x+dp) .* cr) .* gr;
  sl = sin(x) .* cl .* gl;
  c(j) = mean(((sr + sl)).^2);
  dps(j) = dp;
  j = j+1;
end

gains(np,1) = gl;
gains(np,2) = gr;
gains(np,3) = gb;
plot(dps,b,'color',colors{np});
hold on;

h(np) = plot(dps,c,':','color',colors{np});
labels{np} = sprintf('co L%.2f R%.2f gain B%.2f R%.2f',cl,cr,gb,gr);
np = np+1;
end
end
legend(h,labels);
