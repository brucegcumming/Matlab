function spkp = FitSpikeShp(spike, impulse)
%fit a simple descriptive spike shape such than when convolved with the 
%impulse function it produces the observed mean spike

x(1) = 8;
x(2) = 1;
x(3) = 0;
x(4) = 20;
x(5) = 5;
options = optimset('TolFun',1e-6);
spkp = fminsearch(@TrySpike,x,options, spike, impulse);
hold off; TrySpike(spkp, spike, impulse, 'plot');


function [err, insp] = TrySpike(x, spike, impulse, varargin)

insp = zeros(size(spike));
imp = insp;
imp(10) = 1;
%insp(round(x(1))) = x(2);
npts = round(x(4));
npost = round(x(5));
if npts > 1
    insp(round(x(1))+[1:npts]) = interp1([1 npts],[x(2) x(3)],[1:npts]);
else
    insp(round(x(1))+1) = x(2);
    insp(round(x(1))+2) = x(3);
end
    
if round(x(1))+npts+npost+1 > 46
    npost = 46-(round(x(1))+npts+1)
end
insp(round(x(1))+npts-1+[1:npost]) = interp1([1 npost],[x(3) 0],[1:npost]);
p = conv(insp,impulse);
sp = conv(spike,imp);
err = sum((p-sp).^2);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plot',3)
        plot(p);
        hold on;
        plot(sp,'r');
        plot(insp,'g');
    end
    j = j+1;
end

