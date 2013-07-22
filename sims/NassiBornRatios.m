function NassiBornRatios(rate,ratio)
%effect of noise in calulation of figure 7.
%rate is the firing rate for large diameter stimulu
%for this simulation assumed to be the same underlying rate
%as for the control.   At optimum size response = rate * ratio
if nargin == 0
    rate = 5;
    ratio = 4;
end

peak = poissrnd(rate*ratio,1,100); % = CONTROL_C
cool = poissrnd(rate,1,100);  %=COOL_L 
control = poissrnd(rate,1,100); %=CONTROL_L

ss = (peak-control)./peak;
mi = (cool-control)./(cool+control);

plot(ss, mi,'o');
xlabel('SS%control');
ylabel('Effect magnitude index');
%
%This ratio does not show the bias
%plot(control./peak,cool./peak,'o');

