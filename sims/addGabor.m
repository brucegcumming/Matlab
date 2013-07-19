%takes two gabors, raises to different powers, and adds
%

if ~exist('gmode')
    gmode = 1;
end
if gmode ==1 %%make a signle hump
    x = Gabor([0.5 0.8]) + 1;
    y = 1 - (x-1);
    z = x.^2 + y*1.6;
    az = x*1.6 + (y.^2);

    hold off;
    plot(x);
    hold on;
    plot(y,'r');
    plot(z,'k');
    plot(az,'k:');
elseif gmode ==2 %% try to make corr odd and AC even.
    x = Gabor([0.25 0.8 2.5 1 0]);  %%f sd phase A pos
    ax = -x;
    y = Gabor([0.4 0.8 1.57 1 0.25]);
    ay = -y;
    gains = [1 1];
    powers = [2.5 2.5];
    subplot(2,1,1);
    hold off;
    plot(x+1);
    hold on;
    xp = (x+1).^powers(1);
    yp = (y+1).^powers(2);
    axp = (ax+1).^powers(1);
    ayp = (ay+1).^powers(2);
    plot(gains(1).* xp,':')
    plot(y+1,'r');
    plot(gains(2).* yp,'r:');
    csum = (gains(2) * yp+ gains(1) * xp)/sum(gains);
    plot(csum,'k','linewidth',2);
    subplot(2,1,2);
    hold off;
    asum = (gains(2) * ayp+ gains(1) * axp)/sum(gains);
    plot(gains(2).* ayp,'r:');
    hold on;
    plot(gains(1).* axp,':')
    plot(asum,'k','linewidth',2);
end