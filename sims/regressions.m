function regressions(type, varargin)
%type 0 
%


if type == 3

    for j = 1:1000
        x = 1:10;
        rnd = randn(size(x));
        y = x + rnd .*sqrt(x);
        fit = polyfit(x,y,1);
        slopes(j) = fit(1);
    end
    hist(slopes);
    fprintf('Mean %.3f\n',mean(slopes));
    return;
end

x = 1:10;
rnd = randn(size(x));
y = x + rnd;
if type == 1
    y = x;
end
fit = polyfit(x,y,1);
GetFigure('Normal Regression');
plot(x,y,'o');
yfit = fit(2)  + fit(1).*x;
hold on;
plot(x,yfit,'-');
title(sprintf('slope %.2f',fit(1)));
x = x + randn(size(x));
plot(x,y,'ro');
fit = polyfit(x,y,1);
yfit = fit(2)  + fit(1).*x;
plot(x,yfit,'r-');
%GetFigure('Noize on X')
