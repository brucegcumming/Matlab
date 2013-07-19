function MarkSlope(x,y,slope)


holdstat = ishold;
hold on;
plot([min(x) max(x)],[mean(y) - (mean(x)-min(x)) * slope mean(y)-(mean(x)-max(x)) * slope]);
if ~holdstat
    hold off;
end