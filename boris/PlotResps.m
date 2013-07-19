function PlotResps(stripresps)

stripresps(1,:) = stripresps(1,:)./max(stripresps(1,:));
stripresps(2,:) = stripresps(2,:)./max(stripresps(2,:));


for a = 0.05:0.01:0.5
    plot(stripresps(1,:).^a,'b');
    hold on;
    plot(stripresps(2,:).^a,'r');
end