function getbug


x = rand(100,2);
figure;
subplot(1,1,1);
hold off;
plot(x(:,1),x(:,2));
hold on;
plot(x(:,2),x(:,1),'r');
refline(1);
legend({'test1' 'test2'});
set(gca,'xlim',[0 1],'ylim',[0 1]);
refline(1);
h = get(gcf,'Children')
for j = 1:length(h)
    p = get(h(j));
    if isfield(p,'LineWidth')
        set(h(j),'LineWidth',2);
    end
end