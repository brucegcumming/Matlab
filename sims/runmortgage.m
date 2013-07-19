if ~exist('taxrate','var');
    taxrate = 0.25;
end
if exist('figa','var');
    figure(figa);
else
    figa = figure;
end
hold off;
j = 1;
rates = 0.0:0.005:0.1;
plotrates = [0 0.05 0.1];

interest = 5.625;
cost = 400000;

for growthrate = rates
capsum(1,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','income','capital',cost,'borrow',cost);
capsum(2,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','income','capital',0,'borrow',0);
if ismember(growthrate,plotrates)
    plot(capsum(1,:),'g');
    hold on;
    plot(capsum(2,:),'r');
    text(30,mean(capsum(:,30)),sprintf('%.2f%%',growthrate*100));
end
gain(j) = capsum(1,31) - capsum(2,31);
j = j+1;
end
labels{1} = sprintf('30yr, %.3f%% %.0f%% Tax on each years growth',interest,taxrate * 100);
title(labels{1});
if exist('figb','var');
    figure(figb);
else
    figb = figure;
end
hold off;
plot(rates,gain,'b');
xlabel('investment growth');
ylabel('Gain from mortgage');

%Compare 15 year and 30 year mortgate
%
%
%
figure(figb+1);
hold off;
j = 1;
[a,b,income] = amortize(cost, 4.95, 15);
for growthrate = rates
%set taxrate ot 0 to remove tax break on mortgage interest.
%gain(j) = mortgage(growthrate,'taxrate',0,'taxmodel','none');
capsum(1,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','capital','capital',0,'borrow',cost,'term',30,'interest',5.65,'income',income);
capsum(2,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','capital','capital',0,'borrow',cost,'term',15,'interest',4.95,'income',income);
if ismember(growthrate,plotrates)
    plot(capsum(1,:),'g');
    hold on;
    plot(capsum(2,:),'r');
    text(30,mean(capsum(:,30)),sprintf('%.2f%%',growthrate*100));
end
gain(j) = capsum(1,31) - capsum(2,31);
j = j+1;
end
capsum(1,:) = mortgage(0.05,'taxrate',taxrate,'taxmodel','capital','capital',0,'borrow',cost,'term',30,'interest',5.65,'income',income,'noise',1);
plot(capsum(1,:),'b');

labels{2} = sprintf('%.0f%% Captial Gains. 30 yr 5.625%% vs 15 yr 4.95%%',taxrate*100);
title(labels{2});
legend(sprintf('30 yr, 5.625%%'), sprintf('15yr, 4.95%%'),2);
figure(figb);
hold on;
plot(rates,gain,'r');


%Compare 30 years with mortgage/investment, assuming investment
%growth is taxed only at the end, as Capital gains.
%
%
figure(figb+2);
hold off;
j = 1;
for growthrate = rates
capsum(1,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','capital','capital',cost,'borrow',cost);
capsum(2,:) = mortgage(growthrate,'taxrate',taxrate,'taxmodel','capital','capital',0,'borrow',0);
if ismember(growthrate,plotrates)
    plot(capsum(1,:),'g');
    hold on;
    plot(capsum(2,:),'r');
    text(30,mean(capsum(:,30)),sprintf('%.2f%%',growthrate*100));
end
gain(j) = capsum(1,31) - capsum(2,31);
j = j+1;
end
labels{3} = sprintf('%.0f%% Tax on Capital Gain',taxrate*100);
title(labels{3});
figure(figb);
hold on;
plot(rates,gain,'g');


legend(labels);
plot([rates(1) rates(end)],[0  0],'k:');

