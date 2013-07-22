j = 1;
rates = 0.05:0.01:0.1;
for growthrate = rates
gain(j) = mortgage(growthrate);
j = j+1;
end

figure;
plot(rates,gain);