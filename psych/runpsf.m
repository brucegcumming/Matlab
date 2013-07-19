data = readpsychsum('fle288.uadt');

fit = fitpsf(data,'twomean');
idx = find(data.expno == 1)
plotpsych(fit.data(idx),fit.fit(1),fit.fit(2));
hold on;
idx = find(data.expno == 2)
plotpsych(fit.data(idx),fit.fit(3),fit.fit(2));

