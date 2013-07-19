function fig = addhist(file)

x = dlmread('test.em',' ');
hold on
hist(x);
hold off
