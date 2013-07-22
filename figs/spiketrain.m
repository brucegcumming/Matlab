function volts = spiketrain(th)

showplot = 0;

volts = rand(1000,1);
spk = find(volts > th);
volts(spk) = 10;
if showplot
    plot(volts);
end