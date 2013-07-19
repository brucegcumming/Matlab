function sdf = Spks2Sdf(spikes, nt, times, width, flag)

rates = [];
if(strmatch(flag,'exp'));
for t = [times]
  rates = [rates sum(exp((spikes(find(spikes <t)) -t)/width))];
end
  sdf = [rates .* 10000/(width * nt)];
else
for t = [times]
  rates = [rates sum(exp(-((spikes-t).^2)/(2 * width^2)))];
end
sdf = [rates .* 10000/(sqrt(2 * pi) * width * nt)];
end
