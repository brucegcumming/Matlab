function PlotSdfs(sdf, lb)


colors = 'rgbykcrgbykc';
for j = 1:length(sdf.n);
  plot(sdf.times,sdf.sdf(j,:),colors(j));
  hold on;
  labels{j} = sprintf('%s = %.2f (%d)',lb,sdf.val(j),sdf.n(j));
end
legend(labels);
