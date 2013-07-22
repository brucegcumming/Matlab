colors = ['b','g','k','m','c','y','r','b','g','k'];
if(~exist('len'))
  len = 1000;
end

clear val;
hold off;
j = 1;
val(1) = 0;
clear ft;
clear gft;
for pcorr = [0 0.03 0.06 0.1]
for loop = 1:2
  for k = 2:len;
    pval = 0.5 - val(k-1) * pcorr;
    val(k) = val(k-1) + binornd(1,pval,1,1)-0.5;
  end
  ft(loop,:) = abs(fft(val - mean(val)));
  gpath(1) = 0;
  for k = 2:len;
    gpath(k)  = gpath(k-1) + randn(1,1) - gpath(k-1) * pcorr;
  end;
  gft(loop,:) = abs(fft(gpath - mean(gpath)));
end
mft = mean(ft) / trapz(mean(ft));
plot(mft(1:100),'Color',colors(j));
hold on;
  labels{j} = sprintf('%d Binomial Corr %.3f steps',len,pcorr);
  j = j+1;

if(showgauss)
  mft = mean(gft) / trapz(mean(gft));
  line  = plot(mft(1:100),'Color',colors(j));
  labels{j} = sprintf('%d Gaussian steps',len);
  j = j+1;
end;
end;
legend(labels);