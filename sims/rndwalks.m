colors = ['b','g','k','m','c','y','r','b','g','k'];
if(~exist('len'))
  len = 1000;
end

j = 1;
%for len = [100 500 1000 5000];
for pstep = [1 0.5 0.1 0.01]
  clear ft;
  clear gft;
  for loop = 1:50
    steps = binornd(1,0.5,1,len) - 0.5;
    ifstep = binornd(1,pstep,1,len);
    path = cumsum(steps .* ifstep);
    ft(loop,:) = abs(fft(path - mean(path)));
    gpath = cumsum(randn(1,len));
    gft(loop,:) = abs(fft(gpath - mean(gpath)));
  end
mft = mean(ft) / trapz(mean(ft));
plot(mft(1:100),'Color',colors(j));
hold on;
labels{j} = sprintf('%d fixed steps pstep %.2f',len,pstep);
if(showgauss)
mft = mean(gft) / trapz(mean(gft));
j = j+1;
line  = plot(mft(1:100),'Color',colors(j));
labels{j} = sprintf('%d Gaussian steps',len);
end;
j = j+1;
end
legend(labels);