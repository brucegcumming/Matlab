nsubj = 1;

if ~exist('rrmeans','var')
    load('rrresample.mat');
end
if ~exist('rlmeans','var')
    load('rlresample.mat');
end
rufdiffs = rrmeans - rlmeans;
stddiffs(nsubj,:) = rufdiffs ./ std(rufdiffs);
rawdiffs(nsubj,:) = rufdiffs;

nsubj = nsubj+1'
if ~exist('bgct','var')
    load('bgcresample.mat');
end
bgcdiffs = bgct(:,1)-bgct(:,2);
stddiffs(nsubj,:) = bgcdiffs' ./ std(bgcdiffs);
rawdiffs(nsubj,:) = bgcdiffs';

if exist('jcrresample.mat','file') 
nsubj = nsubj+1
if ~exist('jcrt','var')
    load('jcrresample.mat');
end
jcrdiffs = jcrt(:,1)-jcrt(:,2);
stddiffs(nsubj,:) = jcrdiffs' ./ std(jcrdiffs);
rawdiffs(nsubj,:) = jcrdiffs';
end

if exist('drresample.mat','file') & exist('dlresample.mat','file')
nsubj = nsubj+1
if ~exist('drmeans','var')
    load('drresample.mat');
end
if ~exist('dlmeans','var')
    load('dlresample.mat');
end

dufdiffs = dlmeans - drmeans;
stddiffs(nsubj,:) = dufdiffs ./ std(dufdiffs);
rawdiffs(nsubj,:) = dufdiffs;
end
subplot(2,1,1);
hist(sum(stddiffs));
title('Mean Z-Transformed Delays');
subplot(2,1,2);
hist(sum(rawdiffs));
title(sprintf('Mean Raw Delays %d subjects',nsubj));
xlabel('Delay (ms)');
ylabel('N');

figure
subplot(2,2,1);
hist(rufdiffs);
title(sprintf('Rufus'));
subplot(2,2,2);
hist(dufdiffs);
title(sprintf('Dufus'));
subplot(2,2,3);
hist(bgcdiffs);
title(sprintf('BGC'));
subplot(2,2,4);
hist(rufdiffs);
title(sprintf('JCR'));
