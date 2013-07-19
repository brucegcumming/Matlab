%
% compare t-test and ANOVA false positive rates for data like 
% Doung and Freeman.
%

meanrate = 10;
for j = 1:1000
resps = meanrate  + randn(5,12) .* sqrt(meanrate);  %12 phase differences, 5 repetitions
p(j) = anova1(resps,[],'off');

%sresp = 1000 .* sin([1:12].*pi/(6)+rand(1,1).*pi.*2);
F = fft(mean(resps));
phase(j) = -angle(F(2)); %peak for 1st stim in list = 0 phase
maxi = mod(round(12 .* phase(j)/(2.*pi)),12)+1; 
mini = mod(maxi+6,12)+1;
[h, t(j), a,b] = ttest2(resps(:,maxi),resps(:, mini));

%hold off; plot(mean(resps)); hold on; plot([maxi maxi], [0 10]);
%pause;
maxs(j) = maxi;
end
