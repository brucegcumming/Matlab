function ResampleTest(nreps)
%
%  Effect of sample size on variability of estimates produced
% by bootstrap, compared with Guassian SEM measues.

if nargin == 0
    nreps = 1000;  %number of runs with new data
end
nresample = 100; %number of resample runs for each data set. 
%Since we only calculate SD of the resample means, 100 is ample
%Efron and Tibshirani  p. 52


samplesz = [5:20 25:5:100];

for j = 1:length(samplesz)
    S = randn(nreps,samplesz(j));  %generate nreps observation sets
%for each set of observations at this samples size, look at the SD of the
%mean of the redamples. Compare this with SEM calculated with the usualr 
%Gaussian assumptions. 
    for k = 1:size(S,1)
        x = S(k,:);
        newindices = unidrnd(samplesz(j),samplesz(j),nresample);
        resample = x(newindices);
        resamplesd(k) = std(mean(resample)); %Obversed SD of the means of the resamples
    end
    resamplesem(j) = mean(resamplesd); %mean value of SEM from resamples
    parametricsem(j) = mean(std(S))./sqrt(samplesz(j)); %SEM calculated from the same samples
    truesem(j) = 1./sqrt(samplesz(j));
end

plot(samplesz,resamplesem./truesem);
hold on;
plot(samplesz,parametricsem./truesem,'r');
ylabel('SD of sample mean/theoretical SEM');
xlabel('Sample size');
legend('bootstrap','parametric');

