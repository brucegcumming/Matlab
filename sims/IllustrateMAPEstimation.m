function IllustrateMAPEstimation(varargin)

close all
clc

% Constants
alpha      = 4;                                                            % Alpha = 4 is very close to median V1; lower values give stronger excitability fluctuations
conditions = 15;                                                           % Number of stimulus conditions
repeats    = 50;                                                           % Number of repeats per condition
runs       = 1;                                                            % Number of simulation runs (1 to get example; 1,000 to estimate mean correlation)
stimResp   = logspace(0, 1.5, conditions);                                 % In spikes per second
stimDur    = 1;                                                            % In seconds
type       = 1;                                                            % p(E) follows (1) gamma, (2) uniform, (3) log-Normal distribution, or (4) rectified Normal

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'type',4)
        j = j+1;
        type = varargin{j};
    end
    j = j+1;
end
% Initial computations
v    = logspace(0, 1.5, 256);
VMR  = v + (1/alpha)*v.^2;
varE = 1/alpha;

% Simulation loop
for i = 1:runs
    
    if round(i/25) == i/25;
        fprintf('\n \t Analysis completed for %d of %d conditions...', i, runs);
    end
    
    % Computations
    % I. Generate spike counts
    if type == 1                                                           % Gamma distributed
        excitability = gamrnd(alpha, 1./alpha, [conditions, repeats]);     
        excPdf       = gampdf(v, alpha, 1./alpha);
    elseif type == 2
        spread       = sqrt(12*varE);                                      % Uniform distributed
        excitability = 1 + spread*(rand([conditions, repeats]) - 0.5);     
        excPdf       = zeros(size(v));
        excPdf((v > 1-spread/2) & (v < 1+spread/2)) = 1/spread;
    elseif type == 3                                                       % lognormal distribution
        mu           = log(1/sqrt(varE + 1));
        sigma        = sqrt(log(varE/1 + 1));
        excitability = lognrnd(mu, sigma, [conditions, repeats]);          
        excPdf       = lognpdf(v, mu, sigma);
    elseif type == 4                                                       % rectified normal distribution
        excitability = max(.01, 1 + sqrt(varE)*randn([conditions, repeats])); 
        excPdf       = normpdf(v, 1, sqrt(varE));
    end
    
    genRatePoiss = excitability.*repmat(stimResp', [1 repeats]);           % The effective trial-by-trial firing rate (stim-drive * excitability) 
    spikeCounts  = poissrnd(genRatePoiss*stimDur);
    
    trial.blockID    = repmat([1:conditions]', [1 repeats]);               % trial is a structure used to estimate Var[E]
    trial.duration   = stimDur + zeros(conditions, repeats);
    trial.spikeCount = spikeCounts;
    
    
    % II. Get Var[E] out
    % Estimate the function relating mean spike count to spike count variance
    options        = optimset('Display', 'off', 'Maxiter', 10000, 'MaxFuneval', 20000);
    blockIDs       = 1:conditions;
    blockMeanRates = mean(spikeCounts, 2)/stimDur;
    obFun          = @(param) GiveBof(param, blockIDs, blockMeanRates, trial.blockID(:), trial.duration(:), trial.spikeCount(:));
    startValue     = 100;
    [alphaEst, x]  = fminsearch(obFun, startValue, options);
    
    % Define prior on excitability
    excRange = logspace(-4, 3, 1028);                                      % Excitability range
    priorExc = gampdf(excRange, alphaEst, 1/alphaEst);                     % Prior on excitability
    rangeCE  = repmat(excRange, [conditions, 1]).*repmat(stimDur*blockMeanRates, [1, length(excRange)]);  % Appropriate range for prior on conditional expectation
    rangeMat = rangeCE(trial.blockID(:), :);
    
    % III. Get likelihood
    llh = poisspdf(repmat(trial.spikeCount(:), [1, length(excRange)]), rangeMat);
    
    % IV. Get posterior
    posterior = repmat(priorExc, [conditions*repeats, 1]).*llh;
    
    % V. Get MAP-estimate and mean estimate
    [x, ind] = max(posterior');
    MAPexc   = excRange(ind);
    meanEst  = sum(posterior.*repmat(excRange, [length(ind) 1]), 2)./sum(posterior, 2);       
    E = reshape(MAPexc,size(trial.spikeCount));
end

adjspikecounts = spikeCounts./excitability; %Underlying renewal variability from estimated gain

% VI. Plot outcome
figure
loglog(excitability(:), MAPexc, 'ko', 'markersize', 6, 'markerfacecolor', [.75 .75 .75])
hold on, box off, axis square
loglog([.1 10], [.1 10], 'k--')
axis([.1 10 .1 10])
xlabel('Actual excitability', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
ylabel('Map-estimate excitability', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);


figure
loglog(v, VMR, 'k-', 'linewidth', 2, 'color', [0 0 0])
hold on, box off, axis square
loglog(mean(spikeCounts, 2), var(spikeCounts'), 'wo', 'markersize', 12, 'markerfacecolor', [1 0 0])
loglog([.1 1000], [.1 1000] , 'k--', 'linewidth', 1)
axis([.1 1000 .1 1000])
loglog(mean(spikeCounts, 2), var(adjspikecounts'), 'wo', 'markersize', 12, 'markerfacecolor', [0 0 1])
xlabel('Mean (spikes)', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
ylabel('Variance (spikes^2)', 'fontName', 'Helvetica', 'fontAngle', 'oblique', 'fontSize', 16);
legend('ground truth', 'observed data', 'deduced renewal','location', 'NorthWest')

figure;
plot(spikeCounts',E','o');
xlabel('Spike Count');
ylabel('Excitability');
end


%% Functions used in the main script
%% GiveBof
function [NLL] = GiveBof(params, blockIDs, blockMeanRates, passBlockID, passDurations, passSpikeCounts)

alpha = params(1);                                                         % The alpha parameter of the gamma-distribution

% Calculate the predicted spike count distributions under the saturated model
[x, loc] = ismember(passBlockID, blockIDs);
meanRate = blockMeanRates(loc);
mu       = max(.001, passDurations.*meanRate);                             % The predicted mean spike count
var      = mu + ((mu.^2)/alpha);
r        = (mu.^2)./(var - mu);                                            % The parameters r and p of the negative binomial distribution
p        = r./(r + mu);

% Evaluate the likelihood of each spike count under the current parameterization of the model
Llh = nbinpdf(passSpikeCounts, r, p);                                      % The likelihood for each pass
NLL = sum(-log(Llh));                                                      % The negative log-likelihood of the whole data-set
end
%%

