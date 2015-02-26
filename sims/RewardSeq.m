function allcounts = Rewardseq(p, varargin)
%counts = Rewardseq(p)
%Calculate frequency distribution of consecutive correct triasl
%give a p(correct) of p
%generates a sequence of binary values with p, and also adds to
% the sequence a fraction of trials where p = 0.5 (modelling 0 dispatity
% trials) this fraction is 0.1 by default
% Rewardseq(p, 'pzero', x) set this fraction to x
%
% By default, runs 100 blocks each with 100 trials
% Rewardseq(p, 'trials', T, 'blocks' B) runs B blocks of T trials
% if p is a vector, runs the simulation for each p[ value and retunes
% a matrix with one row of counts for each p.
% So, to look at how the distribution looks if its made up by summing to p
% values:
% x = Rewardseq([0.9 0.7]); bar(sum(x));
%
% Also simulates simple reward modulated behaviour
% x = Rewardseq(p, 'rewardcause', [T q]);
% build a sequence where p(ccrrect) = p when onsecutive correct < T.
% Once consecutive correct reaches T, then p(correct) = q;
% N.B. to compare the two possiblitles need to adjst p,q to produce th%e
% overall miss rate in the comparison  E.G:
%x = RewardSeq(0.738];
%y = RewardSeq(0.7,'rewardcause',[3 0.1]); both produce miss rates of 28.5%
%
% RewardSeq('compare')  runs one pair of simulations and compares the
% distributions. 
pzero = 0.1; %fraction of trials rewarded at random
ntrials = 100;
nblocks = 100;
rewardcause = 0;
rewardtrial = 3;

if strcmp(p,'compare')
    x = RewardSeq([0.738 ],'blocks',10000);
    y = RewardSeq([0.7 ],'blocks',10000,'rewardcause',[3 0.1]); 
    if length(y) > length(x)
        x(length(y)) = 0;
    end
    plot(0:length(x)-1,y-x);
    xlabel('Number of consecutive correct responses');
    ylabel('N (Causal) - N(Indep)');
    return;
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'blocks',5)
        j = j+1;
        nblocks = varargin{j};
    elseif strncmpi(varargin{j},'pzero',5)
        j = j+1;
        pzero = varargin{j};
    elseif strncmpi(varargin{j},'reward',5)
        j = j+1;
        if length(varargin{j}) == 2
            rewardtrial = varargin{j}(1);
            rewardcause = varargin{j}(2);
        else
        rewardcause = varargin{j};
        end
    elseif strncmpi(varargin{j},'trials',5)
        j = j+1;
        trials = varargin{j};
    end
    j = j+1;
end

q = p + rewardcause; %p(correct) when reward is high
for k = 1:length(p)
    nzero = round(ntrials * pzero); %number of trials with random reward
    nc  =0;
    counts = zeros(1,20);
    nt = 0;
    for b = 1:nblocks
        y = rand(1,ntrials-nzero) < q(k);
        x = rand(1,ntrials-nzero) < p(k);
        z = rand(1,nzero) < 0.5; 
        x = [x z]; %add random trials and p(correct) = p trials
        x = x(randperm(length(x))); %shuffle the trials in the block
        y = [y z];
        y = y(randperm(length(y)));
        for j = 1:length(x)
            %set p according to nc, if rewardcause has been set
            if nc >= rewardtrial && rewardcause > 0
                resp = y(j);
            else
                resp = x(j);
            end
            if resp == 0 %now count how many were correct
                if nc > length(counts) -1
                    counts(nc+1) = 0;
                end
                counts(nc+1) = counts(nc+1)+1;
                nc = 0;
            else
                if nc == 0 % prev trial was a miss
                    counts(1) = counts(1)+1;
                end
                nc = nc+1;
            end
            n = counts(1) + sum(counts(2:end).* [1:length(counts)-1]);
            if n > j + (b-1) * ntrials
                n = counts(1) + sum(counts(2:end).* [1:length(counts)-1]);
            end
        end
        n = counts(1) + sum(counts(2:end).* [1:length(counts)-1]);
        if nc > length(counts) -1
            counts(nc+1) = 0;
        end
        nt = nt+length(x);
        score(b) = mean(x); %actual %correct
    end
    if(nc > 0) %last trial was good
        counts(nc) = counts(nc)+1;
    end
    n = counts(1) + sum(counts(2:end).* [1:length(counts)-1]);
    bar(counts(2:end));
    n = counts(1) + sum(counts(2:end).* [1:length(counts)-1]);
    fprintf('%.2f%% misses %d/%d\n',100* counts(1)./nt,counts(1),nt);
    allcounts(k,1:length(counts)) = counts;
end
if length(p) > 1
    plot(allcounts');
end
