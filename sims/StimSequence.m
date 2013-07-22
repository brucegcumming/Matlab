% This simulation is to see whether the stimulus in a given trial is less 
% likely to be the same sign as the previous stimulus when stimuli are randomized
% within blocks with the constraint that all stimuli must be presented the
% same number of time within the block.
function res = sim(niter)

if nargin == 0
    niter = 10000;
end

n = [4,6,10,20,40,100]; % number of stimuli
m = [2,5,10,20,40,100]; % number of repetitions

res = [];
for nidx=1:length(n)
    for midx=1:length(m)
        n_stim = n(nidx);
        n_rep = m(midx);
        p = NaN(1,niter);
        for i=1:niter
            % generate stimuli
            s = repmat(sort([-1*[1:(n_stim/2)],1:(n_stim/2)]),n_rep,1);
            s = s(randperm(n_stim*n_rep));
            shift_prd = s(1:n_stim*n_rep-1).*s(2:n_stim*n_rep);
            p(i) = length(find(shift_prd<0))/(n_stim*n_rep-1);
        end
        res = [res;...
            n_stim, n_rep, mean(p), std(p)];
    end
end

clf
c = colormap(jet(length(n)));
for i=1:length(n)
    str{i} = sprintf('n stim = %d',n(i));
    x = res(res(:,1)==n(i),2);
    y = res(res(:,1)==n(i),3);
    line(x,y,'Marker','o','Color',c(i,:))
end
legend(str)
xlabel('Number of Repetition')
ylabel('P(sign(trial i) ~= sign(trial i-1))')
set(gca,'TickDir','out','XScale','log','Box','off')
return