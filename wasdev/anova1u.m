function p = anova1u(disptuning)
% ANOVA test whether it is tuned to disparity
% See help for anovan. Although I am only doing one-way anova, I have to use anovan because we may have unequal reps
X = []; group{1} = [];
nd = length([disptuning.x]);
for jj=1:nd
    X = [ X sqrt(disptuning.counts{jj}) ];
    group{1} = [ group{1} [disptuning.x(jj)].*ones(1,[disptuning.n(jj)])];
end

[p,table] = anovan(X,group,[],[],[],'off');

