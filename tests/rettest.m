function [vector, matrix, spk] = rettest(Trials)

vector = [Trials.Start];
spk = [Trials.Count];

for j = 1:length(Trials)
  matrix(:,j) = ones(20000,1);
end


