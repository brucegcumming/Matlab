function r = sem(x)
%r = sem(x) std(x)./sqrt(length(x))
r = std(x)./sqrt(length(x));
