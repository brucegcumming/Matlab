function sparseness(counts)

sp = (1 - mean(counts)^2/mean(counts.^2))/(1 - 1/length(counts))
hist(counts);