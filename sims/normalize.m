function normalize(varargin)
% test normalization
% calcuate two sets of random pairs. No true difference, but different
% means. Look at vairance of normalized response for raw (divided my mean
% raw), normalized, and sqrt
n = 100000;


a = 10 + randn(n,2) .* sqrt(10); 
b = 100 + randn(n,2).* sqrt(100);

nb = b ./repmat(mean(b,2),1,2);
na = a ./repmat(mean(a,2),1,2);
diffs = (diff(nb,[],2) + diff(na,[],2));
%mean of normalized differeces. denominator is 1 + 1, just for consistency
diffs = (diff(nb,[],2) + diff(na,[],2))./(mean(nb,2)+mean(na,2));
%ratios is a mean ratio for each element in na, nb. 
ratios = (nb(1,:) + na(1,:))./(nb(2,:) + na(2,:));

ratio(1) = mean(ratios);
stds(1) = std(diffs);
rstd(1) = std(ratios);
%mean difference divided by mean resp
diffs = (diff(b,[],2) + diff(a,[],2))./(mean(b,2)+mean(a,2));
ratios = (b(1,:) + a(1,:))./(b(2,:) + a(2,:));
ratio(2) = mean(ratios);
stds(2) = std(diffs);
rstd(2) = std(ratios);
nb = sqrt(b);
na = sqrt(a);
%mean difference divided by mean resp after sqrt
diffs = (diff(nb,[],2) + diff(na,[],2))./(mean(nb,2)+mean(na,2));
ratios = (nb(1,:) + na(1,:))./(nb(2,:) + na(2,:));
rstd(3) = std(ratios)
stds(3) = std(diffs)
ratio(3) = mean(ratios)
