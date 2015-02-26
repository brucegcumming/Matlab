function [cp, details] = CalcCP(acounts, bcounts, varargin)

%function [cp, details] = CalcCP(acounts, bcounts, varargin)
%Calcs a choice prob between two vectors of spikecounts. 
%CP > 0.5 when acounts tends to be > bcounts
%
%needs row vectors
nresample = 0;
j = 1;
while j < nargin-1
    if strncmpi('resample',varargin{j},5)
        j = j+1;
        nresample = varargin{j};
    end
    j = j+1;
end

n = length(acounts);
m = length(bcounts);
if m * n == 0
    cp = NaN;
    details.cp = NaN;
    details.n = n+m;
    details.pval= NaN;
    return;
end
%this step needs row vectors
if diff(size(acounts)) < 0
    zcounts = [acounts' bcounts'];
else
    zcounts = [acounts bcounts];
end

detect = [ones(1,n) zeros(1,length(zcounts)-n)];
fpos = ~detect;

[a, idx] = sort(zcounts);
d = 1 - cumsum(detect(idx))./n;
c = 1 -cumsum(fpos(idx))./m;
% Be careful of ties. By using sort, the data are treated as if you can put
% a criterion between two identical counts. If these have different choices
% associated, the fpos/detect rates depend on the order of the two
% identical counts. So only use these arrays where the counts have changed.
steps = find(diff(zcounts(idx)) > 0);
%reverse steps - if trapz is given positive monotonically reducing values for X and
%Y, it returns a negative number!!

getcrit = 1;
if getcrit
    [minx, id] = min(abs(d-(1-c)));
    details.critperf = mean([d(id) 1-c(id)]);
end
rix = length(steps):-1:1;
steps = steps(rix);
d = [0 d(steps) 1];
c = [0 c(steps) 1];
nullcp = trapz(d,c);
cp = trapz(c,d);
cpdetails.cp = cp;
details.detect = d;
details.fpos = c;
details.n = [n m];
details.sorted = a;

if nresample & n > 1 & m > 1
    for j = 1:nresample
        rd = detect(randperm(n+m));
        rf = ~rd;
        d = 1- cumsum(rd)./n;
        c = 1- cumsum(rf)./m;
        d = [0 d(steps) 1];
        c = [0 c(steps) 1];
        pcp(j) = trapz(d,c);
    end
    details.resamp = pcp;
    details.pval = 2 * min([length(find(pcp > cp)) length(find(pcp < cp))])/(nresample);
else
    details.pval = NaN;
end
