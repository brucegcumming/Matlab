function xc = CompareAllExpt(a,b, varargin)

na = length(a.Header);
nb = length(b.Header);

n = min([na nb]);
for j = 1:n
    c = a.Header(j).cellnumber;
    Ea = All2Expt(a,c);
    Eb = All2Expt(b,c);
    if length(Ea.Trials) == length(Eb.Trials)
        x = corrcoef([Ea.Trials.count],[Eb.Trials.count]);
        overlap(j) = 1;
    else
        [x, ta, tb] = intersect([Ea.Trials.id],[Eb.Trials.id]);
        x = corrcoef([Ea.Trials(ta).count],[Eb.Trials(tb).count]);
        overlap(j) = length(ta)./max([length(Ea.Trials) length(Eb.Trials)]);
    end
    xc(j) = x(1,2);
end
min(xc);