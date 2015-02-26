function [dips, errs] = GetDips(C)
errs = [];

if isfield(C,'fitdprime')
    dips(1) = C.fitdprime(1);
end
if isfield(C,'mahal')
    dips(2) = C.mahal(4); %1D GM fit
    dips(3) = C.mahal(1); %2D Gm fit
end
if isfield(C,'mydip')
    dips(4) = C.mydip(3);
else
    dips(4) = NaN;
    errs.s = sprintf('Missing Mydip in C%d\n',C.probe(1));
    errs.p = C.probe(1);
    if nargout == 1
    fprintf('%s\n',errs.s);
    end
end
if isfield(C,'hdip') %hartigans
    dips(5) = C.hdip;
else
    dips(5) = NaN;
end

