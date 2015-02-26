function truet = ConvertTime(Expt, t)
%truet = ConvertTime(Expt, t) Converts timestamps into matlab datenums

t = double(t);
truet = t; %default
if isfield(Expt,'Header')
    if isfield(Expt.Header,'CreationDate')
        truet = Expt.Header.CreationDate + t ./(24 .* 60 .* 60 .* 10000);
    elseif isfield(Expt.Header,'CreationDate')
        truet = Expt.Header.builddate + t ./(24 .* 60 .* 60 .* 10000);
    end
end
